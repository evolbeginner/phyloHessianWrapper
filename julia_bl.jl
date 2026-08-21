#! /bin/env julia


######################################################################
const LIB = "lib-julia"
push!(LOAD_PATH, joinpath(@__DIR__, LIB))

include(joinpath(LIB, "readArg.jl"))
include(joinpath(LIB, "Q.jl"))
include(joinpath(LIB, "util.jl"))

######################################################################
const GET_BRANCH_ORDER_IN_APE = joinpath(@__DIR__, "additional_scripts", "get_branch_order_in_ape.R")

using Dates
using LinearAlgebra
BLAS.set_num_threads(1)

using StatsFuns: logsumexp
using DelimitedFiles
using Base.Threads

######################################################################
struct Nb
    node::Int
    tip::Int
    branch::Int
end

Base.@kwdef struct RunCtx
    nb::Nb
    nl::Int
    sub_model::String

    rs::Vector{Float64}
    props::Vector{Float64}
    Fs
    Qrs::Vector{Float64}
    freqs::Vector{Float64}

    q_pis
    q_pis_sites
    is_pmsf::Bool

    inv_info::Dict{Symbol, Any}
    pattern
    all_children_vec::Vector{Vector{Int}}
    descendants

    # runtime options (captured once, no globals in hot code)
    hessian_outfile::Union{Nothing, String}
    hessian_type::String
	hessian_infile::Union{Nothing, String}
	fd_scheme::Symbol
    cache_mode::Symbol
    is_compare::Bool
    transform_method
    treefile::String
    std_outfh::IO
    branchout_matrix::String
    bs_branchout_matrix::Union{Nothing, String}
end

######################################################################
# Utilities
######################################################################
@inline choose_q(ctx::RunCtx, i::Int) = ctx.is_pmsf ? ctx.q_pis_sites[i] : ctx.q_pis

function dict_children_to_vec(all_children::Dict, n_nodes::Int)
    out = [Int[] for _ in 1:n_nodes]
    for (k, v) in all_children
        out[Int(k)] = Int.(v)
    end
    return out
end

# +I helper: site is invariant if, after removing gap tips (all-ones columns),
# all remaining observed tips have the same state.
@inline function is_invariant_site_from_liks(
    liks_ori::AbstractMatrix{Float64},
    nb::Nb;
    atol::Float64 = 1e-12
)::Bool
    tips = @view liks_ori[:, 1:nb.tip]

    seen_state = 0
    has_nongap = false

    @inbounds for t in 1:nb.tip
        col = @view tips[:, t]

        # gap coding in this pipeline: all entries are 1.0
        is_gap = true
        for s in eachindex(col)
            if abs(col[s] - 1.0) > atol
                is_gap = false
                break
            end
        end
        if is_gap
            continue
        end

        # expected one-hot observed state
        st = 0
        for s in eachindex(col)
            if col[s] > 1.0 - atol
                st = s
                break
            end
        end
        if st == 0
            return false
        end

        has_nongap = true
        if seen_state == 0
            seen_state = st
        elseif st != seen_state
            return false
        end
    end

    return has_nongap
end

######################################################################
@inline function eigendecom_pt!(
    out::AbstractVector{Float64},
    tmp::AbstractVector{Float64},
    U::AbstractMatrix{Float64},
    Lambda::AbstractVector{Float64},
    U_inv::AbstractMatrix{Float64},
    v::AbstractVector{Float64},
    scale::Float64 = 1.0
)
    mul!(tmp, U_inv, v)
    @inbounds @simd for k in eachindex(Lambda, tmp)
        tmp[k] *= exp(Lambda[k] * scale)
    end
    mul!(out, U, tmp)
    return out
end

######################################################################
# Caches for one fixed branch-length vector.
# :diag cache stores exp(Lambda * b * r) per edge/state.
# :full cache stores full P(t)=U*diag(exp(Lambda*t))*U_inv per edge.
 # Cache each edge independently.  This is important during finite differences:
 # a perturbed branch must not invalidate P(t) for all of the other edges.
const PTDiagCache = Dict{Tuple{UInt64, Float64, Int, Float64}, Vector{Float64}} # (q,rate,edge,bl) -> exp(lambda*t)
const PTFullCache = Dict{Tuple{UInt64, Float64, Int, Float64}, Matrix{Float64}} # (q,rate,edge,bl) -> P(t)
const PTAnyCache = Union{PTDiagCache, PTFullCache}

@inline qid(q) = UInt(objectid(q.U))
@inline canon_rate(r::Float64) = (r == 0.0 ? 0.0 : r)
@inline logaddexp2(a::Float64, b::Float64) = a == -Inf ? b : (b == -Inf ? a : max(a, b) + log1p(exp(-abs(a - b))))

@inline function eigendecom_pt_diag!(
    out::AbstractVector{Float64},
    tmp::AbstractVector{Float64},
    U::AbstractMatrix{Float64},
    U_inv::AbstractMatrix{Float64},
    v::AbstractVector{Float64},
    diagexp_row
)
    mul!(tmp, U_inv, v)
    @inbounds @simd for k in eachindex(tmp)
        tmp[k] *= diagexp_row[k]
    end
    mul!(out, U, tmp)
    return out
end

function get_diagexp_cache!(cache::PTDiagCache, q, bl::Float64, r::Float64, edge::Int, nl::Int)
    key = (qid(q), canon_rate(r), edge, bl)
    get!(cache, key) do
        M = Vector{Float64}(undef, nl)
        Λ = q.Lambda
        t = bl * r
        @inbounds for s in 1:nl
            M[s] = exp(Λ[s] * t)
        end
        M
    end
end

function get_fullpt_cache!(cache::PTFullCache, q, bl::Float64, r::Float64, edge::Int, nl::Int)
    key = (qid(q), canon_rate(r), edge, bl)
    get!(cache, key) do
        U = q.U
        U_inv = q.U_inv
        Λ = q.Lambda

        # r == 0 => P = I for this edge
        if r == 0.0
            return Matrix{Float64}(I, nl, nl)
        end

        tmp = Matrix{Float64}(undef, nl, nl)  # row-scaled U_inv
        t = bl * r
        for i in 1:nl
            s = exp(Λ[i] * t)
            @simd for j in 1:nl
                tmp[i, j] = U_inv[i, j] * s
            end
        end
        P = Matrix{Float64}(undef, nl, nl)
        mul!(P, U, tmp)
        P
    end
end

######################################################################
function do_phylo_log_lk(
    ctx::RunCtx,
    q_pis,
    bl::Vector{Float64},
    liks_ori::Matrix{Float64};
    index::Int,
    pt_cache::Union{Nothing, PTAnyCache}=nothing,
    cache_mode::Symbol=:diag    # :diag or :full
)::Float64
    nb = ctx.nb
    nl = ctx.nl
    sub_model = ctx.sub_model
    rs = ctx.rs
    props = ctx.props
    Qrs = ctx.Qrs
    freqs = ctx.freqs
    is_do_inv = ctx.inv_info[:is_do_inv]
    inv_prop = ctx.inv_info[:inv_prop]
    all_children = ctx.all_children_vec

    props_norm = props ./ sum(props)
    is_lg4 = (sub_model == "LG4M" || sub_model == "LG4X")

    # IMPORTANT: +I should only use invariant mixture for invariant sites
    is_inv_site = is_do_inv ? is_invariant_site_from_liks(liks_ori, nb) : false

    liks = similar(liks_ori)
    comp = zeros(Float64, nb.tip + nb.node)
    childbuf = zeros(Float64, nl)
    m_prod = ones(Float64, nl)
    tmp = zeros(Float64, nl)

    # compute one component log-likelihood given q and scalar rate r
    function comp_loglk(q, r::Float64)
        Q = q.Q
        Pi = q.Pi
        U = q.U
        Lambda = q.Lambda
        U_inv = q.U_inv

        copyto!(liks, liks_ori)
        fill!(comp, 0.0)

        use_eig = (eltype(U) <: Float64 && eltype(Lambda) <: Float64 && eltype(U_inv) <: Float64)

        j = 0
        @inbounds for anc in (nb.tip + nb.node):-1:(1 + nb.tip)
            children = all_children[anc]
            fill!(m_prod, 1.0)

            for c in children
                j += 1
                @views childlik = liks[:, c]

                if r == 0.0
                    # P=I shortcut for invariant component
                    copyto!(childbuf, childlik)
                elseif !use_eig
                    t = bl[j] * r
                    childbuf .= exp(Q * t) * childlik
                else
                    if pt_cache !== nothing && cache_mode === :full
                        mul!(childbuf, get_fullpt_cache!(pt_cache::PTFullCache, q, bl[j], r, j, nl), childlik)
                    elseif pt_cache !== nothing && cache_mode === :diag
                        eigendecom_pt_diag!(childbuf, tmp, U, U_inv, childlik,
                                            get_diagexp_cache!(pt_cache::PTDiagCache, q, bl[j], r, j, nl))
                    else
                        t = bl[j] * r
                        eigendecom_pt!(childbuf, tmp, U, Lambda, U_inv, childlik, t)
                    end
                end

                @simd for s in 1:nl
                    m_prod[s] *= childbuf[s]
                end
            end

            comp[anc] = anc == (1 + nb.tip) ? dot(Pi, m_prod) : sum(m_prod)
            invc = 1.0 / max(comp[anc], 1e-300)
            @simd for s in 1:nl
                liks[s, anc] = m_prod[s] * invc
            end
        end

        lk = 0.0
        @inbounds for anc in (nb.tip + 1):(nb.tip + nb.node)
            lk += log(max(comp[anc], 1e-300))
        end
        return lk
    end

    if is_lg4
        qs = q_pis isa AbstractVector ? q_pis : (q_pis,)
        nq = length(qs)

        # LG4 component weights
        wk = length(props_norm) == nq ? props_norm :
             (length(freqs) == nq ? freqs ./ sum(freqs) : fill(1.0 / nq, nq))

        # variable part: sum_k wk * L(T, r_k Q_k)
        lks_var = Vector{Float64}(undef, nq)
        @inbounds for k in 1:nq
            rk = rs[min(k, length(rs))]
            qscale = length(Qrs) == nq ? Qrs[k] : 1.0
            lks_var[k] = comp_loglk(qs[k], qscale * rk)
        end
        log_var = logsumexp(log.(wk) .+ lks_var)

        if is_do_inv
            if is_inv_site
                # invariant part only for invariant sites
                lks_inv = Vector{Float64}(undef, nq)
                @inbounds for k in 1:nq
                    lks_inv[k] = comp_loglk(qs[k], 0.0)
                end
                log_inv = logsumexp(log.(wk) .+ lks_inv)

                a = log(max(inv_prop, 1e-300)) + log_inv
                b = log(max(1.0 - inv_prop, 1e-300)) + log_var
                return logsumexp([a, b])
            else
                # non-invariant sites: no invariant component
                return log(max(1.0 - inv_prop, 1e-300)) + log_var
            end
        else
            return log_var
        end
    else
        # non-LG4 behavior (cross product over rs and q), with corrected +I handling
        rs2 = rs
        props2 = props_norm
        if is_do_inv
            if is_inv_site
                rs2 = vcat(0.0, rs)
                props2 = vcat(inv_prop, props_norm .* (1 - inv_prop))
            else
                rs2 = rs
                props2 = props_norm .* (1 - inv_prop)
            end
        end

        qs = q_pis isa AbstractVector ? q_pis : (q_pis,)
        nq = length(qs)

        weights = vec((props2 .* freqs')')
        lks = zeros(Float64, length(weights))
        counter = 0

        @inbounds for r_i in eachindex(rs2), q_i in 1:nq
            counter += 1
            q = qs[q_i]
            r = (length(Qrs) == nq ? Qrs[q_i] : 1.0) * rs2[r_i]
            lks[counter] = comp_loglk(q, r)
        end

        return logsumexp(log.(weights) .+ lks)
    end
end

######################################################################
function get_liks_pattern(pattern, nl::Int, nb::Nb)
    out = Vector{Tuple{Matrix{Float64}, Float64}}(undef, length(pattern))
    for (k, x) in enumerate(pattern)
        liks = zeros(Float64, nl, nb.tip + nb.node)
        obs = x[1]
        for i in 1:nb.tip
            if obs[i] == 999
                @views liks[:, i] .= 1.0
            else
                liks[obs[i], i] = 1.0
            end
        end
        out[k] = (liks, Float64(x[2]))
    end
    return out
end

struct PatternBlock
    liks::Array{Float64, 3}      # state × node × pattern
    invmask::BitVector
    indices::UnitRange{Int}
end

function make_pattern_block(pattern2, nl::Int, nb::Nb, indices::UnitRange{Int}, is_do_inv::Bool)
    n = length(indices)
    liks = zeros(Float64, nl, nb.tip + nb.node, n)
    invmask = falses(n)

    @inbounds for (dst, src) in enumerate(indices)
        site_liks = pattern2[src][1]
        for t in 1:nb.tip, s in 1:nl
            liks[s, t, dst] = site_liks[s, t]
        end
        invmask[dst] = is_do_inv ? is_invariant_site_from_liks(site_liks, nb) : false
    end

    return PatternBlock(liks, invmask, indices)
end

function make_pattern_blocks(ctx::RunCtx, pattern2; block_size::Int=256)
    n_pat = length(pattern2)
    blocks = PatternBlock[]
    first_i = 1
    while first_i <= n_pat
        last_i = min(n_pat, first_i + block_size - 1)
        push!(blocks, make_pattern_block(pattern2, ctx.nl, ctx.nb, first_i:last_i, ctx.inv_info[:is_do_inv]))
        first_i = last_i + 1
    end
    return blocks
end

function comp_loglk_block!(
    out::Vector{Float64},
    ctx::RunCtx,
    q,
    bl::Vector{Float64},
    block::PatternBlock,
    r::Float64,
    pt_cache::Union{Nothing, PTAnyCache},
    cache_mode::Symbol
)
    nb = ctx.nb
    nl = ctx.nl
    all_children = ctx.all_children_vec
    n_pat = size(block.liks, 3)

    Q = q.Q
    Pi = q.Pi
    U = q.U
    Lambda = q.Lambda
    U_inv = q.U_inv

    liks = copy(block.liks)
    comp = zeros(Float64, nb.tip + nb.node, n_pat)
    childbuf = Matrix{Float64}(undef, nl, n_pat)
    tmp = Matrix{Float64}(undef, nl, n_pat)
    m_prod = ones(Float64, nl, n_pat)

    use_eig = (eltype(U) <: Float64 && eltype(Lambda) <: Float64 && eltype(U_inv) <: Float64)

    j = 0
    @inbounds for anc in (nb.tip + nb.node):-1:(1 + nb.tip)
        children = all_children[anc]
        fill!(m_prod, 1.0)

        for c in children
            j += 1
            childlik = @view liks[:, c, :]

            if r == 0.0
                copyto!(childbuf, childlik)
            elseif !use_eig
                t = bl[j] * r
                mul!(childbuf, exp(Q * t), childlik)
            elseif pt_cache !== nothing && cache_mode === :full
                mul!(childbuf, get_fullpt_cache!(pt_cache::PTFullCache, q, bl[j], r, j, nl), childlik)
            elseif pt_cache !== nothing && cache_mode === :diag
                mul!(tmp, U_inv, childlik)
                de = get_diagexp_cache!(pt_cache::PTDiagCache, q, bl[j], r, j, nl)
                for p in 1:n_pat, s in 1:nl
                    tmp[s, p] *= de[s]
                end
                mul!(childbuf, U, tmp)
            else
                t = bl[j] * r
                mul!(tmp, U_inv, childlik)
                for p in 1:n_pat, s in 1:nl
                    tmp[s, p] *= exp(Lambda[s] * t)
                end
                mul!(childbuf, U, tmp)
            end

            @simd for idx in eachindex(m_prod, childbuf)
                m_prod[idx] *= childbuf[idx]
            end
        end

        if anc == (1 + nb.tip)
            for p in 1:n_pat
                v = 0.0
                for s in 1:nl
                    v += Pi[s] * m_prod[s, p]
                end
                comp[anc, p] = v
            end
        else
            for p in 1:n_pat
                v = 0.0
                for s in 1:nl
                    v += m_prod[s, p]
                end
                comp[anc, p] = v
            end
        end

        for p in 1:n_pat
            invc = 1.0 / max(comp[anc, p], 1e-300)
            for s in 1:nl
                liks[s, anc, p] = m_prod[s, p] * invc
            end
        end
    end

    fill!(out, 0.0)
    @inbounds for p in 1:n_pat
        lk = 0.0
        for anc in (nb.tip + 1):(nb.tip + nb.node)
            lk += log(max(comp[anc, p], 1e-300))
        end
        out[p] = lk
    end

    return out
end

function pattern_loglk_vec_matrix(
    ctx::RunCtx,
    bl::Vector{Float64},
    blocks::Vector{PatternBlock};
    cache_mode::Symbol = :diag,
    caches = nothing
)
    n_pat = isempty(blocks) ? 0 : last(blocks[end].indices)
    out = Vector{Float64}(undef, n_pat)
    rs = ctx.rs
    props_norm = ctx.props ./ sum(ctx.props)
    freqs = ctx.freqs
    Qrs = ctx.Qrs
    is_do_inv = ctx.inv_info[:is_do_inv]
    inv_prop = ctx.inv_info[:inv_prop]
    is_lg4 = (ctx.sub_model == "LG4M" || ctx.sub_model == "LG4X")
    qs = ctx.q_pis isa AbstractVector ? ctx.q_pis : (ctx.q_pis,)
    nq = length(qs)

    nslots = Threads.maxthreadid()
    caches = caches === nothing ? [cache_mode === :full ? PTFullCache() : PTDiagCache() for _ in 1:nslots] : caches

    Threads.@threads for bi in eachindex(blocks)
        tid = Threads.threadid()
        cache = caches[tid]
        block = blocks[bi]
        n_block = length(block.indices)
        lks_var = fill(-Inf, n_block)
        lks_inv = fill(-Inf, n_block)
        comp = Vector{Float64}(undef, n_block)

        if is_lg4
            wk = length(props_norm) == nq ? props_norm :
                 (length(freqs) == nq ? freqs ./ sum(freqs) : fill(1.0 / nq, nq))

            for k in 1:nq
                rk = rs[min(k, length(rs))]
                qscale = length(Qrs) == nq ? Qrs[k] : 1.0
                comp_loglk_block!(comp, ctx, qs[k], bl, block, qscale * rk, cache, cache_mode)
                logw = log(max(wk[k], 1e-300))
                @inbounds for p in 1:n_block
                    lks_var[p] = logaddexp2(lks_var[p], logw + comp[p])
                end
            end

            if is_do_inv && any(block.invmask)
                for k in 1:nq
                    comp_loglk_block!(comp, ctx, qs[k], bl, block, 0.0, cache, cache_mode)
                    logw = log(max(wk[k], 1e-300))
                    @inbounds for p in 1:n_block
                        if block.invmask[p]
                            lks_inv[p] = logaddexp2(lks_inv[p], logw + comp[p])
                        end
                    end
                end
            end

            @inbounds for p in 1:n_block
                if is_do_inv
                    if block.invmask[p]
                        out[block.indices[p]] = logaddexp2(
                            log(max(inv_prop, 1e-300)) + lks_inv[p],
                            log(max(1.0 - inv_prop, 1e-300)) + lks_var[p]
                        )
                    else
                        out[block.indices[p]] = log(max(1.0 - inv_prop, 1e-300)) + lks_var[p]
                    end
                else
                    out[block.indices[p]] = lks_var[p]
                end
            end
        else
            if is_do_inv && any(block.invmask)
                for q_i in 1:nq
                    q = qs[q_i]
                    qscale = length(Qrs) == nq ? Qrs[q_i] : 1.0
                    comp_loglk_block!(comp, ctx, q, bl, block, 0.0, cache, cache_mode)
                    logw = log(max(inv_prop * freqs[q_i], 1e-300))
                    @inbounds for p in 1:n_block
                        if block.invmask[p]
                            lks_inv[p] = logaddexp2(lks_inv[p], logw + comp[p])
                        end
                    end
                end
            end

            for r_i in eachindex(rs), q_i in 1:nq
                q = qs[q_i]
                r = (length(Qrs) == nq ? Qrs[q_i] : 1.0) * rs[r_i]
                comp_loglk_block!(comp, ctx, q, bl, block, r, cache, cache_mode)
                logw = log(max(props_norm[r_i] * freqs[q_i] * (is_do_inv ? 1.0 - inv_prop : 1.0), 1e-300))
                @inbounds for p in 1:n_block
                    lks_var[p] = logaddexp2(lks_var[p], logw + comp[p])
                end
            end

            @inbounds for p in 1:n_block
                out[block.indices[p]] = (is_do_inv && block.invmask[p]) ?
                    logaddexp2(lks_inv[p], lks_var[p]) :
                    lks_var[p]
            end
        end
    end

    return out
end

function pattern_loglk_vec_matrix(
    ctx::RunCtx,
    bl::Vector{Float64},
    pattern2;
    cache_mode::Symbol = :diag,
    block_size::Int = 256,
    caches = nothing
)
    if ctx.is_pmsf
        return pattern_loglk_vec_batched(ctx, bl, pattern2; cache_mode=cache_mode, caches=caches)
    end
    return pattern_loglk_vec_matrix(ctx, bl, make_pattern_blocks(ctx, pattern2; block_size=block_size);
                                    cache_mode=cache_mode, caches=caches)
end

######################################################################
# IQ-TREE-style directional partials for the STK score matrix.  For edge
# parent--child, inside[:,child] describes the child side and outside[:,child]
# describes everything on the parent side.  Changing that edge therefore only
# requires outside' * P(t) * inside; no tree traversal is repeated.
function edge_component_logliks!(
    l0::Matrix{Float64}, lp::Matrix{Float64}, lm::Matrix{Float64},
    ctx::RunCtx, q, bl::Vector{Float64}, block::PatternBlock, r::Float64,
    cache::PTFullCache
)
    nb, nl = ctx.nb, ctx.nl
    nn = nb.tip + nb.node
    np = size(block.liks, 3)
    children = ctx.all_children_vec

    parent = zeros(Int, nn)
    edge = zeros(Int, nn)
    ej = 0
    @inbounds for anc in nn:-1:(nb.tip + 1)
        for c in children[anc]
            ej += 1; parent[c] = anc; edge[c] = ej
        end
    end
    ej == nb.branch || error("Tree traversal produced $ej edges, expected $(nb.branch)")

    inside = copy(block.liks)
    inscale = zeros(Float64, nn, np)
    work = Matrix{Float64}(undef, nl, np)
    prod = Matrix{Float64}(undef, nl, np)

    # child -> parent messages
    @inbounds for anc in nn:-1:(nb.tip + 1)
        fill!(prod, 1.0)
        scales = zeros(Float64, np)
        for c in children[anc]
            e = edge[c]
            P = get_fullpt_cache!(cache, q, bl[e], r, e, nl)
            mul!(work, P, @view(inside[:, c, :]))
            prod .*= work
            scales .+= @view(inscale[c, :])
        end
        for p in 1:np
            z = max(sum(@view(prod[:, p])), 1e-300)
            @views inside[:, anc, p] .= prod[:, p] ./ z
            inscale[anc, p] = scales[p] + log(z)
        end
    end

    # outside[:,node] is the complementary likelihood expressed at node's own
    # state. edgeoutside[:,child] is expressed at parent(child)'s state and is
    # therefore the left side of edgeoutside' * P * inside.
    outside = zeros(Float64, nl, nn, np)
    outscale = zeros(Float64, nn, np)
    edgeoutside = similar(outside)
    edgescale = similar(outscale)
    root = nb.tip + 1
    @inbounds for p in 1:np
        outside[:, root, p] .= q.Pi
    end

    tmp = Matrix{Float64}(undef, nl, np)
    @inbounds for anc in root:nn
        ch = children[anc]
        for focal in ch
            prod .= @view(outside[:, anc, :])
            scales = copy(@view(outscale[anc, :]))
            for sib in ch
                sib == focal && continue
                e2 = edge[sib]
                P2 = get_fullpt_cache!(cache, q, bl[e2], r, e2, nl)
                mul!(tmp, P2, @view(inside[:, sib, :]))
                prod .*= tmp
                scales .+= @view(inscale[sib, :])
            end
            for p in 1:np
                z = max(sum(@view(prod[:, p])), 1e-300)
                @views edgeoutside[:, focal, p] .= prod[:, p] ./ z
                edgescale[focal, p] = scales[p] + log(z)
            end
            # Move the complementary likelihood across anc--focal so it can
            # seed the next preorder level in focal's state space.
            Pf = get_fullpt_cache!(cache, q, bl[edge[focal]], r, edge[focal], nl)
            mul!(tmp, transpose(Pf), @view(edgeoutside[:, focal, :]))
            for p in 1:np
                z = max(sum(@view(tmp[:, p])), 1e-300)
                @views outside[:, focal, p] .= tmp[:, p] ./ z
                outscale[focal, p] = edgescale[focal, p] + log(z)
            end
        end
    end

    # One baseline contraction plus an analytic branch score per edge. Since
    # P(b)=exp(Q*r*b), dP/db = r*Q*P(b); no P(b+eps)/P(b-eps) matrices are
    # constructed here. `l0` receives the component log likelihood and `lp`
    # receives its component d(log L)/db score. `lm` is retained as scratch.
    @inbounds for child in 1:nn
        e = edge[child]
        e == 0 && continue
        P0 = get_fullpt_cache!(cache, q, bl[e], r, e, nl)
        mul!(work, P0, @view(inside[:, child, :]))
        if r != 0.0
            mul!(tmp, q.Q, work)
        else
            fill!(tmp, 0.0)
        end
        for p in 1:np
            v = max(dot(@view(edgeoutside[:, child, p]), @view(work[:, p])), 1e-300)
            dv = dot(@view(edgeoutside[:, child, p]), @view(tmp[:, p]))
            l0[p, e] = edgescale[child, p] + inscale[child, p] + log(v)
            lp[p, e] = r == 0.0 ? 0.0 : r * dv / v
            lm[p, e] = 0.0
        end
    end
    return nothing
end

function directional_score_matrix(bl::Vector{Float64}, pattern2, ctx::RunCtx,
                                  blocks::Vector{PatternBlock}; step_scale::Float64=1e6)
    npat, ne = length(pattern2), ctx.nb.branch
    L0 = fill(-Inf, npat, ne)
    D = zeros(Float64, npat, ne)
    props = ctx.props ./ sum(ctx.props)
    qs = ctx.q_pis isa AbstractVector ? ctx.q_pis : (ctx.q_pis,)
    nq = length(qs)
    is_lg4 = ctx.sub_model == "LG4M" || ctx.sub_model == "LG4X"
    caches = [PTFullCache() for _ in 1:Threads.maxthreadid()]

    Threads.@threads for bi in eachindex(blocks)
        block = blocks[bi]
        cache = caches[Threads.threadid()]
        empty!(cache)
        rows = block.indices
        a0 = fill(-Inf, length(rows), ne)
        score = zeros(Float64, length(rows), ne)
        comp_l = similar(a0); comp_d = similar(a0); scratch = similar(a0)
        if is_lg4
            wk = length(props) == nq ? props : (length(ctx.freqs) == nq ? ctx.freqs ./ sum(ctx.freqs) : fill(1/nq, nq))
            for k in 1:nq
                r = (length(ctx.Qrs) == nq ? ctx.Qrs[k] : 1.0) * ctx.rs[min(k, length(ctx.rs))]
                wt = wk[k] * (ctx.inv_info[:is_do_inv] ? 1-ctx.inv_info[:inv_prop] : 1.0)
                edge_component_logliks!(comp_l, comp_d, scratch, ctx, qs[k], bl, block, r, cache)
                logw = log(max(wt, 1e-300))
                @inbounds for p in axes(a0, 1), e in 1:ne
                    term = logw + comp_l[p, e]; old = a0[p, e]; new = logaddexp2(old, term)
                    score[p, e] = old == -Inf ? exp(term - new) * comp_d[p, e] :
                        exp(old - new) * score[p, e] + exp(term - new) * comp_d[p, e]
                    a0[p, e] = new
                end
            end
        else
            for ri in eachindex(ctx.rs), qi in 1:nq
                r = (length(ctx.Qrs) == nq ? ctx.Qrs[qi] : 1.0) * ctx.rs[ri]
                wt = props[ri] * ctx.freqs[qi] * (ctx.inv_info[:is_do_inv] ? 1-ctx.inv_info[:inv_prop] : 1.0)
                edge_component_logliks!(comp_l, comp_d, scratch, ctx, qs[qi], bl, block, r, cache)
                logw = log(max(wt, 1e-300))
                @inbounds for p in axes(a0, 1), e in 1:ne
                    term = logw + comp_l[p, e]; old = a0[p, e]; new = logaddexp2(old, term)
                    score[p, e] = old == -Inf ? exp(term - new) * comp_d[p, e] :
                        exp(old - new) * score[p, e] + exp(term - new) * comp_d[p, e]
                    a0[p, e] = new
                end
            end
        end
        # Invariable components have r=0, hence exactly the same value at +/-.
        if ctx.inv_info[:is_do_inv] && any(block.invmask)
            for qi in 1:nq
                iw = is_lg4 ? (length(props) == nq ? props[qi] : (length(ctx.freqs) == nq ? ctx.freqs[qi]/sum(ctx.freqs) : 1/nq)) : ctx.freqs[min(qi,end)]
                wt = ctx.inv_info[:inv_prop] * iw
                edge_component_logliks!(comp_l, comp_d, scratch, ctx, qs[qi], bl, block, 0.0, cache)
                logw = log(max(wt, 1e-300))
                for p in eachindex(block.invmask), e in 1:ne
                    block.invmask[p] || continue
                    term = logw + comp_l[p, e]; old = a0[p, e]; new = logaddexp2(old, term)
                    score[p, e] = old == -Inf ? exp(term - new) * comp_d[p, e] :
                        exp(old - new) * score[p, e] + exp(term - new) * comp_d[p, e]
                    a0[p, e] = new
                end
            end
        end
        L0[rows,:] .= a0; D[rows,:] .= score
    end
    return D, vec(L0[:,1])
end

# PMSF variant: Q differs by pattern, so patterns cannot be combined into one
# block. We still avoid the expensive branch-by-branch full-tree pruning: each
# pattern/rate component gets one inside pass, one outside pass, and all edge
# perturbations are evaluated from those cached messages.
function pmsf_directional_pattern(
    ctx::RunCtx, q_pis, bl::Vector{Float64}, liks_ori::Matrix{Float64},
    steps::Vector{Float64}, edge_for_child::Vector{Int},
    postorder::Vector{Int}, root::Int, cache::PTFullCache
)
    nb, nl = ctx.nb, ctx.nl
    nn = nb.tip + nb.node
    children = ctx.all_children_vec
    qs = q_pis isa AbstractVector ? q_pis : (q_pis,)
    nq = length(qs)
    is_inv = ctx.inv_info[:is_do_inv] ? is_invariant_site_from_liks(liks_ori, nb) : false
    props = ctx.props ./ sum(ctx.props)
    freqs = ctx.freqs
    steps2 = steps

    l0 = fill(-Inf, nb.branch); lp = similar(l0); lm = similar(l0)
    fill!(lp, -Inf); fill!(lm, -Inf)
    inside = Matrix{Float64}(undef, nl, nn)
    outside = Matrix{Float64}(undef, nl, nn)
    edgeoutside = Matrix{Float64}(undef, nl, nn)
    inscale = zeros(Float64, nn); outscale = zeros(Float64, nn); edgescale = zeros(Float64, nn)
    prod = zeros(Float64, nl); childbuf = zeros(Float64, nl); tmp = zeros(Float64, nl); work = zeros(Float64, nl)

    function one_component!(dest0, destp, destm, q, r, logweight)
        copyto!(inside, liks_ori); fill!(inscale, 0.0)
        U, Ui, Lam = q.U, q.U_inv, q.Lambda
        use_eig = eltype(U) <: Float64 && eltype(Ui) <: Float64 && eltype(Lam) <: Float64
        fill!(prod, 1.0)
        for anc in postorder
            fill!(prod, 1.0); scale = 0.0
            for c in children[anc]
                e = edge_for_child[c]
                if r == 0.0
                    copyto!(childbuf, @view(inside[:, c]))
                elseif use_eig
                    mul!(childbuf, get_fullpt_cache!(cache, q, bl[e], r, e, nl), @view(inside[:, c]))
                else
                    mul!(childbuf, exp(q.Q * (bl[e] * r)), @view(inside[:, c]))
                end
                @simd for s in 1:nl prod[s] *= childbuf[s] end
                scale += inscale[c]
            end
            z = max(sum(prod), 1e-300)
            @views inside[:, anc] .= prod ./ z
            inscale[anc] = scale + log(z)
        end

        fill!(outside, 0.0); fill!(outscale, 0.0)
        @views outside[:, root] .= q.Pi
        for anc in root:nn
            for focal in children[anc]
                copyto!(prod, @view(outside[:, anc])); scale = outscale[anc]
                for sib in children[anc]
                    sib == focal && continue
                    e = edge_for_child[sib]
                    P = get_fullpt_cache!(cache, q, bl[e], r, e, nl)
                    mul!(childbuf, P, @view(inside[:, sib]))
                    @simd for s in 1:nl prod[s] *= childbuf[s] end
                    scale += inscale[sib]
                end
                z = max(sum(prod), 1e-300)
                @views edgeoutside[:, focal] .= prod ./ z
                edgescale[focal] = scale + log(z)
                P = get_fullpt_cache!(cache, q, bl[edge_for_child[focal]], r, edge_for_child[focal], nl)
                mul!(tmp, transpose(P), @view(edgeoutside[:, focal]))
                z = max(sum(tmp), 1e-300)
                @views outside[:, focal] .= tmp ./ z
                outscale[focal] = edgescale[focal] + log(z)
            end
        end

        for child in 1:nn
            e = edge_for_child[child]; e == 0 && continue
            P0 = get_fullpt_cache!(cache, q, bl[e], r, e, nl)
            Pp = get_fullpt_cache!(cache, q, bl[e] + steps2[e], r, e, nl)
            Pm = get_fullpt_cache!(cache, q, bl[e] - steps2[e], r, e, nl)
            for (dest, P) in ((dest0, P0), (destp, Pp), (destm, Pm))
                mul!(work, P, @view(inside[:, child]))
                v = max(dot(@view(edgeoutside[:, child]), work), 1e-300)
                dest[e] = logaddexp2(dest[e], logweight + edgescale[child] + inscale[child] + log(v))
            end
        end
    end

    for qi in 1:nq
        q = qs[qi]
        for ri in eachindex(ctx.rs)
            r = (length(ctx.Qrs) == nq ? ctx.Qrs[qi] : 1.0) * ctx.rs[ri]
            wt = props[ri] * freqs[qi] * (ctx.inv_info[:is_do_inv] ? 1 - ctx.inv_info[:inv_prop] : 1.0)
            one_component!(l0, lp, lm, q, r, log(max(wt, 1e-300)))
        end
        if is_inv
            wt = ctx.inv_info[:inv_prop] * freqs[qi]
            one_component!(l0, lp, lm, q, 0.0, log(max(wt, 1e-300)))
        end
    end
    (l0, lp, lm)
end

function pmsf_directional_score_matrix(bl::Vector{Float64}, pattern2, ctx::RunCtx; step_scale::Float64=1e6)
    npat, ne = length(pattern2), ctx.nb.branch
    steps = (1 .+ bl) ./ step_scale
    D = Matrix{Float64}(undef, npat, ne); l0 = Vector{Float64}(undef, npat)
    nn = ctx.nb.tip + ctx.nb.node; edge_for_child = zeros(Int, nn); e = 0
    for anc in nn:-1:(ctx.nb.tip + 1), child in ctx.all_children_vec[anc]
        e += 1; edge_for_child[child] = e
    end
    postorder = collect((ctx.nb.tip + ctx.nb.node):-1:(ctx.nb.tip + 1))
    root = ctx.nb.tip + 1
    caches = [PTFullCache() for _ in 1:Threads.maxthreadid()]
    Threads.@threads for i in eachindex(pattern2)
        cache = caches[Threads.threadid()]
        empty!(cache)
        a0, ap, am = pmsf_directional_pattern(ctx, choose_q(ctx, i), bl, pattern2[i][1], steps,
                                               edge_for_child, postorder, root, cache)
        l0[i] = a0[1]
        @inbounds for j in 1:ne
            D[i, j] = (ap[j] - am[j]) / (2 * steps[j])
        end
    end
    D, l0
end

######################################################################
function delta_log_lk!(
    bls_work::Vector{Float64},
    bls_base::Vector{Float64},
    indices::Vector{Int},
    nums::Vector{Float64},
    i::Int,
    pattern2,
    ctx::RunCtx
)
    # reuse preallocated work vector
    copyto!(bls_work, bls_base)
    @inbounds for k in eachindex(indices, nums)
        bls_work[indices[k]] += nums[k]
    end

    site_pat, times_of_pattern = pattern2[i]
    q = choose_q(ctx, i)

    return do_phylo_log_lk(ctx, q, bls_work, site_pat; index = i) * times_of_pattern
end

######################################################################
function hessian_STK2004_fast(
    bls::Vector{Float64},
    pattern2,
    ctx::RunCtx;
    cache_mode::Symbol = :diag,
    fd_scheme::Symbol = :central,   # :central (default) or :forward
    step_scale::Float64 = 1e6,
    pattern_blocks::Union{Nothing, Vector{PatternBlock}} = nothing
)
    nb = ctx.nb
    n_pat = length(pattern2)
    w = Float64[p[2] for p in pattern2]

    shared_caches = [cache_mode === :full ? PTFullCache() : PTDiagCache() for _ in 1:Threads.maxthreadid()]
    loglk_vec = pattern_blocks === nothing ?
        ((x) -> pattern_loglk_vec_matrix(ctx, x, pattern2; cache_mode=cache_mode, caches=shared_caches)) :
        ((x) -> pattern_loglk_vec_matrix(ctx, x, pattern_blocks; cache_mode=cache_mode, caches=shared_caches))

    # baseline per-pattern lnL (unweighted)
    l0 = loglk_vec(bls)

    # D[k,i] = d lnL_k / d b_i  (unweighted per pattern)
    D = Matrix{Float64}(undef, n_pat, nb.branch)

    for i in 1:nb.branch
        step = (1 + bls[i]) / step_scale

        if fd_scheme === :forward
            blp = copy(bls); blp[i] += step
            fp = loglk_vec(blp)

            @inbounds @simd for k in 1:n_pat
                D[k, i] = (fp[k] - l0[k]) / step
            end

        elseif fd_scheme === :central
            blp = copy(bls); blp[i] += step
            blm = copy(bls); blm[i] -= step

            fp = loglk_vec(blp)
            fm = loglk_vec(blm)

            @inbounds @simd for k in 1:n_pat
                D[k, i] = (fp[k] - fm[k]) / (2 * step)
            end
        else
            error("Unknown fd_scheme=$(fd_scheme). Use :forward or :central.")
        end
    end

    # gradient of total lnL: g_i = Σ_k w_k * D[k,i]
    g = vec(transpose(D) * w)

    # STK OPG-style Hessian: H = - Σ_k w_k * d_k d_k'
    Dw = D .* reshape(w, :, 1)
    h = -(transpose(D) * Dw)

    lnL0 = dot(w, l0)
    return (h, g, lnL0, l0)
end

function hessian_Yang2000(
    bls::Vector{Float64}, pattern2, ctx::RunCtx;
    pattern_blocks::Union{Nothing, Vector{PatternBlock}}=nothing
)
    w = Float64[p[2] for p in pattern2]
    D, l0 = ctx.is_pmsf ?
        pmsf_directional_score_matrix(bls, pattern2, ctx) :
        directional_score_matrix(bls, pattern2, ctx,
            something(pattern_blocks, make_pattern_blocks(ctx, pattern2; block_size=256)))

    # The directional passes above are parallel over pattern blocks/sites.
    # Assemble the weighted OPG explicitly as well, because `julia -t` does
    # not control the separate BLAS thread pool used by D' * (w .* D).
    npat, ne = size(D)
    g = Vector{Float64}(undef, ne)
    h = Matrix{Float64}(undef, ne, ne)
    Threads.@threads for i in 1:ne
        gi = 0.0
        @inbounds @simd for k in 1:npat
            gi += w[k] * D[k, i]
        end
        g[i] = gi

        @inbounds for j in i:ne
            hij = 0.0
            @simd for k in 1:npat
                hij -= w[k] * D[k, i] * D[k, j]
            end
            h[i, j] = hij
        end
    end
    @inbounds for i in 2:ne, j in 1:(i - 1)
        h[i, j] = h[j, i]
    end
    return (h, g, dot(w, l0), l0)
end


######################################################################
function hessian_fd(bls::Vector{Float64}, sum_phylo_log_lk2::Function, nb_branch::Int)
    h = ones(nb_branch, nb_branch)
    function calculate_lnL(bls0, indices, nums)
        a = deepcopy(bls0)
        a[indices] .+= nums
        loglk = sum_phylo_log_lk2(a)
    end

    f0 = calculate_lnL(bls, [1], (0))
    f = calculate_lnL

    @inbounds Threads.@threads for i in 1:nb_branch
        @inbounds for j in 1:nb_branch
            if h[i, j] != 1.0
                println([i, j])
                h[i, j] = h[j, i]
                continue
            end
            bls0 = copy(bls) # bls is mutable
            △bls = max.((bls0[i], bls0[j]) ./ 1e6, 1e-6)
            (bl1, bl2) = △bls
            if i == j
                v1 = (f(bls, [i], (bl1)) - 2 * f0 + f(bls, [i], (-bl1))) / (reduce(*, △bls))
            else
                v1 = (f(bls0, [i, j], (bl1, bl2)) + f(bls0, [i, j], (-bl1, -bl2)) - f(bls0, [i, j], (-bl1, bl2)) - f(bls0, [i, j], (bl1, -bl2))) / (4 * reduce(*, △bls))
            end
            h[i, j] = v1
        end
    end

    println()
    # h .= -h

    return h
end

######################################################################
@inline function get_new_order(bl_order::Dict{String, Dict{Any, Any}})
    sorted_keys = sort(collect(keys(bl_order["mcmctree2ape"])))
    [bl_order["mcmctree2ape"][k] for k in sorted_keys]
end

@inline function get_new_order_rev(bl_order::Dict{String, Dict{Any, Any}})
    sorted_keys = sort(collect(keys(bl_order["ape2mcmctree"])))
    [bl_order["ape2mcmctree"][k] for k in sorted_keys]
end

function pattern_loglk_vec_batched(
    ctx::RunCtx,
    bl::Vector{Float64},
    pattern2;
    cache_mode::Symbol = :diag,
    caches = nothing
)
    n_pat = length(pattern2)
    out = Vector{Float64}(undef, n_pat)
    nslots = Threads.maxthreadid()
    caches = caches === nothing ? [cache_mode === :full ? PTFullCache() : PTDiagCache() for _ in 1:nslots] : caches

    # Schedule contiguous blocks.  This keeps each worker on a compact range of
    # patterns and, importantly, retains one transition cache for the whole
    # block instead of rebuilding a cache per pattern/evaluation.
    block_size = max(1, min(256, cld(n_pat, max(1, 4 * nslots))))
    ranges = UnitRange{Int}[]
    first_i = 1
    while first_i <= n_pat
        last_i = min(n_pat, first_i + block_size - 1)
        push!(ranges, first_i:last_i)
        first_i = last_i + 1
    end

    Threads.@threads for ri in eachindex(ranges)
        tid = Threads.threadid()
        cache = caches[tid]
        for i in ranges[ri]
            site_pat = pattern2[i][1]
            q = choose_q(ctx, i)
            out[i] = do_phylo_log_lk(ctx, q, bl, site_pat;
                                     index=i, pt_cache=cache,
                                     cache_mode=cache_mode)
        end
    end
    return out
end

function sum_phylo_log_lk_fast(
    ctx::RunCtx,
    bl::Vector{Float64},
    pattern2;
    cache_mode::Symbol = :diag,
    pattern_blocks::Union{Nothing, Vector{PatternBlock}} = nothing
)
    lks = pattern_blocks === nothing ?
        pattern_loglk_vec_matrix(ctx, bl, pattern2; cache_mode=cache_mode) :
        pattern_loglk_vec_matrix(ctx, bl, pattern_blocks; cache_mode=cache_mode)
    w = Float64[p[2] for p in pattern2]
    return dot(w, lks)
end

function sum_phylo_log_lk(ctx::RunCtx, param::Vector{Float64}, pattern2; cache_mode::Symbol=:diag)
    s = 0.0
    pt_cache = cache_mode === :full ? PTFullCache() : PTDiagCache()  # valid for this fixed param vector
    for i in eachindex(pattern2)
        q = choose_q(ctx, i)
        s += do_phylo_log_lk(ctx, q, param, pattern2[i][1];
                             index=i, pt_cache=pt_cache, cache_mode=cache_mode) * pattern2[i][2]
    end
    return s
end

function calculate_gradient(bls::Vector{Float64}, sum_phylo_log_lk2::Function)
    delta_loglks = Vector{Float64}(undef, length(bls))
    bls2 = similar(bls)
    copyto!(bls2, bls)
    @inbounds for i in eachindex(bls2)
        if bls2[i] <= 1e-6
            bls2[i] = 0.0
        end
    end

    for i in eachindex(bls2)
        step = (1 + bls2[i]) / 1e6
        a = copy(bls2)
        a[i] += step
        b = copy(bls2)
        b[i] -= step
        delta_loglks[i] = (sum_phylo_log_lk2(a) - sum_phylo_log_lk2(b)) / (2 * step)
    end
    return delta_loglks
end

######################################################################
function transform_bl(bls, g, h)
    u = sqrt.(bls)
    db_over_du = 2 .* u
    d2b_over_du2 = fill(2.0, length(u))

    gu = g .* db_over_du
    hu = Diagonal(g .* d2b_over_du2) .+ h .* (db_over_du * db_over_du')
    return (u, gu, hu)
end

function compare_true_vs_approx_lnL(bls_vec::Vector{Float64}, gs, hs, lnL0, sum_phylo_log_lk2, transform_method)
    hs2 = [x for x in hs if x !== nothing]
    bls = bls_vec
    covs = map(x -> -inv(Hermitian(x + Diagonal(ones(size(x, 1))) * 1e-6)) * 2, hs2)

    # NOTE: uses MvNormal/Uniform if you bring Distributions back
    for i in 1:100
        new_bls = copy(bls)
        lnL = sum_phylo_log_lk2(new_bls)

        lnL_approxes = Float64[]
        for (g, h) in zip(gs, hs2)
            if transform_method != nothing
                bls2 = sqrt.(bls)
                new_bls2, gu, hu = transform_bl(new_bls, g, h)
                delta_bls = new_bls2 - bls2
                lnL_approx = lnL0 + (delta_bls' * gu) + (delta_bls' * hu * delta_bls / 2)
            else
                delta_bls = new_bls - bls
                lnL_approx = lnL0 + (delta_bls' * g) + (delta_bls' * h * delta_bls / 2)
            end
            push!(lnL_approxes, lnL_approx)
        end
        println(join(vcat("TREE$i", lnL, lnL_approxes...), "\t"))
    end
end

function compare_true_vs_approx_lnL(bls_vec::Vector{Vector{Float64}}, gs, hs, lnL0, sum_phylo_log_lk2, transform_method)
    hs2 = [x for x in hs if x !== nothing]
    bls = popfirst!(bls_vec)

    for (idx, bs_bls) in enumerate(bls_vec[1:min(100, length(bls_vec))])
        new_bls = copy(bs_bls)
        lnL = sum_phylo_log_lk2(new_bls)

        lnL_approxes = Float64[]
        for (g, h) in zip(gs, hs2)
            if transform_method != nothing
                bls2 = sqrt.(bls)
                new_bls2, gu, hu = transform_bl(new_bls, g, h)
                delta_bls = new_bls2 - bls2
                lnL_approx = lnL0 + (delta_bls' * gu) + (delta_bls' * hu * delta_bls / 2)
            else
                delta_bls = new_bls - bls
                lnL_approx = lnL0 + (delta_bls' * g) + (delta_bls' * h * delta_bls / 2)
            end
            push!(lnL_approxes, lnL_approx)
        end
        println(join(vcat("TREE$idx", lnL, lnL_approxes...), "\t"))
    end
end

######################################################################
function read_hessian_infile(hessian_infile)
    in_fh = open(hessian_infile, "r")
    is_read_hessian = false
    hessian_data = Float64[]

    for line in readlines(in_fh)
        line = strip(line)
        if is_read_hessian
            if !occursin(r"\w", line)
                continue
            else
                nums = parse.(Float64, split(line, r"\s+"))
                append!(hessian_data, nums)
            end
        end
        if occursin("Hessian", line)
            is_read_hessian = true
        end
    end
    close(in_fh)

    n = Int(sqrt(length(hessian_data)))
    if n^2 != length(hessian_data)
        @warn "Hessian's length is not a square of a integer!"
        return reshape(hessian_data, :, length(hessian_data) ÷ n)
    end
    return reshape(hessian_data, n, n)
end

function get_h_from_inBV(hessian_infile, bl_order)
    h_from_inBV_mcmctree = hessian_infile === nothing ? nothing : read_hessian_infile(hessian_infile)
    if h_from_inBV_mcmctree === nothing
        return nothing
    end
    new_bl_order_rev = get_new_order_rev(bl_order)
    return h_from_inBV_mcmctree[new_bl_order_rev, new_bl_order_rev]
end

######################################################################
function main(ctx::RunCtx)
    nb = ctx.nb
    pattern2 = get_liks_pattern(ctx.pattern, ctx.nl, nb)
    pattern_blocks = ctx.is_pmsf ? nothing : make_pattern_blocks(ctx, pattern2; block_size=256)

    bls_vec = Vector{Vector{Float64}}()
    bls = read_bls_from_branchout_matrix(ctx.branchout_matrix)[1]
    if ctx.bs_branchout_matrix !== nothing
        bls_vec = read_bls_from_branchout_matrix(ctx.bs_branchout_matrix)
    end

    bl_order = get_bl_order(ctx.branchout_matrix)
    new_bl_order = get_new_order(bl_order)

    println(now())

    cache_mode = ctx.cache_mode
    println("cache_mode = ", cache_mode)
    sum_phylo_log_lk2 = (x::Vector{Float64}) -> sum_phylo_log_lk_fast(ctx, x, pattern2;
        cache_mode=cache_mode, pattern_blocks=pattern_blocks)

    lnL0 = sum_phylo_log_lk2(bls)
    println("Starting Hessian and gradient calculation.")
    println(now())

    if ctx.hessian_outfile !== nothing
        # keep your existing gradient behavior
        g = zeros(length(bls))

        lnL_for_compare = lnL0
        htype = lowercase(ctx.hessian_type)
        if occursin(r"^(finite_difference|fd|2nd_order|2nd_order_derivative|2nd)$", htype)
            h = hessian_fd(bls, sum_phylo_log_lk2, nb.branch)
            _, g_stk, lnL0_stk, _ = hessian_STK2004_fast(bls, pattern2, ctx;
                cache_mode=cache_mode, fd_scheme=ctx.fd_scheme, pattern_blocks=pattern_blocks)
            g = g_stk
            # g = calculate_gradient(bls, sum_phylo_log_lk2)
        elseif htype == "stk2004"
            h, g_stk, lnL0_stk, _ = hessian_STK2004_fast(bls, pattern2, ctx;
                cache_mode=cache_mode, fd_scheme=ctx.fd_scheme, pattern_blocks=pattern_blocks)
            g = g_stk
            lnL_for_compare = lnL0_stk
        elseif htype == "yang2000"
            h, g_yang, lnL0_yang, _ = hessian_Yang2000(bls, pattern2, ctx;
                pattern_blocks=pattern_blocks)
            g = g_yang
            lnL_for_compare = lnL0_yang
        else
            error("Unknown --hessian_type $(ctx.hessian_type). Use STK2004, Yang2000, or fd.")
        end
        g_mcmctree = g[new_bl_order]

        out_fh = open(ctx.hessian_outfile, "w")
        println(out_fh, "")
        println(out_fh, string(nb.tip))
        println(out_fh, "")
        println(out_fh, read(ctx.treefile, String))
        println(out_fh, join(bls[new_bl_order], " "))
        println(out_fh, "")
        println(out_fh, join(g_mcmctree, " "))
        write(out_fh, "\n\n")
        println(out_fh, "Hessian")
        println(out_fh, "")
        close(out_fh)

        h_mcmctree = h[new_bl_order, new_bl_order]
        h_mcmctree_out = format_number.(h_mcmctree)
        open(ctx.hessian_outfile, "a") do fh
            writedlm(fh, h_mcmctree_out, ' ')
        end

        if ctx.is_compare
            h_from_inBV = get_h_from_inBV(ctx.hessian_infile, bl_order)
            if isempty(bls_vec)
                compare_true_vs_approx_lnL(bls, [g], [h, h_from_inBV], lnL_for_compare, sum_phylo_log_lk2, ctx.transform_method)
            else
                gs = [g, zeros(Float64, length(g))]
                compare_true_vs_approx_lnL(vcat([bls], bls_vec), gs, [h, h_from_inBV], lnL_for_compare, sum_phylo_log_lk2, ctx.transform_method)
            end
        end
    end

    println(join(["best lnL", lnL0, "\n"], "\t"))
    println(ctx.std_outfh, join(["best lnL", lnL0, "\n"], "\t"))

    close(ctx.std_outfh)
    println(now())
    println()
end

######################################################################
function get_args_refactored()
    println()
    println(now())

    # from your existing arg/config environment (readArg.jl)
    rs, props, Fs, Qrs, freqs, inv_info = get_iqtree_params(iqtree_file, phyml_file)

    println(now())
    nb, pattern, all_children, cherry_nodes, descendants, site2pattern = read_basics(basics_indir)
    nl = st == "AA" ? 20 : 4

    Qrs, q_pis, q_pis_sites, freqs = get_Qrs_freqs(
        Fs, Qrs, freqs, is_pmsf, pmsf_file, site2pattern, sub_model, mix_freq_model;
        seq_type=type, iqtree_file=iqtree_file
    )

    all_children_vec = dict_children_to_vec(all_children, nb.tip + nb.node)

    ctx = RunCtx(
        nb = nb,
        nl = nl,
        sub_model = sub_model,
        rs = rs,
        props = props,
        Fs = Fs,
        Qrs = Qrs,
        freqs = freqs,
        q_pis = q_pis,
        q_pis_sites = q_pis_sites,
        is_pmsf = is_pmsf,
        inv_info = inv_info,
        pattern = pattern,
        all_children_vec = all_children_vec,
        descendants = descendants,

        hessian_outfile = hessian_outfile,
        hessian_type = hessian_type,
        hessian_infile = hessian_infile,
		fd_scheme = fd_scheme,
        cache_mode = cache_mode,
        is_compare = is_compare,
        transform_method = transform_method,
        treefile = treefile,
        std_outfh = std_outfh,
        branchout_matrix = branchout_matrix,
        bs_branchout_matrix = bs_branchout_matrix
    )

    return ctx
end

######################################################################
if abspath(PROGRAM_FILE) == @__FILE__
    ctx = get_args_refactored()
    main(ctx)
end
