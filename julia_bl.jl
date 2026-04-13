#! /bin/env julia

######################################################################
# v0.24.1 (fix +I log-likelihood: apply invariant mixture only to invariant sites)

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
function do_phylo_log_lk(
    ctx::RunCtx,
    q_pis,
    bl::Vector{Float64},
    liks_ori::Matrix{Float64};
    index::Int
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

        j = 0
        @inbounds for anc in (nb.tip + nb.node):-1:(1 + nb.tip)
            children = all_children[anc]
            fill!(m_prod, 1.0)

            for c in children
                j += 1
                t = bl[j] * r
                @views childlik = liks[:, c]

                if !(eltype(U) <: Float64 && eltype(Lambda) <: Float64 && eltype(U_inv) <: Float64)
                    childbuf .= exp(Q * t) * childlik
                else
                    eigendecom_pt!(childbuf, tmp, U, Lambda, U_inv, childlik, t)
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
function hessian_STK2004(
    bls::Vector{Float64},
    pattern2,
    ctx::RunCtx
)
    nb = ctx.nb
    n_pat = length(pattern2)

    lnLs = zeros(Float64, n_pat)
    h_local = [zeros(Float64, nb.branch, nb.branch) for _ in 1:nthreads()]
    g_local = [zeros(Float64, nb.branch) for _ in 1:nthreads()]

    # thread-local bls work buffers
    bls_work = [similar(bls) for _ in 1:nthreads()]

    Threads.@threads for k in 1:n_pat
        tid = threadid()
        h = h_local[tid]
        g = g_local[tid]
        bw = bls_work[tid]

        lnLs[k] = delta_log_lk!(bw, bls, [1], [0.0], k, pattern2, ctx)

        thetas = zeros(Float64, nb.branch)
        @inbounds for i in 1:nb.branch
            step = (1 + bls[i]) / 1e6
            fplus  = delta_log_lk!(bw, bls, [i], [ step], k, pattern2, ctx)
            fminus = delta_log_lk!(bw, bls, [i], [-step], k, pattern2, ctx)
            thetas[i] = (fplus - fminus) / (2 * step)
        end

        g .+= thetas
        alpha = 1.0 / pattern2[k][2]
        BLAS.ger!(alpha, thetas, thetas, h)
    end

    hsum = zeros(Float64, nb.branch, nb.branch)
    gsum = zeros(Float64, nb.branch)
    for t in 1:nthreads()
        hsum .+= h_local[t]
        gsum .+= g_local[t]
    end

    h = -hsum
    lnL0 = sum(lnLs)
    return (h, gsum, lnL0, lnLs)
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

function sum_phylo_log_lk(ctx::RunCtx, param::Vector{Float64}, pattern2)
    s = 0.0
    for i in eachindex(pattern2)
        q = choose_q(ctx, i)
        s += do_phylo_log_lk(ctx, q, param, pattern2[i][1]; index = i) * pattern2[i][2]
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

    bls_vec = Vector{Vector{Float64}}()
    bls = read_bls_from_branchout_matrix(ctx.branchout_matrix)[1]
    if ctx.bs_branchout_matrix !== nothing
        bls_vec = read_bls_from_branchout_matrix(ctx.bs_branchout_matrix)
    end

    bl_order = get_bl_order(ctx.branchout_matrix)
    new_bl_order = get_new_order(bl_order)

    println(now())

    sum_phylo_log_lk2 = (x::Vector{Float64}) -> sum_phylo_log_lk(ctx, x, pattern2)

    lnL0 = sum_phylo_log_lk2(bls)

    if ctx.hessian_outfile !== nothing
        h, g_stk, lnL0_stk, _ = hessian_STK2004(bls, pattern2, ctx)
        g = calculate_gradient(bls, sum_phylo_log_lk2)  # keep your original behavior
        g_mcmctree = g[new_bl_order]

        if occursin(r"^(finite_difference|fd|2nd_order|2nd_order_derivative)$", ctx.hessian_type)
            h = hessian_fd(bls, sum_phylo_log_lk2, nb.branch)
        end

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
                compare_true_vs_approx_lnL(bls, [g], [h, h_from_inBV], lnL0_stk, sum_phylo_log_lk2, ctx.transform_method)
            else
                gs = [g, zeros(Float64, length(g))]
                compare_true_vs_approx_lnL(vcat([bls], bls_vec), gs, [h, h_from_inBV], lnL0_stk, sum_phylo_log_lk2, ctx.transform_method)
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
    nl = type == "AA" ? 20 : 4

    Qrs, q_pis, q_pis_sites, freqs = get_Qrs_freqs(
        Fs, Qrs, freqs, is_pmsf, pmsf_file, site2pattern, sub_model, mix_freq_model
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
