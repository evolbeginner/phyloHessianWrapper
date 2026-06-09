#!/usr/bin/env julia
# nex_split_export.jl
#
# Usage:
#   julia nex_split_export.jl --model_nex models.nex --outdir out
#   julia nex_split_export.jl --model_nex models.nex --outdir out --force
#   julia nex_split_export.jl --model_nex models.nex --outdir out EX2 C10 LG4X
#
# Behavior:
# - Exports one .dat per top-level model (or selected top-level models).
# - Models WITH exchangeability matrix -> outdir/regular
# - Models WITHOUT exchangeability matrix -> outdir/mfm
# - --force removes and recreates outdir if it exists.

using Printf

struct Def
    kind::Symbol      # :model or :frequency
    name::String
    rhs::String
    pos::Int
end

const RESERVED = Set([
    "MIX","FMIX","POISSON","G4","R4",
    "JTT","LG","WAG","DAYHOFF","BLOSUM62",
    "empirical"
])

# -----------------------------
# CLI parsing
# -----------------------------
function parse_cli(args::Vector{String})
    model_nex = nothing
    outdir = nothing
    force = false
    targets = String[]

    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--model_nex"
            i == length(args) && error("Missing value for --model_nex")
            model_nex = args[i+1]
            i += 2
        elseif a == "--outdir"
            i == length(args) && error("Missing value for --outdir")
            outdir = args[i+1]
            i += 2
        elseif a == "--force"
            force = true
            i += 1
        elseif startswith(a, "--")
            error("Unknown option: $a")
        else
            # positional model target
            push!(targets, a)
            i += 1
        end
    end

    model_nex === nothing && error("Required option missing: --model_nex")
    outdir === nothing && error("Required option missing: --outdir")
    return (String(model_nex), String(outdir), force, targets)
end

# -----------------------------
# NEXUS parsing helpers
# -----------------------------
function strip_nexus_comments(s::AbstractString)
    io = IOBuffer()
    depth = 0
    for c in s
        if c == '['
            depth += 1
        elseif c == ']'
            depth = max(depth - 1, 0)
        elseif depth == 0
            write(io, c)
        end
    end
    return String(take!(io))
end

function extract_models_block(s::AbstractString)
    m = match(r"(?is)\bbegin\s+models\s*;\s*(.*?)\s*\bend\s*;", s)
    return m === nothing ? "" : m.captures[1]
end

function parse_defs(block::AbstractString)
    rx = r"(?is)\b(model|frequency)\s+([A-Za-z_][A-Za-z0-9_]*)\s*=\s*(.*?);"
    defs = Dict{String,Def}()
    for m in eachmatch(rx, block)
        kind = lowercase(m.captures[1]) == "model" ? :model : :frequency
        name = String(m.captures[2])
        rhs  = strip(String(m.captures[3]))
        defs[name] = Def(kind, name, rhs, first(m.offset))
    end
    return defs
end

ident_tokens(s::AbstractString) = [String(m.match) for m in eachmatch(r"[A-Za-z_][A-Za-z0-9_]*", s)]

function build_dep_graph(defs::Dict{String,Def})
    deps = Dict{String,Vector{String}}()
    for (name, d) in defs
        refs = String[]
        for t in ident_tokens(d.rhs)
            if t == name || t in RESERVED
                continue
            end
            if haskey(defs, t)
                push!(refs, t)
            end
        end
        deps[name] = unique(refs)
    end
    return deps
end

# top-level models = models not referenced by another model
function top_level_models(defs::Dict{String,Def}, deps::Dict{String,Vector{String}})
    model_names = [n for (n,d) in defs if d.kind == :model]
    referenced = Set{String}()
    for (_, vs) in deps
        for v in vs
            if haskey(defs, v) && defs[v].kind == :model
                push!(referenced, v)
            end
        end
    end
    roots = [m for m in model_names if !(m in referenced)]
    sort!(roots, by = n -> defs[n].pos)
    return roots
end

function topo_for_target(target::String, deps::Dict{String,Vector{String}})
    seen = Set{String}()
    temp = Set{String}()
    out  = String[]

    function dfs(u::String)
        if u in seen
            return
        end
        if u in temp
            error("Cycle detected involving '$u'")
        end
        push!(temp, u)
        for v in get(deps, u, String[])
            dfs(v)
        end
        delete!(temp, u)
        push!(seen, u)
        push!(out, u)
    end

    dfs(target)
    return out
end

# -----------------------------
# Classification logic
# -----------------------------
# "Numeric model rhs" means model body is essentially just numbers/separators
# => exchangeability matrix + frequencies style model block.
function is_numeric_model_rhs(rhs::AbstractString)
    # remove all numeric tokens
    t = replace(rhs, r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?" => " ")
    # remove separators/operators
    t = replace(t, r"[\s,;:{}()\[\]*+\-\/]" => "")
    # if any letters remain -> not pure numeric model
    return !occursin(r"[A-Za-z_]", t)
end

function has_exchangeability_matrix(closure_names::Vector{String}, defs::Dict{String,Def})
    for n in closure_names
        d = defs[n]
        if d.kind == :model && is_numeric_model_rhs(d.rhs)
            return true
        end
    end
    return false
end

# -----------------------------
# Writer
# -----------------------------
function write_dat(path::String, names::Vector{String}, defs::Dict{String,Def})
    open(path, "w") do io
        for n in names
            d = defs[n]
            println(io, string(d.kind), " ", d.name, " = ", d.rhs, ";")
        end
    end
end

# -----------------------------
# Main
# -----------------------------
function main()
    model_nex, outdir, force, requested = parse_cli(ARGS)

    isfile(model_nex) || error("Input file not found: $model_nex")

    if isdir(outdir) || isfile(outdir)
        if force
            rm(outdir; recursive=true, force=true)
        else
            error("Outdir already exists. Use --force to recreate: $outdir")
        end
    end

    mkpath(outdir)
    regular_dir = joinpath(outdir, "regular")
    mfm_dir     = joinpath(outdir, "mfm")
    mkpath(regular_dir)
    mkpath(mfm_dir)

    raw = read(model_nex, String)
    txt = strip_nexus_comments(raw)
    block = extract_models_block(txt)
    isempty(strip(block)) && error("No 'begin models; ... end;' block found in $model_nex")

    defs = parse_defs(block)
    isempty(defs) && error("No model/frequency definitions parsed.")

    deps = build_dep_graph(defs)
    roots = top_level_models(defs, deps)
    isempty(roots) && error("No top-level models found.")

    targets = if isempty(requested)
        roots
    else
        for t in requested
            haskey(defs, t) || error("Requested '$t' not found.")
            defs[t].kind == :model || error("Requested '$t' is not a model.")
            t in roots || error("Requested '$t' is not top-level (internal component).")
        end
        requested
    end

    println("Export targets: ", join(targets, ", "))

    n_regular = 0
    n_mfm = 0

    for t in targets
        closure_names = topo_for_target(t, deps)
        is_regular = has_exchangeability_matrix(closure_names, defs)
        outpath = joinpath(is_regular ? regular_dir : mfm_dir, "$t.dat")
        write_dat(outpath, closure_names, defs)

        if is_regular
            n_regular += 1
            @info "Wrote regular model" model=t path=outpath
        else
            n_mfm += 1
            @info "Wrote mfm model" model=t path=outpath
        end
    end

    println("\nDone.")
    println("  regular: $n_regular")
    println("  mfm:     $n_mfm")
    println("  outdir:  $outdir")
end

main()
