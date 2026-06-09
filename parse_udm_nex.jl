#!/usr/bin/env julia

# parse_udm_nex.jl
#
# Usage:
#   julia parse_udm_nex.jl --model_nex UDM_clr_iqtree_merged.nex --outdir UDM_clr --force

module UDMNexParser

using Printf
using Logging

# ---------------------------
# Data structures
# ---------------------------
struct FmixEntry
    component::String
    scale::Union{Nothing, Float64}
    weight::Union{Nothing, Float64}
end

mutable struct ModelBlock
    source_file::Union{Nothing, String}
    profiles::Dict{String, Vector{Float64}}               # frequency NAME = numbers...
    fmix::Dict{String, Vector{FmixEntry}}                 # frequency NAME = FMIX{...}
end

ModelBlock(source_file::Union{Nothing, String}=nothing) = ModelBlock(
    source_file,
    Dict{String, Vector{Float64}}(),
    Dict{String, Vector{FmixEntry}}()
)

# ---------------------------
# Helpers
# ---------------------------
strip_trailing_semicolon(s::AbstractString) = endswith(strip(s), ";") ? strip(chop(strip(s))) : strip(s)

function parse_source_file_line(line::AbstractString)
    m = match(r"^\[\s*SOURCE_FILE:\s*(.*?)\s*\]$", strip(line))
    return m === nothing ? nothing : m.captures[1]
end

function parse_float_token(tok::AbstractString)
    t = strip(tok)
    t = replace(t, "," => "")  # just in case
    return tryparse(Float64, t)
end

function parse_profile_rhs(rhs::AbstractString)
    vals = Float64[]
    toks = split(strip(rhs))
    for t in toks
        x = parse_float_token(t)
        if x === nothing
            @warn "Could not parse numeric token in profile RHS" token=t rhs=rhs
            return nothing
        end
        push!(vals, x)
    end
    return vals
end

function parse_fmix_rhs(rhs::AbstractString)
    # expects FMIX{...}
    r = strip(rhs)
    if !startswith(uppercase(r), "FMIX{") || !endswith(r, "}")
        @warn "Malformed FMIX RHS" rhs=rhs
        return nothing
    end

    inside = r[firstindex(r)+5:lastindex(r)-1]  # remove "FMIX{" and trailing "}"
    items = split(inside, ",")
    entries = FmixEntry[]

    for raw in items
        token = strip(raw)
        isempty(token) && continue

        parts = split(token, ":")
        component = strip(parts[1])

        scale::Union{Nothing, Float64} = nothing
        weight::Union{Nothing, Float64} = nothing

        if length(parts) >= 2
            scale = tryparse(Float64, strip(parts[2]))
            if length(parts) >= 3
                weight = tryparse(Float64, strip(parts[3]))
            end
        end

        push!(entries, FmixEntry(component, scale, weight))
    end

    return entries
end

function parse_frequency_statement(stmt::AbstractString)
    # stmt without trailing ';'
    # frequency NAME = RHS
    m = match(r"^frequency\s+([A-Za-z0-9_]+)\s*=\s*(.+)$"i, strip(stmt))
    if m === nothing
        return nothing
    end
    name = m.captures[1]
    rhs  = strip(m.captures[2])

    if startswith(uppercase(rhs), "FMIX{")
        fmix_entries = parse_fmix_rhs(rhs)
        return (:fmix, name, fmix_entries)
    else
        vals = parse_profile_rhs(rhs)
        return (:profile, name, vals)
    end
end

# ---------------------------
# Core parser
# ---------------------------
function parse_nex_models(path::String)
    blocks = ModelBlock[]
    pending_source::Union{Nothing, String} = nothing
    current_block::Union{Nothing, ModelBlock} = nothing

    in_models = false
    stmt_buf = IOBuffer()

    function flush_statement!()
        s = String(take!(stmt_buf)) |> strip
        isempty(s) && return
        s = strip_trailing_semicolon(s)
        parsed = parse_frequency_statement(s)
        if parsed === nothing
            return
        end
        kind, name, payload = parsed
        if current_block === nothing
            @warn "Found frequency statement outside a models block; ignoring" statement=s
            return
        end
        if kind == :profile
            if payload === nothing
                @warn "Skipping malformed profile" name=name statement=s
            else
                current_block.profiles[name] = payload
            end
        elseif kind == :fmix
            if payload === nothing
                @warn "Skipping malformed FMIX" name=name statement=s
            else
                current_block.fmix[name] = payload
            end
        end
    end

    open(path, "r") do io
        for rawline in eachline(io)
            line = strip(rawline)

            # Track SOURCE_FILE marker
            sf = parse_source_file_line(line)
            if sf !== nothing
                pending_source = sf
                continue
            end

            # Skip comments/empty lines
            isempty(line) && continue
            startswith(line, "#") && continue

            # Begin / End models
            if occursin(r"^begin\s+models\s*;"i, line)
                in_models = true
                current_block = ModelBlock(pending_source)
                pending_source = nothing
                continue
            elseif occursin(r"^end\s*;"i, line) && in_models
                # flush any dangling statement (best-effort)
                if position(stmt_buf) > 0
                    flush_statement!()
                end
                push!(blocks, current_block)
                current_block = nothing
                in_models = false
                continue
            end

            # Ignore everything outside models blocks
            if !in_models
                continue
            end

            # Accumulate statement chunks until ';'
            if position(stmt_buf) > 0
                print(stmt_buf, " ")
            end
            print(stmt_buf, line)

            if occursin(";", line)
                # Could contain multiple ';' in one line
                s_all = String(take!(stmt_buf))
                parts = split(s_all, ';')

                # Completed statements
                for i in 1:length(parts)-1
                    stmt = strip(parts[i]) * ";"
                    isempty(strip(stmt)) && continue
                    tmp = strip_trailing_semicolon(stmt)
                    parsed = parse_frequency_statement(tmp)
                    if parsed === nothing
                        continue
                    end

                    kind, name, payload = parsed
                    if current_block === nothing
                        @warn "Found statement with no active block; skipping" statement=tmp
                        continue
                    end
                    if kind == :profile
                        payload === nothing ? (@warn "Skipping malformed profile" name=name) :
                                              (current_block.profiles[name] = payload)
                    else
                        payload === nothing ? (@warn "Skipping malformed FMIX" name=name) :
                                              (current_block.fmix[name] = payload)
                    end
                end

                # Remainder (if any)
                remainder = strip(parts[end])
                if !isempty(remainder)
                    print(stmt_buf, remainder)
                end
            end
        end
    end

    # If file ended while in models, close gracefully
    if in_models && current_block !== nothing
        @warn "File ended before 'end;' of models block. Saving partial block."
        if position(stmt_buf) > 0
            flush_statement!()
        end
        push!(blocks, current_block)
    end

    return blocks
end

# ---------------------------
# Validation / reporting
# ---------------------------
function check_profile_lengths(blocks::Vector{ModelBlock}; expected::Int=20)
    bad = 0
    total = 0
    for (bi, b) in enumerate(blocks)
        for (name, vals) in b.profiles
            total += 1
            if length(vals) != expected
                bad += 1
                @warn "Profile length mismatch" block=bi source=b.source_file name=name length=length(vals) expected=expected
            end
        end
    end
    return total, bad
end

function print_summary(blocks::Vector{ModelBlock})
    println("Parsed model blocks: ", length(blocks))
    total_profiles = sum(length(b.profiles) for b in blocks)
    total_fmix = sum(length(b.fmix) for b in blocks)
    println("Total numeric profiles : ", total_profiles)
    println("Total FMIX definitions : ", total_fmix)
    tot, bad = check_profile_lengths(blocks)
    println("Profile length check   : $(tot - bad)/$tot valid (expected length = 20)")
end

# ---------------------------
# CSV writers
# ---------------------------
function csv_escape(s::AbstractString)
    if occursin(',', s) || occursin('"', s) || occursin('\n', s)
        return "\"" * replace(s, "\"" => "\"\"") * "\""
    else
        return s
    end
end

function write_profiles_long(blocks::Vector{ModelBlock}, outpath::String)
    open(outpath, "w") do io
        println(io, "block_idx,source_file,profile_name,aa_index,value")
        for (bi, b) in enumerate(blocks)
            src = b.source_file === nothing ? "" : b.source_file
            for (pname, vals) in sort!(collect(b.profiles); by=x->x[1])
                for (i, v) in enumerate(vals)
                    @printf(io, "%d,%s,%s,%d,%.17g\n",
                        bi,
                        csv_escape(src),
                        csv_escape(pname),
                        i,
                        v
                    )
                end
            end
        end
    end
end

function write_fmix_entries(blocks::Vector{ModelBlock}, outpath::String)
    open(outpath, "w") do io
        println(io, "block_idx,source_file,fmix_name,entry_idx,component,scale,weight")
        for (bi, b) in enumerate(blocks)
            src = b.source_file === nothing ? "" : b.source_file
            for (fname, entries) in sort!(collect(b.fmix); by=x->x[1])
                for (ei, e) in enumerate(entries)
                    scale_s = e.scale === nothing ? "" : @sprintf("%.17g", e.scale)
                    weight_s = e.weight === nothing ? "" : @sprintf("%.17g", e.weight)
                    println(io,
                        string(
                            bi, ",",
                            csv_escape(src), ",",
                            csv_escape(fname), ",",
                            ei, ",",
                            csv_escape(e.component), ",",
                            scale_s, ",",
                            weight_s
                        )
                    )
                end
            end
        end
    end
end

function write_blocks_summary(blocks::Vector{ModelBlock}, outpath::String)
    open(outpath, "w") do io
        println(io, "block_idx,source_file,n_profiles,n_fmix")
        for (bi, b) in enumerate(blocks)
            src = b.source_file === nothing ? "" : b.source_file
            println(io, string(bi, ",", csv_escape(src), ",", length(b.profiles), ",", length(b.fmix)))
        end
    end
end

# ---------------------------
# CLI parsing
# ---------------------------
Base.@kwdef mutable struct CLIArgs
    model_nex::Union{Nothing,String} = nothing
    outdir::Union{Nothing,String} = nothing
    force::Bool = false
    help::Bool = false
end

function usage()
    return """
Usage:
  julia parse_udm_nex.jl --model_nex <input.nex> --outdir <output_dir> [--force]

Options:
  --model_nex <path>   Input merged NEX file (required)
  --outdir <dir>       Output directory (required)
  --force              Allow writing into non-empty output directory
  -h, --help           Show this help
"""
end

function parse_cli(args::Vector{String})
    cli = CLIArgs()
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "-h" || a == "--help"
            cli.help = true
            i += 1
        elseif a == "--force"
            cli.force = true
            i += 1
        elseif a == "--model_nex"
            i += 1
            i > length(args) && error("Missing value for --model_nex")
            cli.model_nex = args[i]
            i += 1
        elseif a == "--outdir"
            i += 1
            i > length(args) && error("Missing value for --outdir")
            cli.outdir = args[i]
            i += 1
        else
            error("Unknown argument: $a")
        end
    end
    return cli
end

function ensure_outdir(outdir::String; force::Bool=false)
    if isdir(outdir)
        entries = readdir(outdir)
        if !isempty(entries) && !force
            error("Output directory '$outdir' is not empty. Use --force to continue.")
        end
    else
        mkpath(outdir)
    end
end

# ---------------------------
# Entry point
# ---------------------------
function run_cli(args::Vector{String})
    cli = try
        parse_cli(args)
    catch e
        @error sprint(showerror, e)
        println(usage())
        return 2
    end

    if cli.help
        println(usage())
        return 0
    end

    if cli.model_nex === nothing || cli.outdir === nothing
        @error "Both --model_nex and --outdir are required."
        println(usage())
        return 2
    end

    input_path = cli.model_nex::String
    out_dir = cli.outdir::String

    if !isfile(input_path)
        @error "Input file not found" input_path
        return 2
    end

    try
        ensure_outdir(out_dir; force=cli.force)
    catch e
        @error sprint(showerror, e)
        return 2
    end

    blocks = parse_nex_models(input_path)
    print_summary(blocks)

    profiles_csv = joinpath(out_dir, "profiles_long.csv")
    fmix_csv     = joinpath(out_dir, "fmix_entries.csv")
    summary_csv  = joinpath(out_dir, "blocks_summary.csv")

    write_profiles_long(blocks, profiles_csv)
    write_fmix_entries(blocks, fmix_csv)
    write_blocks_summary(blocks, summary_csv)

    println("Wrote:")
    println("  ", profiles_csv)
    println("  ", fmix_csv)
    println("  ", summary_csv)

    return 0
end

end # module UDMNexParser

if abspath(PROGRAM_FILE) == @__FILE__
    exit(UDMNexParser.run_cli(ARGS))
end
