#!/usr/bin/env julia

using Dates
println(now(), " Starting ......")

using ArgParse
using FASTX
using Phylo
println(now())

const AA_INDEX = Dict(c => i for (i, c) in enumerate(collect("ARNDCQEGHILKMFPSTWYV")))
const DNA_INDEX = Dict('A' => 1, 'C' => 2, 'G' => 3, 'T' => 4)

mutable struct SimpleNode
    name::String
    length::Float64
    children::Vector{SimpleNode}
    id::Int
end

function read_tree(path::String)
    tree = open(parsenewick, path)
    tip_order = String.(getleafnames(tree))
    tip_rank = Dict(name => i for (i, name) in enumerate(tip_order))
    phylo_leaves(node) = isleaf(tree, node) ? [String(getnodename(tree, node))] :
        reduce(vcat, phylo_leaves.(getchildren(tree, node)))
    function copy_node(node)
        raw_len = isroot(tree, node) ? 0.0 : getlength(tree, getinbound(tree, node))
        len = ismissing(raw_len) ? 0.0 : Float64(raw_len)
        raw_kids = collect(getchildren(tree, node))
        sort!(raw_kids, by=child -> minimum(tip_rank[name] for name in phylo_leaves(child)))
        kids = [copy_node(child) for child in raw_kids]
        SimpleNode(String(getnodename(tree, node)), len, kids, 0)
    end
    copy_node(getroot(tree))
end

function read_rooted_species_tree(path::String)
    # MCMCtree species.trees files commonly start with a metadata line such as
    # "5 1" and may end the Newick with calibration text such as >8<12;.
    for raw in eachline(path)
        line = strip(raw)
        (occursin('(', line) && occursin(')', line)) || continue
        line = replace(line, r">[^;]*" => "")
        endswith(line, ';') || (line *= ";")
        tmp = tempname()
        try
            write(tmp, line * "\n")
            return read_tree(tmp)
        finally
            isfile(tmp) && rm(tmp)
        end
    end
    error("Could not find a rooted Newick tree in $path")
end

leaves(node::SimpleNode) = isempty(node.children) ? [node.name] : reduce(vcat, leaves.(node.children))

function unroot_topology!(root::SimpleNode)
    if length(root.children) == 2
        split_i = findfirst(child -> !isempty(child.children), root.children)
        split_i === nothing && error("Cannot unroot a two-tip tree")
        split = root.children[split_i]
        length(split.children) == 2 || error("The selected root child is not binary")
        other = root.children[3 - split_i]
        # The suppressed root edge is immaterial for a topology-only reference.
        root = SimpleNode("root", 0.0, [split.children[1], split.children[2], other], 0)
    end
    length(root.children) == 3 || error("An unrooted reference must have three root children")
    # Match paml_order_unroot.R: retain child order within subtrees, but put
    # internal root children before singleton tips. This makes tree_order in
    # branch_out.matrix compatible with the historical R/Ruby pipeline.
    sort!(root.children, by=child -> isempty(child.children) ? 1 : 0, alg=Base.Sort.MergeSort)
    root
end

function write_newick(path::String, root::SimpleNode; lengths::Bool=false)
    function emit(n)
        label = isempty(n.children) ? n.name : "(" * join(emit.(n.children), ", ") * ")"
        lengths && n !== root ? label * ":" * string(n.length) : label
    end
    open(path, "w") do io
        println(io, emit(root), ";")
    end
end

function read_fasta(path::String)
    seqs = Dict{String, String}()
    open(FASTA.Reader, path) do reader
        for record in reader
            name = String(identifier(record))
            isempty(name) && error("Empty FASTA identifier")
            haskey(seqs, name) && error("Duplicate FASTA sequence name: $name")
            seqs[name] = uppercase(sequence(String, record))
        end
    end
    isempty(seqs) && error("No sequences found in $path")
    lengths = unique(length.(values(seqs)))
    length(lengths) == 1 || error("FASTA sequences have different lengths")
    seqs
end

function assign_ids!(root::SimpleNode, tip_order::Vector{String})
    tip_id = Dict(name => i for (i, name) in enumerate(tip_order))
    internals = SimpleNode[]
    function visit(n)
        if isempty(n.children)
            n.id = get(tip_id, n.name, 0)
            n.id == 0 && error("Tree tip $(n.name) is absent from the tip order")
        else
            push!(internals, n)
            visit.(n.children)
        end
    end
    visit(root)
    for (i, n) in enumerate(internals)
        n.id = length(tip_order) + i
    end
    internals
end

function descendant_ids(n::SimpleNode; tips_only::Bool=false)
    if !tips_only
        out = Int[]
        queue = copy(n.children)
        while !isempty(queue)
            child = popfirst!(queue)
            push!(out, child.id)
            append!(queue, child.children)
        end
        return out
    end
    out = Int[]
    function visit(x)
        for child in x.children
            (!tips_only || isempty(child.children)) && push!(out, child.id)
            visit(child)
        end
    end
    visit(n)
    out
end

function encode_patterns(seqs::Dict{String,String}, tips::Vector{String}, seqtype::String)
    missing = setdiff(tips, collect(keys(seqs)))
    extra = setdiff(collect(keys(seqs)), tips)
    isempty(missing) || error("Alignment is missing tree tips: $(join(missing, ", "))")
    isempty(extra) || error("Alignment has tips absent from tree: $(join(extra, ", "))")
    nsite = length(seqs[first(tips)])
    alphabet = seqtype == "AA" ? AA_INDEX : DNA_INDEX
    patterns = Vector{Vector{Int}}()
    counts = Int[]
    site2pattern = Vector{Int}(undef, nsite)
    lookup = Dict{Tuple{Vararg{Int}},Int}()
    for site in 1:nsite
        states = [get(alphabet, seqs[name][site], 999) for name in tips]
        key = Tuple(states)
        pi = get(lookup, key, 0)
        if pi == 0
            push!(patterns, states); push!(counts, 0)
            pi = length(patterns); lookup[key] = pi
        end
        counts[pi] += 1
        site2pattern[site] = pi
    end
    patterns, counts, site2pattern
end

split_key(side, all_tips) = begin
    a = sort!(collect(side), by=lowercase)
    b = sort!(setdiff(all_tips, a), by=lowercase)
    sa, sb = join(a, "-"), join(b, "-")
    lowercase(sa) <= lowercase(sb) ? (sa * "," * sb) : (sb * "," * sa)
end

display_split(side, all_tips) = join(sort!(collect(side), by=lowercase), "-") * "," *
                                join(sort!(setdiff(all_tips, side), by=lowercase), "-")

function edge_records(root::SimpleNode, all_tips::Vector{String}; julia_order::Bool=false)
    records = NamedTuple[]
    if julia_order
        nodes = SimpleNode[]
        function collect_internal(n)
            isempty(n.children) || push!(nodes, n)
            collect_internal.(n.children)
        end
        collect_internal(root)
        for parent in reverse(nodes), child in parent.children
            side = leaves(child)
            push!(records, (key=split_key(side, all_tips), display=display_split(side, all_tips),
                            length=child.length, parent=parent.id, child=child.id))
        end
    else
        function visit(parent)
            for child in parent.children
                side = leaves(child)
                push!(records, (key=split_key(side, all_tips), display=display_split(side, all_tips),
                                length=child.length, parent=parent.id, child=child.id))
                isempty(child.children) || visit(child)
            end
        end
        visit(root)
    end
    records
end

function write_basics(outdir, root, tips, seqs, seqtype)
    mkpath(outdir)
    internals = assign_ids!(root, tips)
    patterns, counts, site2pattern = encode_patterns(seqs, tips, seqtype)
    open(joinpath(outdir, "basics"), "w") do io
        println(io, "nb.node\t", length(internals)); println(io, "nb.tip\t", length(tips))
    end
    open(joinpath(outdir, "all_children"), "w") do io
        for n in internals println(io, join([n.id; getfield.(n.children, :id)], '\t')) end
    end
    open(joinpath(outdir, "descendants"), "w") do io
        for n in internals println(io, join([n.id; descendant_ids(n; tips_only=true)], '\t')) end
    end
    open(joinpath(outdir, "all"), "w") do io
        for n in internals println(io, join([n.id; descendant_ids(n)], '\t')) end
    end
    open(joinpath(outdir, "cherry"), "w") do io
        for n in internals all(isempty(c.children) for c in n.children) && println(io, n.id) end
    end
    open(joinpath(outdir, "pattern"), "w") do io
        for i in eachindex(patterns) println(io, join([patterns[i]; counts[i]], '\t')) end
    end
    open(joinpath(outdir, "site2pattern"), "w") do io
        for i in eachindex(site2pattern) println(io, i, '\t', site2pattern[i]) end
    end
end

function write_branches(outdir, fitted, reference, tips)
    mkpath(outdir)
    assign_ids!(fitted, tips); assign_ids!(reference, leaves(reference))
    fit_records = edge_records(fitted, tips; julia_order=true)
    fit_by_split = Dict(r.key => r.length for r in fit_records)
    ref_records = edge_records(reference, tips)
    ref_keys = [r.key for r in ref_records]
    Set(keys(fit_by_split)) == Set(ref_keys) || error("Fitted and reference tree topologies differ")
    open(joinpath(outdir, "branch_out"), "w") do io
        for r in ref_records println(io, r.display, '\t', fit_by_split[r.key]) end
    end
    tree_order = Dict(key => i for (i, key) in enumerate(ref_keys))
    ape_order = Dict(r.key => i for (i, r) in enumerate(fit_records))
    # Match R's locale-aware ordering used by align_bl.R: splits beginning
    # with the same tip place the multi-tip side before the singleton side.
    sorted_records = sort(ref_records, by=r -> begin
        side = split(first(split(r.display, ',')), '-')
        (lowercase(side[1]), -length(side), lowercase(join(side, '-')))
    end)
    open(joinpath(outdir, "branch_out.matrix"), "w") do io
        println(io, "branch\tbl\ttree_order\tape_order")
        for r in sorted_records
            println(io, tree_order[r.key], '\t', r.display, '\t', fit_by_split[r.key], '\t', tree_order[r.key], '\t', ape_order[r.key])
        end
        println(io, join(getfield.(fit_records, :length), ','))
    end
end

function parse_cli()
    settings = ArgParseSettings(description="Prepare julia_bl.jl tree, pattern, and branch-order inputs using Phylo.jl")
    @add_arg_table! settings begin
        "--seq", "-s"; required=true; arg_type=String; help="FASTA alignment"
        "--tree", "-t"; required=true; arg_type=String; help="fitted Newick tree with branch lengths"
        "--species-tree"; required=true; arg_type=String; help="rooted MCMCtree species tree used to generate ref.tre"
        "--outdir", "-o"; required=true; arg_type=String; help="output directory containing julia/, bl/, and ref.tre"
        "--type"; arg_type=String; default="AA"; help="AA or DNA"
        "--force"; action=:store_true; help="replace julia/ and bl/ output directories"
    end
    parse_args(settings)
end

function main()
	println(now())
    args = parse_cli()
    seqtype = uppercase(args["type"])
    seqtype in ("AA", "DNA") || error("--type must be AA or DNA")
    outdir = args["outdir"]
    julia_dir, bl_dir = joinpath(outdir, "julia"), joinpath(outdir, "bl")
    if !args["force"] && (ispath(julia_dir) || ispath(bl_dir))
        error("Output julia/ or bl/ already exists; use --force")
    end
	println(now())
    args["force"] && foreach(path -> ispath(path) && rm(path; recursive=true), (julia_dir, bl_dir))
    fitted = read_tree(args["tree"])
    # Equivalent to: Rscript additional_scripts/paml_order_unroot.R
    # species.trees ref.tre. The fitted tree supplies branch lengths only.
    reference = unroot_topology!(read_rooted_species_tree(args["species-tree"]))
    fit_tips, ref_tips = leaves(fitted), leaves(reference)
    Set(fit_tips) == Set(ref_tips) || error("Fitted and reference trees have different tip sets")
    ref_path = joinpath(outdir, "ref.tre")
    mkpath(outdir); write_newick(ref_path, reference)
    seqs = read_fasta(args["seq"])
	println(now())
    write_basics(julia_dir, fitted, fit_tips, seqs, seqtype)
    write_branches(bl_dir, fitted, reference, fit_tips)
    println("Wrote ", ref_path)
    println("Wrote ", julia_dir)
    println("Wrote ", joinpath(bl_dir, "branch_out.matrix"))
	println(now(), " Done!")
end

main()
