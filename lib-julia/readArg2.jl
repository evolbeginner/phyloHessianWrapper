#! /bin/env julia

#####################################
dir = dirname(@__FILE__)
include(joinpath(dir, "util.jl"))

#####################################
# Default values
defaults = Dict(
    "type" => "DNA",
    "model" => "LG",
    "hessian_type" => "SKT2004"
)

flags = Set(["force", "tolerate", "transform"])

aliases = Dict(
    "t" => "type",
    "type" => "type",
    "indir" => "indir",
    "basics_indir" => "basics_indir",
    "tree" => "tree",
    "te" => "tree",
    "iqtree" => "iqtree",
    "phyml" => "phyml",
    "b" => "branchout_matrix",
    "branchout_matrix" => "branchout_matrix",
    "m" => "model",
    "model" => "model",
    "mix_freq_model" => "mix_freq_model",
    "mix_freq" => "mix_freq_model",
    "mfm" => "mix_freq_model",
    "pmsf" => "pmsf",
    "outdir" => "outdir",
    "force" => "force",
    "tolerate" => "tolerate",
    "hessian" => "read_hessian",
    "inBV" => "read_hessian",
    "in_BV" => "read_hessian",
    "hessian_type" => "hessian_type",
    "read_hessian" => "read_hessian",
    "read_inBV" => "read_hessian",
    "read_in_BV" => "read_hessian",
    "transform" => "transform"
)

function parse_args(args)
    opts = Dict{String, Any}()

    i = 1
    while i <= length(args)
        arg = args[i]
        if startswith(arg, "--")
            key = replace(arg, "--" => "")
        elseif startswith(arg, "-")
            key = replace(arg, "-" => "")
        else
            i += 1
            continue
        end

		if ! haskey(aliases, key)
			error(arg * "is not a valid argument")
		end

        canonical = get(aliases, key, key)

        if canonical in flags
            opts[canonical] = true
            i += 1
        elseif i < length(args) && !startswith(args[i+1], "-")
            opts[canonical] = args[i+1]
            i += 2
        else
            opts[canonical] = true
            i += 1
        end
    end

    return opts
end

function error_msg(sent)
    println(stderr, sent)
    exit(1)
end

#####################################
opt = parse_args(ARGS)

# Merge defaults
for (k, v) in defaults
    if !haskey(opt, k)
        opt[k] = v
    end
end

# Validate input
if !(opt["type"] in ["AA", "DNA"])
    error_msg("type has to be AA or DNA")
end

# Assign variables
type = opt["type"]
indir = get(opt, "indir", nothing)
basics_indir = get(opt, "basics_indir", "julia")
treefile = get(opt, "tree", nothing)
iqtree_file = get(opt, "iqtree", nothing)
phyml_file = get(opt, "phyml", nothing)
branchout_matrix = get(opt, "branchout_matrix", "branch_out.matrix")
sub_model = opt["model"]
mix_freq_model = get(opt, "mix_freq_model", nothing)
pmsf_file = get(opt, "pmsf", nothing)
outdir = get(opt, "outdir", nothing)
is_force = get(opt, "force", false)
is_tolerate = get(opt, "tolerate", false)
hessian_type = opt["hessian_type"]
hessian_infile = get(opt, "read_hessian", nothing)
transform_method = get(opt, "transform", false)

# Output setup
if outdir == nothing
    error("outdir not specified. exiting ......")
else
    mkdir_with_force(outdir, is_force, is_tolerate)
    std_outfh = open(joinpath(outdir, "info"), "w")
    hessian_outfile = joinpath(outdir, "in.BV")
end

# Tree file logic
if treefile !== nothing && !occursin(r"phy_phyml_tree\.txt$", treefile)
    if iqtree_file === nothing
        iqtree_file = joinpath(dirname(treefile), getCorename(treefile) * ".iqtree")
    end
elseif treefile !== nothing
    if phyml_file === nothing
        phyml_file = replace(treefile, "phy_phyml_tree" => "phy_phyml_stats")
    end
end

is_pmsf = pmsf_file !== nothing
mix_freq_model = mix_freq_model == "nothing" ? nothing : mix_freq_model

