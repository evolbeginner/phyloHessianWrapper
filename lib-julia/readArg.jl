#! /bin/env julia


#####################################
dir = dirname(@__FILE__)
include(joinpath(dir, "util.jl"))

using ArgParse


#####################################
type, indir = nothing, nothing
iqtree_file = nothing
basics_indir = nothing
branchout_matrix = nothing
bs_branchout_matrix = nothing

hessian_outfile = nothing
hessian_type = "SKT2004"
hessian_infile = nothing
is_compare = false

pmsf_file = nothing
is_pmsf = false


#####################################
function parse_commandline()
    opt = ArgParseSettings()

    @add_arg_table opt begin
		"--type", "-t"
			help = "type of seq"
			arg_type = String
			default = "DNA"
        "--indir"
			help = "indir"
            arg_type = String
            default = nothing
		"--basics_indir"
			help = "basics_indir"
			arg_type = String
			default = nothing
		"--tree", "--te"
			help = "treefile"
			arg_type = String
			default = nothing
		"--iqtree"
			help = ".iqtree"
			arg_type = String
			default = nothing
		"--phyml"
			help = ".phy_phymls.txt"
			arg_type = String
			default = nothing
		"--branchout_matrix", "-b"
			help = "branch_out.matrix"
			arg_type =  String
			default = nothing
		"--bs_branchout_matrix", "--b2"
			help = "bs_branch_out.matrix if any"
			arg_type =  String
			default = nothing
		"--model", "-m"
			help = "substitution matrix"
			arg_type =  String
			default = "LG"
		"--mix_freq_model", "--mix_freq", "--mfm"
			help = "site heterogeneous equilibrium freq mixture model (e.g., C10-C60, CF4)"
			arg_type = String
			default = nothing
		"--pmsf"
			help = "PMSF file"
			arg_type = String
			default = nothing
		"--outdir"
			help = "outdir"
			arg_type = String
			default = nothing
		"--force"
			help = "force"
			action = :store_true
		"--tolerate"
			help = "tolerate"
			action = :store_true
		"--hessian", "--inBV", "--in_BV"
			help = "hessian file in the mcmctree order"
			arg_type = String
			default = nothing
		"--hessian_type"
			help = "SKT2004 (default) or fd (finite_difference)"
			arg_type = String
			default = "SKT2004"
		"--read_hessian", "--read_inBV", "--read_in_BV"
			help = "read the in.BV file"
			arg_type = String
			default = nothing
		"--compare"
			help = "to compare btwn extact and approx lnL"
			action = :store_true
		"--transform"
			help = "transform method"
			arg_type = String
			default = nothing
    end

    return(parse_args(opt))
end


function check_input(opt)
	if ! haskey(opt, "type")
		error("type not given")
	end
	if ! (opt["type"] in ["AA", "DNA"])
		error("type has to be AA or DNA")
	end
end



#####################################
opt = parse_commandline()

check_input(opt)

type = opt["type"]
indir = opt["indir"]
basics_indir = opt["basics_indir"]
treefile = opt["tree"]
iqtree_file = opt["iqtree"]
phyml_file = opt["phyml"]
branchout_matrix = opt["branchout_matrix"]
bs_branchout_matrix = opt["bs_branchout_matrix"]
sub_model = opt["model"]

outdir = opt["outdir"]
is_force = opt["force"]
is_tolerate = opt["tolerate"]
#hessian_outfile = opt["hessian"]
hessian_type = opt["hessian_type"]
hessian_infile = opt["read_hessian"]

is_compare = opt["compare"]
if hessian_infile != nothing
	is_compare = true
end


if outdir == nothing
	error("outdir not specified. exiting ......")
else
	mkdir_with_force(outdir, is_force, is_tolerate)
	std_outfh = open(joinpath(outdir, "info"), "w")
	hessian_outfile = joinpath(outdir, "in.BV")
end

if basics_indir === nothing
	basics_indir = "julia"
end

if ! occursin(r"phy_phyml_tree\.txt$", treefile)
	if iqtree_file === nothing
		iqtree_file = joinpath( dirname(treefile), getCorename(treefile)*".iqtree" )
	end
else
	if phyml_file === nothing
		phyml_file = replace(treefile, "phy_phyml_tree" => "phy_phyml_stats")
	end
end

if branchout_matrix === nothing
	branchout_matrix = "branch_out.matrix"
end

mix_freq_model = opt["mix_freq_model"]
mix_freq_model = (mix_freq_model == "nothing") ? nothing : mix_freq_model

pmsf_file = opt["pmsf"]
is_pmsf = (pmsf_file != nothing) ? true : false

transform_method = opt["transform"]


