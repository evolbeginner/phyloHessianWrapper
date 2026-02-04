module PhyloParser

# Load dependencies
using ArgParse

# Include utility functions
dir = @__DIR__
include(joinpath(dir, "util.jl"))

# Exported functions
export parse_commandline, check_input, error_msg, run_parser

# -----------------------------
# Function Definitions
# -----------------------------

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
        "--model", "-m"
            help = "substitution matrix"
            arg_type =  String
            default = "LG"
        "--mix_freq_model", "--mix_freq", "--mfm"
            help = "site heterogeneous equilibrium freq mixture model"
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
            help = "hessian file"
            arg_type = String
            default = nothing
        "--hessian_type"
            help = "SKT2004 or fd"
            arg_type = String
            default = "SKT2004"
        "--read_hessian", "--read_inBV", "--read_in_BV"
            help = "read in.BV file"
            arg_type = String
            default = nothing
        "--transform"
            help = "transform method"
            action = :store_true
    end

    return parse_args(opt)
end

function check_input(opt)
    if !haskey(opt, "type")
        error_msg("type not given")
    end
    if !(opt["type"] in ["AA", "DNA"])
        error_msg("type has to be AA or DNA")
    end
end

function error_msg(sent)
    println(stderr, sent)
    exit(1)
end

function run_parser()
    opt = parse_commandline()
    check_input(opt)

    # You can move the rest of your logic here or into separate functions
    # For example: setup_output_dir(opt), resolve_tree_files(opt), etc.
end

end # module
