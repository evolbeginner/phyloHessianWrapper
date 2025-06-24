#! /usr/bin/env julia


#############################################
include(joinpath(ENV["MY_JULIA_PATH"], "Dir.jl"))

using ArgParse


#############################################
model_nex = nothing
outdir = nothing
is_force = false

opt = ArgParseSettings()
@add_arg_table opt begin
	"--model_nex"
		help = ""
		arg_type = String
		default = nothing
	"--outdir"
		help = "outdir"
		arg_type = String
		default = nothing
	"--force"
		action = :store_true
end
parsed_opt = parse_args(opt)

infile = expanduser(parsed_opt["model_nex"])
outdir = expanduser(parsed_opt["outdir"])
is_force = parsed_opt["force"]

if outdir != nothing
	mkdir_with_force(outdir, is_force)
else
	error("outdir not specified.")
end


#############################################
function parse_model_nex(infile)
	is_prepare_model_name = false
	is_new_model = false
	is_within_model = false
	model_contents = Vector{String}()

	in_fh = open(infile, "r")
	for line in readlines(in_fh)
		if is_new_model
			#println(line)
			is_new_model = false
		end
		if occursin("[ " * repeat("-",10), line)
			is_new_model = true
			is_within_model = true
		end

	#[ main definition of EX2 with fixed component rates ]
		if occursin(r"^\[ ", line)
			is_prepare_model_name = true
		end

		if is_within_model
			line = replace(line, r"[ ]+;$" => ";")
			line = replace(line, r"^[ ]([0-9]+\.[0-9])" => s"\1")
			push!(model_contents, line)
		end

		m = match(r"^model ([^=]+)\s+?=.+;", line)
		if m != nothing && is_prepare_model_name
			model_name = m.captures[1]
			println(model_name)

			out_file_path = joinpath(outdir, model_name * ".dat")
			open(out_file_path, "w") do file
                write(file, join(model_contents, "\n"))
            end
            empty!(model_contents)

			#println(join(model_contents,"\n"))
			#empty!(model_contents)
			is_prepare_model_name = false
			is_within_model = false
			println(repeat("\n", 3))
		end

	end
	close(in_fh)
end


#############################################
#=
[ ---------------------------------------------------------
    EX2 mixture model of Le, Lartillot & Gascuel (2008) 
 --------------------------------------------------------- ]
.
.
.
[ main definition of EX2 with fixed component rates ]
model EX2 =MIX{BurEX2:0.672020808818762,ExpEX2:1.6413466609931};
=#


#infile = expanduser("~/software/phylo/iqtree/latest/models.nex")

parse_model_nex(infile)


