#! /usr/bin/env julia


######################################################
using ArgParse
using DelimitedFiles

include(joinpath(ENV["MY_JULIA_PATH"], "Dir.jl"))
include("lib/retrieve_model_text.jl")


######################################################
infile = nothing
type = nothing
outdir = nothing


opt = ArgParseSettings()
@add_arg_table opt begin
	"--infile", "--in", "-i"
		help = "infile"
		arg_type = String
		default = nothing
	"--type", "-t"
		help = "AA or DNA"
		arg_type = String
		default = nothing
	"--outdir"
		help = "outdir"
		arg_type = String
		default = nothing
	"--force"
		help = "is_force"
		action = :store_true
	"--tolerate"
		help = "is_tolerate"
		action = :store_true
end

parsed_opt = parse_args(opt)

infile = parsed_opt["infile"]
type = parsed_opt["type"]
outdir = parsed_opt["outdir"]
is_force = parsed_opt["force"]
is_tolerate = parsed_opt["tolerate"]

mkdir_with_force(outdir, is_force, is_tolerate)


######################################################
if type == "AA"
	n_char = 20
elseif type == "DNA"
	n_char = 4
else
	throw("wrong type $type")
	exit()
end


######################################################
r = fill("Inf", n_char-1, n_char-1)
freqs = fill(0.0, n_char)

line_dict = filter_lines(extract_lines(infile))


######################################################
for (k,lines) in line_dict
	#daa[ 1*20+ 0] =   58.00; daa[ 2*20+ 0] =   54.00; daa[ 2*20+ 1] =   45.00; daa[ 3*20+ 0] =   81.00;
	pattern1 = r"daa\[\s*(\d+)\*(\d+)\+\s*(\d+)\]\s*=\s*([^;]+);"
	#f[ 0] = 0.076748; f[ 1] = 0.051691; f[ 2] = 0.042645; f[ 3] = 0.051544;
	pattern2 = r"f\[\s*(\d+)\]\s*=\s*([^;]+);"

	for line in lines
		matches = eachmatch(pattern1, line)
		for m in matches
			row, col = map(x->parse(Int64, x), m.captures[[1,3]])
			col = col + 1
			#r[row, col] = parse(Float64, m.captures[4])
			r[row, col] = m.captures[4]
		end

		matches = eachmatch(pattern2, line)
		for m in matches
			index = parse(Int64, m.captures[1]); index = index + 1
			freq = parse(Float64, m.captures[2])
			freqs[index] = freq
		end
	end

	r[r .== "Inf"] .= ""

	outfile = joinpath(outdir, k * ".dat")
	out_fh = open(outfile, "w")
	writedlm(out_fh, r)
	println(out_fh, "")
	println(out_fh, join(freqs, " "))
	close(out_fh)

	run(`sed -i 's/\t\+/ /g' $outfile`)
	run(`sed -i 's/ \+$/ /g' $outfile`)
	run(`sed -i 's/[ ]\+/ /g' $outfile`)
end


