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
	lines = readlines(infile)

	component_freq = Dict{String, String}()
	model_fmix = Dict{String, String}()
	model_components = Dict{String, Vector{String}}()
	model_order = String[]

	for raw_line in lines
		line = strip(raw_line)

		# e.g. frequency UDM0008CLR_C0000 = ...;
		#m_comp = match(r"^frequency\s+(UDM[0-9]+L?CLR_C[0-9]+)\s*=\s*(.+);$", line)
		m_comp = match(r"^frequency\s+([^ ]+)\s*=\s*(.+);$", line)
		if m_comp != nothing && (! occursin("FMIX", line))
			component_freq[m_comp.captures[1]] = line
			continue
		end

		# e.g. frequency UDM0008CLR = FMIX{UDM0008CLR_C0000:1.0:0.3,...};
		m_mix = match(r"^frequency\s+(\S+)\s*=\s*FMIX\s*\{(.*)\}\s*;?$"i, strip(line))
		if m_mix != nothing
			model_name = m_mix.captures[1]
			fmix_body = m_mix.captures[2]
			model_fmix[model_name] = fmix_body
			push!(model_order, model_name)

			components = String[]
			for token in split(fmix_body, ",")
				comp_name = split(strip(token), ":")[1]
				push!(components, comp_name)
			end
			model_components[model_name] = components
		end
	end

	for model_name in model_order
		println(model_name)
		out_file_path = joinpath(outdir, model_name * ".dat")
		open(out_file_path, "w") do file
			for comp_name in model_components[model_name]
				if haskey(component_freq, comp_name)
					write(file, component_freq[comp_name] * "\n")
				else
					@warn "Component frequency not found for $(model_name): $(comp_name)"
				end
			end
			write(file, "model $(model_name) = POISSON+G4+FMIX{$(model_fmix[model_name])};\n")
		end
	end
end


#############################################
parse_model_nex(infile)
