function extract_lines(filename::String)
    lines = readlines(filename)  # Read all lines from the file
    result = Dict()  # Initialize an empty vector to store the lines
    inside_block = false  # Flag to track if we're inside the block
	model = nothing

    for line in lines
        m = match(r"case PLL_(\w+)", line)
		if m != nothing
			model = m.captures[1]
			result[model] = String[]
            inside_block = true
        end

        if inside_block
            push!(result[model], line)
        end

        if inside_block && occursin("break;", line)
			inside_block = false
        end
    end

    return result
end


function filter_lines(result)
	select_correct_model(result)
end


function select_correct_model(result)
	# correct model: e.g., WAG, JTT, FLU
	line_dict = Dict{String, Vector{String}}()
	for (k,v) in result
		length(v) <= 1 && continue
		if ! any(x->occursin(r"daa\[\d", x), v)
			continue
		end
		line_dict[k] = v
	end
	return(line_dict)
end

#=
filename = ARGS[1]
extracted_lines = extract_lines(filename)

for (k,v) in extracted_lines
	println(k)
	map(println, v)
	println()
end
=#


