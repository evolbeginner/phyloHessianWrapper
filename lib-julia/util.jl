function get_bl_node_rela()
	j = 0
	node2bls = Dict{Int,Vector{Int}}(); bl2node = Dict{Int,Int}()
	for anc in ((nb.node + nb.tip):-1:(1 + nb.tip))
		children = all_children[anc] # used to be "anc-nb.tip"
		for i in 1:length(children)
			j += 1
			push!(get!(node2bls, anc, Vector{Int}()), j)
			bl2node[j] = anc
		end
	end
	return([node2bls, bl2node])
end


function get_bls_from_treefile(script, treefile)
	cmd = `$(script) -t $(treefile)`
	run(cmd)
end


# read the last line of the file branch_out.matrix
function read_bls_from_branchout_matrix(infile)
	lines = readlines(infile)
	bls = split(lines[end], ",")
	bls = [parse(Float64, x) for x in bls]
end


function get_bl_order(infile)
	bl_order = Dict("mcmctree2ape" => Dict(), "ape2mcmctree" => Dict())
	lines = readlines(infile)
	for line in lines[2:end-1]
		line_arr = split(line, "\t")
		order1, order2 = map(x->parse(Int64,x), line_arr[(end-1):end])
		bl_order["mcmctree2ape"][order1] = order2
		bl_order["ape2mcmctree"][order2] = order1
	end
	return(bl_order)
end


function get_params_from_iqtree(infile)
	# returns
	rs = Vector{Float64}()
	props = Vector{Float64}()
	Fs = Vector{Vector{Float64}}()
	Qrs = Vector{Float64}()
	freqs = Vector{Float64}()

	# infile: iqtree_file
	is_rate_prop = false
	is_freq = false
	is_mixture = false
	empty_line_count = 0

	fh = open(infile, "r")
	lines = readlines(fh)
	for line in lines
		(line == "") ? continue : nothing #skip empty lines
		line_arr = split(line)
		# set status
		if match(r"^ Category", line) != nothing
			is_rate_prop = true
			continue
		elseif match(r"State frequencies", line) != nothing
			empty_line_count = 0
			is_freq = true
			push!(Fs, Vector{Float64}())
			continue
		elseif match(r"Mixture model of substitution", line) != nothing
			empty_line_count = 0
			is_mixture = true
			continue
		end

		# push
		if is_rate_prop
			if match(r"^  [0-9]", line) == nothing
				is_rate_prop = false
				continue
			end
			push!(rs, parse(Float64, line_arr[2]))
			push!(props, parse(Float64, line_arr[3]))
		elseif is_freq
			#pi(A) = 0.0600
			if match(r"^.+pi\([A-Z]\)", line) != nothing
				push!(Fs[end], parse(Float64, line_arr[3]))
			else
				is_freq = false
				continue
			end
		elseif is_mixture
   			#  1  BurEX2        0.5269   0.3775   BurEX2
			if match(r"^\s+[0-9]", line) != nothing
				push!(Qrs, parse(Float64, line_arr[3]))
				push!(freqs, parse(Float64, line_arr[4]))
			elseif match(r"^\s+No", line) != nothing
				continue
			else
				is_mixture = false
				continue
			end
		end
	end
	close(fh)

	return(Dict(:rs=>rs, :props=>props, :Fs=>Fs, :Qrs=>Qrs, :freqs=>freqs))
end


function all_rows_identical(matrix::AbstractMatrix)
	# Check if all rows are the same as the first row
	first_row = matrix[1, :]
	all(row -> row == first_row, eachrow(matrix))
end


function format_number(x)
    rounded_number = round(x, sigdigits=4)
    if abs(x) >= 1e4 || abs(x) < 1e-3
        return rounded_number
    else
        return rounded_number
    end
end


function is_complex_matrix(A)
           return eltype(A) <: Complex
end


function find_nonzero_imaginary_parts(A)
    return filter(x -> isa(x,Complex) && imag(x) != 0, A)
end


function convert_complex_to_float_in_matrix(data)
    # Determine if the input is a vector and reshape it to a matrix if necessary
    if ndims(data) == 1
        data = reshape(data, :, 1)  # Reshape vector to a column matrix
    end

    # Get the dimensions of the data
    rows, cols = size(data)
    # Initialize a new matrix of Float64
    new_data = Matrix{Float64}(undef, rows, cols)
    # Iterate over each element
    for i in 1:rows
        for j in 1:cols
            # Assign only the real part as Float64
            new_data[i, j] = real(data[i, j])
        end
    end

    # If the original data was a vector, convert it back to a vector
    if cols == 1
        return vec(new_data)  # Convert single-column matrix back to vector
    else
        return new_data
    end
end


###########################################################
function getCorename(infile, is_strict=false)
  b = basename(infile)
  if !is_strict
    if occursin(".", b)
      m = match(r"(.+)(\.[^.]+)+$", b)
    else
      m = match(r"(.+)", b)
    end 
  else
    m = match(r"^([^.]+)", b)
  end 

  if m != nothing
    c = m[1]
  else
    c = ""
  end
  return(c)
end


function mkdir_with_force(outdir::String, is_force::Bool=false, is_tolerate::Bool=false)
    if !isdir(outdir)
        run(`mkdir -p $outdir`)
    else
        if is_tolerate
            # Do nothing 
        elseif is_force
            run(`rm -rf $outdir`)
            run(`mkdir -p $outdir`)
        else
            error("The outdir $outdir has already existed!")
        end
    end
end


