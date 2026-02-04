using Distributions, QuadGK


###########################################################
function get_iqtree_params(iqtree_file, phyml_file)
	if iqtree_file != nothing
		iqtree_params = get_params_from_iqtree(iqtree_file)
	elseif phyml_file != nothing
		iqtree_params = get_params_from_phyml(phyml_file)
	else
		error("Either .iqtree or .phy_phyml_stats.txt should be provided. exiting ......")
	end

	rs = iqtree_params[:rs]
	props = iqtree_params[:props]
	Fs = iqtree_params[:Fs]
	Qrs = iqtree_params[:Qrs]
	freqs = iqtree_params[:freqs]
	inv_info = Dict(:is_do_inv=>false, :inv_prop=>0.0)

	if phyml_file != nothing
		Qrs .= Qrs / sum(freqs .* Qrs)
	end

	if isempty(rs)
		rs = [1.0]
		props = [1.0]
	elseif any(x->x==0, rs)
		inv_info = Dict(:is_do_inv=>true, :inv_prop=>props[1])
		rs = iqtree_params[:rs][2:end]
		props = iqtree_params[:props][2:end]
		if phyml_file != nothing
			props = props ./ sum(iqtree_params[:props])
			rs = rs ./ sum(props)
		end
	end

	println()
	println("rs:\t", rs)
	println("props:\t", props)
	println("freqs:\t", freqs)
	println("Qrs:\t", Qrs)
	return([rs, props, Fs, Qrs, freqs, inv_info])
end


function read_basics(basics_indir)
	# Get "all_children"
	infile = joinpath(basics_indir, "all_children")
	all_children = Dict{Int, Vector{Int}}()
	for line in eachline(infile)
		parts = split(line, '\t')
		key = parse(Int, parts[1])
		values = [parse(Int, part) for part in parts[2:end]]
		all_children[key] = values
	end
	println(std_outfh, "all_children")
	println(std_outfh, all_children)
	println(std_outfh, "\n")

	# Get the pattern or "v"
	infile = joinpath(basics_indir, "pattern")
	num_patterns = countlines(infile)
	#pattern = []
	pattern0 = Vector{Tuple{Vector{Int}, Int}}(undef, num_patterns)
	for (i, line) in enumerate(eachline(infile))
		parts = split(line, '\t')
		key = parse(Int, parts[end])
		values = [parse(Int, part) for part in parts[1:end-1]]
		#push!(pattern, [values, key])
		pattern0[i] = (values, key)
	end
	pattern = pattern0

	# Get the basic info of the tree
	infile = joinpath(basics_indir, "basics")
	basics = Dict{String,Int}()
	for line in eachline(infile)
		parts = split(line, '\t')
		basics[parts[1]] = parse(Int, parts[2])
	end
	println(std_outfh, "basics")
	println(std_outfh, basics)
	println(std_outfh, "\n")
	flush(std_outfh)

	nb = Nb(basics["nb.node"], basics["nb.tip"], basics["nb.node"]+basics["nb.tip"]-1)

	# Get cherry nodes
	infile = joinpath(basics_indir, "cherry")
	cherry_nodes = Vector{Int}()
	for line in eachline(infile)
		push!(cherry_nodes, parse(Int, line))
	end

	# Get descendants
	infile = joinpath(basics_indir, "descendants")
	#descendants = Dict{Int,Dict}()
	descendants = Dict()
	for line in eachline(infile)
		line_arr = split(line, '\t')
		int_node = parse(Int, line_arr[1])
		descendants[int_node] = [parse(Int, x) for x in line_arr[2:end]]
	end


	infile = joinpath(basics_indir, "site2pattern")
	site2pattern = Dict(); pattern2site = Dict()
	for line in eachline(infile)
		line_arr = map(x->parse(Int,x), split(line,'\t'))
		site2pattern[line_arr[1]] = line_arr[2]
		pattern2site[line_arr[2]] = line_arr[1]
	end

	return(nb, pattern, all_children, cherry_nodes, descendants, site2pattern)
end


function read_hessian_infile(hessian_infile)
	in_fh = open(hessian_infile, "r")
	is_read_hessian = false
	count = 0
	hessian_data = Float64[]

	for line in readlines(in_fh)
		if is_read_hessian
			if ! occursin(r"\w", line)
				continue
			else
				nums = parse.(Float64, split(line, r"\s+"))
				append!(hessian_data, nums)
			end
		end
		if occursin("Hessian", line)
			is_read_hessian = true
		end
	end
	close(in_fh)

	n = Int(sqrt(length(hessian_data)))
    if n^2 != length(hessian_data)
    	@warn "Hessian's length is not a square of a integer!"
        return reshape(hessian_data, :, length(hessian_data) ÷ n)
    end
	h = reshape(hessian_data, n, n)

	return(h)
end


###########################################################
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
	#bls = split(lines[end], ",")
	#bls = [parse(Float64, x) for x in bls]
	lines_with_bls = filter(line -> ! occursin(r"\s", line), lines)
	bls_vec = map(x -> [parse(Float64, i) for i in split(x, ",")], lines_with_bls)
end


function get_bl_order(infile)
	bl_order = Dict("mcmctree2ape" => Dict(), "ape2mcmctree" => Dict())
	lines = readlines(infile)
	lines_with_tab = filter(line -> occursin(r"\t", line), lines[2:end])
	for line in lines_with_tab
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


function get_params_from_phyml(infile)
	rs = Vector{Float64}()
	props = Vector{Float64}()
	Fs = Vector{Vector{Float64}}()
	Qrs = Vector{Float64}()
	freqs = Vector{Float64}()

	# Gamma discrete
	ncat = 1
	alpha = nothing

	# infile: 1.phy_phyml_stats.txt
	is_rate_prop = false
	is_freq = false
	is_mixture = false

	fh = open(infile, "r")
	lines = readlines(fh)
	for line in lines
		line_arr = split(line, r"\s+")
		if occursin("Discrete gamma model", line)
			is_gamma = true
		end
		if occursin("FreeRate model", line)
			is_freerate = true
		end

		m = match(r"Number of (?:classes|categories)[ ]?:\s+(\d+)", line)
		if m != nothing
			ncat = parse(Int32, m[1])
		end

		m = match(r"Gamma shape parameter[ ]?:\s(\S+)", line)
		if m != nothing
			alpha = parse(Float64, m[1])
		end

		m = match(r"Relative rate in class (?:\d+)[ ]?:\s+(\d+\.\d+) \[freq=(.+)\]", line)
		if m != nothing
			push!(rs, parse(Float64, m[1]))
			push!(props, parse(Float64, m[2]))
		end

		m = match(r"Proportion of invariant[ ]?:\s+(.+)", line)
		if m != nothing
			rs = [0.0; rs]
			props = [parse(Float64, m[1]); props]
		end

		if occursin(r"^Model\tweight", line)
			is_mixture = true
			continue
		end
		
		if is_mixture
			#Model	weight		Pinvar		Rate
			#Buried	0.382513	0.107340	0.427673
			#Intermediate	0.436971	0.107340	0.837596
			#HExposed	0.180516	0.107340	1.518636
			if occursin(r"^Model", line)
				continue
			elseif occursin(r"^Confidence", line)
				is_mixture = false
				continue
			elseif match(r"[0-9]+", line) != nothing
				push!(Qrs, parse(Float64, line_arr[4]))
				push!(freqs, parse(Float64, line_arr[2]))
				continue
			else
				is_mixture = false
				continue
			end
		end
	end

	if ncat != 1 && (isempty(rs))
		rs, props = generate_rs_props_from_alpha(alpha, ncat)
 	elseif ncat != 1 && length(rs) == 1 #length(rs) == 1 in case of invariant site
		rs_gamma, props_gamma = generate_rs_props_from_alpha(alpha, ncat)
		rs = [rs; rs_gamma]
		props = [props; props_gamma]
	end

	return(Dict(:rs=>rs, :props=>props, :Fs=>Fs, :Qrs=>Qrs, :freqs=>freqs))
end


function generate_rs_props_from_alpha(α, ncat=4)
	quantiles = quantile.(Gamma(α, 1/α), range(0, 1, length=ncat+1))
	function mean_gamma_rate(a, b, shape, rate)
		dist = Gamma(shape, 1/rate)
		integral, _ = quadgk(x -> x * pdf(dist, x), a, b)
		prob = cdf(dist, b) - cdf(dist, a)
		integral / prob
	end
	# Compute rates for each bin
	rates = [mean_gamma_rate(quantiles[i], quantiles[i+1], α, α) for i in 1:ncat]
	# Normalize
	normalized_rates = rates ./ mean(rates)
	props = fill(1.0/ncat, ncat)
	return(normalized_rates, props)
end


function get_Qrs_freqs(Fs, Qrs, freqs, is_pmsf, pmsf_file, site2pattern, sub_model, mix_freq_model)
	q_pis_sites = Vector()
	q_pis = Vector()
	if(is_pmsf)
		pmsf_pis = read_pmsf_file(pmsf_file)
		pmsf_pis = convert_pmsf_file_by_site(pmsf_pis, site2pattern)
		# Qs q_pis_sites
		#q_pis_sites = generate_Qs(sub_model, pmsf_pis)
		generate_Qs!(q_pis_sites, sub_model, pmsf_pis)
	else
		is_mix_freq = false
		if mix_freq_model != nothing || mix_freq_model == ""
			Fs = isempty(Fs) ? Fs : [ first(Fs) ]
			q_pis, is_mix_freq = generate_Q(mix_freq_model) 
			[ push!(Fs, fs) for fs in map(x->x.Pi, q_pis) ]
			q_pis, _ = generate_Q(sub_model)
		else
			q_pis, _ = generate_Q(sub_model) 
		end
		println(std_outfh, join(["q_pis size:", size(q_pis)], "\t") * "\n")
		
		if (! isempty(Fs)) && (! isempty(Fs[1]))
			q_pis = generate_Q_w_F(q_pis, Fs, is_mix_freq)
		end
	end

	if isempty(Qrs)
		Qrs = isempty(q_pis_sites) ? fill(1.0, length(q_pis)) : fill(1.0, length(q_pis_sites[1]))
	end
	if isempty(freqs)
		freqs = fill(1.0/length(Qrs), length(Qrs))
	end

	return (Qrs, q_pis, q_pis_sites, freqs)
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


#function is_complex_matrix(A)
#           return eltype(A) <: Complex
#end

is_complex_matrix(A) = A isa AbstractMatrix{<:Complex}


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


