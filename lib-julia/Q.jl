#! /usr/bin/env julia


################################################
using DelimitedFiles
using LinearAlgebra
#using DataFrames
#using CSV
#using StaticArrays


################################################
convert_to_zero(x) = isa(x, AbstractString) ? 0.0 : x


################################################
mutable struct Q_Pi
	Q
	Pi
	R
	U
	Lambda
	U_inv
end


################################################
function read_pmsf_file(pmsf_file)
	data = readdlm(pmsf_file)
	pmsf_pis = convert(Array{Float64,2}, data)
	pmsf_pis = pmsf_pis[:, 2:end] # del the 1st col
	return(pmsf_pis)
end


function convert_pmsf_file_by_site(pmsf_pis, site2pattern)
	pmsf_pis2 = similar(pmsf_pis)
	for i in 1:size(pmsf_pis, 1)
		pmsf_pis2[site2pattern[i],:] = pmsf_pis[i,:]
	end
	pmsf_pis2 = pmsf_pis2[1:length(unique(values(site2pattern))),:]
	return(pmsf_pis2)
end


function generate_Qs!(q_pis_sites, sub_model, pmsf_pis; seq_type=nothing, iqtree_file=nothing)
	#Q = [-1 1/3 1/3 1/3; 1/3 -1 1/3 1/3; 1/3 1/3 -1 1/3; 1/3 1/3 1/3 -1]
	#Qs = [Q for _ in 1:num_patterns]
	q_pis0, _ = generate_Q(sub_model, seq_type=seq_type, iqtree_file=iqtree_file) 
	for i in 1:size(pmsf_pis, 1)
		q_pis = generate_Q_w_F(q_pis0, [pmsf_pis[i,:]])
		#get_eigen!(q_pis[1])
		push!(q_pis_sites, q_pis)
	end
	#return(q_pis_sites)
end


const DNA_STATES = ["A", "C", "G", "T"]
const DNA_STATE_INDEX = Dict(state => i for (i, state) in enumerate(DNA_STATES))
const DNA_SUB_MODELS = Set(["JC", "JC69", "F81", "K2P", "K80", "HKY", "HKY85", "TN", "TN93", "TNEF", "TPM2", "TPM3", "TIM", "TIM2", "TIM3", "TVM", "SYM", "GTR"])


function generate_Q(model; seq_type=nothing, iqtree_file=nothing)
	model_name = uppercase(split(model, '+')[1])
	if seq_type == "DNA" || model_name in DNA_SUB_MODELS
		if iqtree_file != nothing && isfile(iqtree_file)
			return generate_DNA_Q_from_iqtree(iqtree_file)
		elseif model_name in ["JC", "JC69"]
			return generate_DNA_Q(fill(1.0, 4, 4), fill(0.25, 4)), false
		else
			error("DNA substitution model $model needs an IQ-TREE .iqtree file with estimated DNA rates/frequencies.")
		end
	end
	q_pis, is_mix_freq = read_paml_matrix(model)
end


function generate_DNA_Q_from_iqtree(iqtree_file)
	R, Pi, Q = read_DNA_params_from_iqtree(iqtree_file)
	q_pi = Q === nothing ? generate_DNA_Q(R, Pi) : generate_DNA_Q_from_matrix(Q, Pi)
	get_eigen!(q_pi)
	return([q_pi], false)
end


function generate_DNA_Q(R, Pi)
	Pi = Pi ./ sum(Pi)
	R = copy(R)
	for i in 1:length(Pi)
		R[i, i] = 0.0
	end
	Q = [R[i, j] * Pi[j] for i in 1:length(Pi), j in 1:length(Pi)]
	q_pi = Q_Pi(Q, Pi, R, nothing, nothing, nothing)
	fill_diag_for_Q!(q_pi)
	return(q_pi)
end


function generate_DNA_Q_from_matrix(Q, Pi)
	Pi = Pi ./ sum(Pi)
	Q = copy(Q)
	q_pi = Q_Pi(Q, Pi, nothing, nothing, nothing, nothing)
	q_pi.Q .= q_pi.Q .- diagm(diag(q_pi.Q))
	fill_diag_for_Q!(q_pi)
	return(q_pi)
end


function read_DNA_params_from_iqtree(iqtree_file)
	R = ones(Float64, 4, 4)
	Pi = fill(NaN, 4)
	Q = Matrix{Float64}(undef, 0, 0)

	is_rate = false
	is_freq = false
	is_q = false
	q_rows = Vector{Vector{Float64}}()

	for line in eachline(iqtree_file)
		if occursin(r"^Rate parameter R:", line)
			is_rate = true
			is_freq = false
			is_q = false
			continue
		elseif occursin(r"^State frequencies", line)
			is_rate = false
			is_freq = true
			is_q = false
			continue
		elseif occursin(r"^Rate matrix Q:", line)
			is_rate = false
			is_freq = false
			is_q = true
			empty!(q_rows)
			continue
		end

		if is_rate
			m = match(r"^\s*([ACGT])-([ACGT]):\s*([-+0-9.eE]+)", line)
			if m != nothing
				i = DNA_STATE_INDEX[m[1]]
				j = DNA_STATE_INDEX[m[2]]
				R[i, j] = R[j, i] = parse(Float64, m[3])
			elseif !isempty(strip(line))
				is_rate = false
			end
		elseif is_freq
			m = match(r"^\s*pi\(([ACGT])\)\s*=\s*([-+0-9.eE]+)", line)
			if m != nothing
				Pi[DNA_STATE_INDEX[m[1]]] = parse(Float64, m[2])
			elseif !isempty(strip(line))
				is_freq = false
			end
		elseif is_q
			m = match(r"^\s*([ACGT])\s+(.+)$", line)
			if m != nothing
				vals = parse.(Float64, split(strip(m[2])))
				if length(vals) == 4
					push!(q_rows, vals)
					if length(q_rows) == 4
						Q = reduce(vcat, reshape.(q_rows, 1, :))
						is_q = false
					end
				end
			elseif !isempty(strip(line))
				is_q = false
			end
		end
	end

	if any(isnan, Pi)
		Pi .= 0.25
	end

	return(R, Pi, isempty(Q) ? nothing : Q)
end


function generate_Q_from_m(m, ncol)
	Pi = m[end, :]
	# deal w/ freq Pi 0.1 ... 0.2;
	if typeof(Pi[end]) in (String, SubString{String}) && occursin(r";$", Pi[end])
		Pi[end] = parse(Float64, replace(Pi[end], r";$" => ""))
	end

	R = m[1:end-1, :]
	R = vcat(reshape(zeros(ncol), 1, ncol), R)
	R = convert_to_zero.(R)
	copy_lower_tri_as_upper_tri!(R)

	Q = [ R[i,j]*Pi[j] for i in 1:ncol, j in 1:ncol ]

	#Pi = fill(1/ncol,ncol)
	#Q = fill(1/ncol, ncol, ncol)
	Q[1:size(Q,1)+1:end] .= 0

	neg_row_sums = -sum(Q, dims=2)
	for i in 1:ncol
		Q[i, i] = neg_row_sums[i, 1]
	end
	Q .= -Q / sum(Pi.*diag(Q))
	Pi = Pi ./ sum(Pi)

	return(Q_Pi(Q, Pi, R, nothing, nothing, nothing))
end


function read_paml_matrix(model; model_dir=joinpath(@__DIR__, "../substitution_model"))
	q_pis = Vector()

	paml_matrix_text = nothing
	for i ∈ ["regular", "mfm", "mfm/UDM"], suffix ∈ [".dat", ".nex"]
		model_dat_file = joinpath(model_dir, i, model * suffix)
		if isfile(model_dat_file)
			paml_matrix_text = model_dat_file
			break
		end
	end
	if ! isfile(paml_matrix_text)
		error("subs model $model does not exist in $model_dir ......")
	end

	m0 = readdlm(paml_matrix_text, ' ')
	m0 = m0[[!all(row .== "") for row in eachrow(m0)], :]
	m0 = m0[ [ isa(row[1], Number) || row[1] == "frequency" for row in eachrow(m0)], :]

	is_mix_freq = false

	if any(==("frequency"), m0)
		is_mix_freq = true
		m0 = map(row -> (row[1] == "frequency") ? row[4:end] : row, eachrow(m0))
		for i in 1:length(m0)
        	m0[i][end] = parse(Float64, replace(m0[i][end], ";" => ""))
		end
		m0 = hcat(m0...)'
	end

	nrow, ncol = size(m0)

	if ! is_mix_freq
		first, last = 1, ncol
		while(first < nrow)
			m = m0[first:last, :]
			first += ncol
			last = first + ncol - 1
			q_pi = generate_Q_from_m(m, ncol)
			q_pi.Q = -q_pi.Q / sum(q_pi.Pi .* diag(q_pi.Q))
			get_eigen!(q_pi)
			push!(q_pis, q_pi)
		end
	else
		for i in 1:size(m0, 1)
			q_pi = Q_Pi(nothing, m0[i,:], nothing, nothing, nothing, nothing)
			push!(q_pis, q_pi)
		end
	end

	return(q_pis, is_mix_freq)
end


function convert_to_numeric(ori_mat)
	numeric_m = similar(ori_mat)
	numeric_m = [
    	tryparse(Float64, "1") for x in ori_mat
	]
end


function check_numeric_string(str)
	return tryparse(Float64, str) !== nothing
end


function copy_lower_tri_as_upper_tri!(original_matrix)
	#new_matrix = similar(original_matrix)
	n = size(original_matrix,1)
	for i in 1:n
		for j in i+1:n
			original_matrix[i, j] = original_matrix[j, i]
		end
	end
end


function get_eigen!(q_pi)
	Q = q_pi.Q
	eigen_decomp = LinearAlgebra.eigen(Q)
	U = eigen_decomp.vectors
	Lambda = eigen_decomp.values
	U_inv = LinearAlgebra.inv(U)

	#U = SMatrix{size(U,1), size(U,2)}(U)
	#U_inv = SMatrix{size(U,1), size(U,2)}(U_inv)
	#Lambda = SVector{size(U,2)}(Lambda)

	q_pi.U = U
	q_pi.Lambda = Lambda
	q_pi.U_inv = U_inv
end


function generate_Q_w_F(q_pis, Pis, is_mix_freq=false)
	#new_q_pis = Array{Any}(undef, length(q_pis) * length(Pis))
	new_q_pis = Vector()

	if is_mix_freq #EX2+C2 => 4 matrices
		for (i, d) in enumerate(q_pis)
			for (j, Pi) in enumerate(Pis)
				push_into_q_pis!(new_q_pis, d, Pis[j])
			end
		end
	else #EX2+C2 => 2 matrices, or EX2+F
		for (i, d) in enumerate(q_pis)
			push_into_q_pis!(new_q_pis, d, Pis[i])
		end
	end
	return(new_q_pis)
end


function push_into_q_pis!(new_q_pis, d, Pi)
	q_pi = deepcopy(d)
	#q_pi = d
	Q = q_pi.Q

	# old_Pi
	old_Pi = q_pi.Pi
	ncol = length(Pi)

	q_pi.Q = [ Q[x,y] / old_Pi[y] for x in 1:ncol, y in 1:ncol ]
	#q_pi.Q ./= old_Pi'

	q_pi.Pi = Pi ./ sum(Pi) # normalize
	q_pi.Q = [ q_pi.Q[x,y]*Pi[y] for x in 1:ncol, y in 1:ncol ] # here q_pi.Q is actually R
	#q_pi.Q .*= Pi'
	fill_diag_for_Q!(q_pi)
	get_eigen!(q_pi)
	push!(new_q_pis, q_pi)
end


function fill_diag_for_Q!(q_pi)
	Q = q_pi.Q
	Pi = q_pi.Pi
	Q .= Q .- diagm(diag(Q))
	neg_row_sums = -sum(Q, dims=2)
	for i in 1:length(q_pi.Pi)
		Q[i, i] = neg_row_sums[i, 1]
	end
	Q .= -Q / sum(Pi.*diag(Q))
	q_pi.Q = Q
end


################################################
#generate_Q("EX2")
