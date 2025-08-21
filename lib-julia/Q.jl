#! /usr/bin/env julia


################################################
using DelimitedFiles
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


function generate_Qs!(q_pis_sites, sub_model, pmsf_pis)
	#Q = [-1 1/3 1/3 1/3; 1/3 -1 1/3 1/3; 1/3 1/3 -1 1/3; 1/3 1/3 1/3 -1]
	#Qs = [Q for _ in 1:num_patterns]
	q_pis0, _ = generate_Q(sub_model) 
	for i in 1:size(pmsf_pis, 1)
		q_pis = generate_Q_w_F(q_pis0, [pmsf_pis[i,:]])
		#get_eigen!(q_pis[1])
		push!(q_pis_sites, q_pis)
	end
	#return(q_pis_sites)
end


function generate_Q(model)
	q_pis, is_mix_freq = read_paml_matrix(model)
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
	for i ∈ ["regular", "mfm"], suffix ∈ [".dat", ".nex"]
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
	eigen_decomp = eigen(Q)
	U = eigen_decomp.vectors
	Lambda = eigen_decomp.values
	U_inv = inv(U)

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


