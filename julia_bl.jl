#! /bin/env julia


######################################################################
# v0.23.2


######################################################################
#DIR = Base.source_path()
const LIB = "lib-julia"
push!(LOAD_PATH, joinpath(@__DIR__, LIB))

include(joinpath(LIB, "readArg.jl"))
include(joinpath(LIB, "Q.jl"))
include(joinpath(LIB, "util.jl"))
#include(joinpath(LIB, "const.jl"))


######################################################################
#CONSTANTS
const GET_BRANCH_ORDER_IN_APE = joinpath(@__DIR__, "additional_scripts", "get_branch_order_in_ape.R")
#const step = cbrt(eps())


######################################################################
using Dates

using LinearAlgebra
using StatsFuns: logsumexp
#using Distributions
#using Calculus
#using FiniteDiff

#using Folds
#using ConcurrentCollections


######################################################################
struct Nb
    node::Int
	tip::Int
	branch::Int
end


######################################################################
function do_phylo_log_lk(q_pis::Vector{Any}, param::Vector{Float64}, liks_ori::Array{Float64,2}; sub_model::String, nb::Nb, nl::Int, all_children, rs::Vector{Float64}=[1.0], props::Vector{Float64}=[1.0], is_do_inv::Bool=false, inv_prop::Float64=0.0, Qrs::Vector{Float64}, freqs::Vector{Float64}, index::Int)
    bl = param

	# gamma rate
	#rs = [1.0]

	# rate's weight
	#props = fill(1/length(rs), length(rs))
	props = props ./ sum(props)

	# for inv
	if(is_do_inv)
		gap_cols = all(liks_ori[:, 1:nb.tip] .== 1.0, dims=2)
		res = liks_ori[:, 1:nb.tip][ .! gap_cols[:], : ]
		is_inv_site = all_rows_identical(res)
		if(is_inv_site)
			rs = [0; rs]
			props = [inv_prop; props * (1-inv_prop)]
		end
	end

	# Qs' freq (or weight)
	#freqs = [0.2, 0.8]

	# final weights
	weights = props .* freqs'
	weights = weights'[:] # matrix converted to vector
	if occursin("LG4M|LG4M", sub_model)
		weights = props
	end

	liks = similar(liks_ori)

	# other related params
	log_lk = 0.0
	lks = zeros(Float64, length(weights))
	counter = 0

	@inbounds for r_i in eachindex(rs), q_i in eachindex(q_pis)
		#LG4M or LG4X
		if occursin(r"^(LG4M|LG4X)$", sub_model) && r_i != q_i
			continue
		end
		counter += 1

		Q = q_pis[q_i].Q
		Pi = q_pis[q_i].Pi
		U = q_pis[q_i].U; Lambda = q_pis[q_i].Lambda; U_inv = q_pis[q_i].U_inv
		r = Qrs[q_i] * rs[r_i] # r defined

		#liks = copy(liks_ori)
		copyto!(liks, liks_ori)
		j = 0

		comp = zeros(Float64, nb.tip+nb.node)

		@simd for anc in ((nb.tip+nb.node):-1:(1+nb.tip))
			children = all_children[anc]
			m = zeros(Float64, nl, length(children))
			@simd for i in 1:length(children)
				j += 1
				if any(is_complex_matrix, [U, Lambda, U_inv])
					m[:,i] = exp(Q*bl[j]*r) * liks[:,children[i]]
				else
					eigendecom_pt!(view(m,:,i), U, Lambda, U_inv, liks[:,children[i]], bl[j]*r)
				end
			end
			m_prod = prod(m, dims=2)
			comp[anc] = ifelse(anc == (1 + nb.tip), dot(Pi, m_prod), sum(m_prod))
			liks[:, anc] = m_prod / comp[anc]
		end
		comp .= clamp.(comp, 1e-300, Inf)
		lks[counter] = sum(log.(comp[nb.tip+1:end]))
	end

	if(is_do_inv)
		if(is_inv_site)
			log_lk = logsumexp(log.(weights) .+ lks)
		else
			lks[isnan.(lks)] .= 0
			inv_weights = weights * (1-inv_prop)
			log_lk = logsumexp(log.(inv_weights) .+ lks)
		end
	else
		log_lk = logsumexp(log.(weights) .+ lks)
	end

	#return (ismissing(log_lk) ? Inf : log_lk)
	return([log_lk])
end


######################################################################
function get_liks_pattern(pattern::Vector)
	pattern2 = Vector()
	for x in pattern
		liks = zeros(Float64, nl, nb.tip + nb.node)
		for i in 1:nb.tip
			if x[1][i] == 999 # gap
				liks[:, i] = ones(nl)
			else
				liks[x[1][i], i] = 1.0
			end
		end
		push!(pattern2, (liks, x[2]))
	end
	return(pattern2)
end


@inline function eigendecom_pt!(out, U::Matrix{Float64}, Lambda::Vector{Float64}, U_inv::Matrix{Float64}, v::Vector{Float64}, scale::Float64=1.0)
	tmp = similar(U[:,1])
	#out = similar(tmp)
	# 1
	mul!(tmp, U_inv, v)
	# 2
    @inbounds @simd for k in eachindex(Lambda, tmp)
        tmp[k] *= exp(Lambda[k] * scale)
    end
	# 3
    mul!(out, U, tmp)
    return out
end


######################################################################
function delta_log_lk(bls::Vector{Float64}, indices::Vector{Int}, nums::Vector{Float64}, i::Int, pattern2::Vector, sub_model::String, nb::Nb, nl::Int, rs, props, inv_info, Qrs, freqs, all_children)
    site_pat, times_of_pattern = pattern2[i]
    q = is_pmsf ? q_pis_sites[i] : q_pis

	#bls_copy = copy(bls)
	bls_copy = similar(bls)
    copyto!(bls_copy, bls)

    bls_copy[indices] .+= nums

    # Pass all constant parameters as keyword arguments
    rv = do_phylo_log_lk(q, bls_copy, site_pat;
						 sub_model=sub_model, nb=nb, nl=nl, all_children=all_children,
                         rs=rs, props=props, 
                         is_do_inv=inv_info[:is_do_inv],
                         inv_prop=inv_info[:inv_prop],
                         Qrs=Qrs, freqs=freqs, index=i)[1] * times_of_pattern
    return rv
end


function hessian_STK2004(bls::Vector{Float64}, pattern2::Vector, nb::Nb, nl::Int, q_pis, q_pis_sites, inv_info::Dict, all_children::Dict, descendants::Dict)
	lnLs = zeros(Float64, length(pattern2))
	lnL0 = Float64(0)

	h = Dict(i => zeros(Float64, nb.branch, nb.branch) for i in 1:length(pattern2))
	gradients = [ zeros(Float64, nb.branch) for _ in 1:length(pattern2) ]

	# pre-assign argu
	#delta_log_lk2 = (bls, indices, nums, k) -> delta_log_lk(bls, indices, nums, k, pattern2)
	delta_log_lk2 = (bls, indices, nums, k) -> delta_log_lk(bls, indices, nums, k, pattern2, sub_model, nb, nl, rs, props, inv_info, Qrs, freqs, all_children)

	@inbounds Threads.@threads for k in 1:length(pattern2)
		lnLs[k] = delta_log_lk2(bls, [1], [0.0], k)
		thetas = zeros(Float64, nb.branch)
		@inbounds @simd for i in 1:nb.branch
			step = (1 + bls[i]) ./ 1e6
			# already considers site_pattern_freqs
			thetas[i] = (delta_log_lk2(bls,[i],[step],k) - delta_log_lk2(bls,[i],[-step],k))/(2*step)
		end
		gradients[k] += thetas

		# calculate h[k]
		#h[k] += thetas .* thetas' / pattern2[k][end]
		alpha = 1.0 / pattern2[k][end]
		LinearAlgebra.BLAS.ger!(alpha, thetas, thetas, h[k])
	end

	h = -sum(values(h)) # negative Hessian, this is a must regardless of whether the sign of lnL is correct
	gradients = sum(values(gradients))
	lnL0 = sum(lnLs)
	println()

	return (h, gradients, lnL0, lnLs)
end


function hessian_fd(bls::Vector{Float64}, sum_phylo_log_lk2::Function)
	h = ones(nb.branch, nb.branch)
	function calculate_lnL(bls0, indices, nums)
		a = deepcopy(bls0)
		a[indices] .+= nums
		loglk = sum_phylo_log_lk2(a)
	end
	
	f0 = calculate_lnL(bls,[1],(0))
	f = calculate_lnL

	@inbounds Threads.@threads for i in 1:nb.branch
		@inbounds for j in 1:nb.branch
			if h[i,j] != 1.0
				println([i,j])
				h[i,j] = h[j,i]
				continue
			end
			bls0 = copy(bls) #bls is mutable
			△bls = max.((bls0[i], bls0[j])./1e6, 1e-6)
			(bl1, bl2) = △bls
			if i == j
				v1 = ( f(bls,[i],(bl1)) - 2*f0 + f(bls,[i],(-bl1)) ) / (reduce(*,△bls))
			else
				v1 = (f(bls0,[i,j],(bl1,bl2)) + f(bls0,[i,j],(-bl1,-bl2)) - f(bls0,[i,j],(-bl1,bl2)) - f(bls0,[i,j],(bl1,-bl2))) / (4*reduce(*,△bls))
			end
			h[i,j] = v1
		end
	end

	println()
	#h .= -h

	return(h)
end


######################################################################
@inline function get_new_order(bl_order::Dict{String, Dict{Any, Any}})
	sorted_keys = sort(collect(keys(bl_order["mcmctree2ape"])))
	new_order = [bl_order["mcmctree2ape"][k] for k in sorted_keys]
end

@inline function get_new_order_rev(bl_order::Dict{String, Dict{Any, Any}})
	sorted_keys = sort(collect(keys(bl_order["ape2mcmctree"])))
	new_order = [bl_order["ape2mcmctree"][k] for k in sorted_keys]
end


function sum_phylo_log_lk(param::Vector{Float64}, pattern2::Vector{Any})
	# pattern2: v[1]=liks, v[2]=num
	sum( map(i -> do_phylo_log_lk((is_pmsf ? q_pis_sites[i] : q_pis), param, pattern2[i][1];
		sub_model = sub_model, nb=nb, nl=nl, all_children=all_children,
		rs = rs, props = props, is_do_inv=inv_info[:is_do_inv], inv_prop=inv_info[:inv_prop],
		Qrs = Qrs, freqs = freqs, index = i
	)[1]*pattern2[i][2], 1:length(pattern2)) )
end


function calculate_gradient(bls::Vector{Float64}, sum_phylo_log_lk2::Function)
	delta_loglks = Vector{Float64}(undef, length(bls))
	bls = map!(x->x<=1e-6 ? 0 : x, bls, bls)
	for i in 1:length(bls)
		zs = zeros(length(bls))
		#step = 2*(1e-3+bls[i]) / 1e6
		zs[i] = step
		delta_loglk = (sum_phylo_log_lk2(bls + zs) - sum_phylo_log_lk2(bls - zs)) / (2*step)
		delta_loglks[i] = delta_loglk
	end
	return(delta_loglks)
end


######################################################################
function compare_true_vs_approx_lnL(bls_vec::Vector{Float64}, gs, hs, lnL0, sum_phylo_log_lk2, transform_method)

	# important
	bls = bls_vec

	hs = [x for x in hs if x !== nothing]

	is_mvn = true
	unif_radius = 2

	covs = map(x -> -inv( Hermitian(x + Diagonal(ones(dim(x)[1]))*1e-6) )*2, hs)

	for i in 1:100
		new_bls = copy(bls)
		if is_mvn
			new_bls = max.(vec(rand(MvNormal(bls,covs[1]),1)), 1e-6)
		elseif unif_radius != nothing
			σ = sqrt.(diag(covs[1]))
			new_bls = [rand(max(Uniform(bls[i] - unif_radius*σ[i], 1e-6), bls[i] + unif_radius*σ[i])) for i in eachindex(bls)]
		end
		lnL = sum_phylo_log_lk2(new_bls)

		lnL_approxes = Float64[]
		for (g,h) in zip(gs, hs)
			if transform_method != nothing
				#println(transform_method, "!!!")
				bls2 = sqrt.(bls) #bls2 := u
				new_bls2, gu, hu = transform_bl(new_bls, g, h)
				delta_bls = new_bls2 - bls2
				lnL_approx = lnL0 + (delta_bls' * gu) + (delta_bls' * hu * delta_bls/2)
			else
				delta_bls = new_bls - bls
				lnL_approx = lnL0 + (delta_bls' * g) + (delta_bls' * h * delta_bls/2)
			end
			push!(lnL_approxes, lnL_approx)
		end
		println(join(vcat("TREE$i", lnL, lnL_approxes...), "\t"))
	end
end


function compare_true_vs_approx_lnL(bls_vec::Vector{Vector{Float64}}, gs, hs, lnL0, sum_phylo_log_lk2, transform_method)
	hs = [x for x in hs if x !== nothing]

	bls = popfirst!(bls_vec)
	#lnL0 = sum_phylo_log_lk2(bls)

	# here is for bs trees
	for bs_bls in bls_vec[1:min(100, length(bls_vec))]
		new_bls = copy(bs_bls)
		lnL = sum_phylo_log_lk2(new_bls)

		lnL_approxes = Float64[]
		for (g,h) in zip(gs,hs)
			if transform_method != nothing
				bls2 = sqrt.(bls) # bls := u
				new_bls2, gu, hu = transform_bl(new_bls, g, h)
				delta_bls = new_bls2 - bls2
				lnL_approx = lnL0 + (delta_bls' * gu) + (delta_bls' * hu * delta_bls/2)
			else
				delta_bls = new_bls - bls
				lnL_approx = lnL0 + (delta_bls' * g) + (delta_bls' * h * delta_bls/2)
			end
			push!(lnL_approxes, lnL_approx)
		end
		count = length(lnL_approxes)
		println(join(vcat("TREE$count", lnL, lnL_approxes...), "\t"))
	end
end


function transform_bl(bls, g, h)
	#Δ(u) = l(u) − l(u_hat) ≈ Δu'*gu + 1/2*Δu'*Hu*Δu,

	u = sqrt.(bls)
	db_over_du = 2 * u
	d2b_over_du2 = fill(2, length(u))
	#u = asin.(sqrt.((3/4) .- (3/4) * exp.(-4bls/3)))
	#db_over_du = cos(u/2)*sin(u/2) / (1 − 4/3*sin(u/2)^2)

	gu = g .* db_over_du
	hu = Diagonal(g.*d2b_over_du2) .+ h.*(db_over_du*db_over_du')

	return([u, gu, hu])
end


function read_hessian_infile(hessian_infile)
	in_fh = open(hessian_infile, "r")
	is_read_hessian = false
	count = 0
	hessian_data = Float64[]

	for line in readlines(in_fh)
		line = strip(line)
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


function get_h_from_inBV(hessian_infile, bl_order)
	h_from_inBV_mcmctree = hessian_infile == nothing ? nothing : read_hessian_infile(hessian_infile)
	h_from_inBV = nothing
	if h_from_inBV_mcmctree != nothing
		new_bl_order_rev = get_new_order_rev(bl_order)
		h_from_inBV = h_from_inBV_mcmctree[new_bl_order_rev, new_bl_order_rev]
	end
	return(h_from_inBV)
end


######################################################################
#function main(rs, props, Fs, Qrs, freqs, pattern)
function main(nb, nl, sub_model, rs, props, Fs, Qrs, freqs, q_pis, q_pis_sites, inv_info, pattern, all_children, descendants)

	pattern2 = get_liks_pattern(pattern)
	bls = fill(0.1, nb.tip+nb.node-1)
	bls_vec = Vector()

	# in case it contains bs trees
	bls = read_bls_from_branchout_matrix(branchout_matrix)[1]
	if bs_branchout_matrix != nothing
		bls_vec = read_bls_from_branchout_matrix(bs_branchout_matrix)
	end

	# reordered according MCMCtree's order
	bl_order = get_bl_order(branchout_matrix)
	new_bl_order = get_new_order(bl_order)

	println(now())

	# pre-assign pattern2
#=
	lnL0 = sum_phylo_log_lk2(bls)
	println(join(["best lnL", lnL0, "\n"], "\t"))
	println(std_outfh, join(["best lnL", lnL0, "\n"], "\t"))
	println(now())
=#

	sum_phylo_log_lk2 = (bls::Vector{Float64}) -> sum_phylo_log_lk(bls, pattern2)

	# Hessian
	if hessian_outfile != nothing
		h, g, lnL0, lnLs = hessian_STK2004(bls, pattern2, nb, nl, q_pis, q_pis_sites, inv_info, all_children, descendants)
		g_mcmctree = g[new_bl_order]
		if occursin(r"^finite_difference|fd$", hessian_type)
			#sum_phylo_log_lk2 = (bls::Vector{Float64}) -> sum_phylo_log_lk(bls, pattern2)
			h = hessian_fd(bls, sum_phylo_log_lk2)
			#h = -FiniteDiff.finite_difference_hessian(sum_phylo_log_lk2, bls)
		end

		out_fh = open(hessian_outfile, "w")
		println(out_fh, "")
		println(out_fh, string(nb.tip))
		println(out_fh, "")
		println(out_fh, read(treefile, String))
		println(out_fh, join(bls[new_bl_order], " "))
		println(out_fh, "")
		println(out_fh, join(g_mcmctree, " "))
		write(out_fh, "\n\n")
		println(out_fh, "Hessian")
		println(out_fh, "")
		close(out_fh)

		#sorted_keys = sort(collect(keys(bl_order["mcmctree2ape"])))
		#new_order = [bl_order["mcmctree2ape"][k] for k in sorted_keys]
		#h_mcmctree = h[new_order, new_order]
		h_mcmctree = h[new_bl_order, new_bl_order]
		h_mcmctree_out = format_number.(h_mcmctree)

		open(hessian_outfile, "a") do out_fh
			writedlm(out_fh, h_mcmctree_out, ' ')
		end
		
		if is_compare
			h_from_inBV = get_h_from_inBV(hessian_infile, bl_order)
			if isempty(bls_vec)
				compare_true_vs_approx_lnL(bls, [g], [h, h_from_inBV], lnL0, sum_phylo_log_lk2, transform_method)
			else
				gs = [g, zeros(Float64, length(g))]
				compare_true_vs_approx_lnL(vcat([bls], bls_vec), gs, [h, h_from_inBV], lnL0, sum_phylo_log_lk2, transform_method)
				#compare_true_vs_approx_lnL(vcat([bls], bls_vec), [g], [h], lnL0, sum_phylo_log_lk2, transform_method)
			end
		end
	end
	
	println(join(["best lnL", lnL0, "\n"], "\t"))
	println(std_outfh, join(["best lnL", lnL0, "\n"], "\t"))

	close(std_outfh)
	println(now())

	println()

end


######################################################################
function get_args()
	println()
	println(now())
	rs, props, Fs, Qrs, freqs, inv_info = get_iqtree_params(iqtree_file, phyml_file)

	println(now())
	nb, pattern, all_children, cherry_nodes, descendants, site2pattern = read_basics(basics_indir)
	nl = type == "AA" ? 20 : 4

	Qrs, q_pis, q_pis_sites, freqs = get_Qrs_freqs(Fs, Qrs, freqs, is_pmsf, pmsf_file, site2pattern, sub_model, mix_freq_model)


	return (nb, nl, sub_model, rs, props, Fs, Qrs, freqs, q_pis, q_pis_sites, inv_info, pattern, all_children, descendants)
end


######################################################################
######################################################################
if abspath(PROGRAM_FILE) == @__FILE__
	nb, nl, sub_model, rs, props, Fs, Qrs, freqs, q_pis, q_pis_sites, inv_info, pattern, all_children, descendants = get_args()
	main(nb, nl, sub_model, rs, props, Fs, Qrs, freqs, q_pis, q_pis_sites, inv_info, pattern, all_children, descendants)
end


