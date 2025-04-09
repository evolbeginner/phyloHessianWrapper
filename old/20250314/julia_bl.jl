#! /bin/env julia

# © Sishuo Wang [2024-2025]. All rights reserved.
println("© Sishuo Wang [2024-2025]. All rights reserved.", "\n")


######################################################################
# v0.18.7
# transform_method --transform


######################################################################
#DIR = Base.source_path()
DIR = script_dir = @__DIR__
#using MKL


lib = "lib-julia"
include(joinpath(lib, "readArg.jl"))
include(joinpath(lib, "Q.jl"))
include(joinpath(lib, "util.jl"))
#include(joinpath(lib, "const.jl"))

GET_BRANCH_ORDER_IN_APE = joinpath(DIR, "additional_scripts", "get_branch_order_in_ape.R")


######################################################################
using ArgParse
using Dates

using LinearAlgebra
#using Calculus
using FiniteDiff
using StatsFuns: logsumexp

using Folds
using ConcurrentCollections

using Distributions


######################################################################
#BLAS.set_num_threads(1)


######################################################################
inv_info = Dict(:is_do_inv=>false, :inv_prop=>0.0)

iqtree_params = get_params_from_iqtree(iqtree_file)
rs = iqtree_params[:rs]
props = iqtree_params[:props]
Fs = iqtree_params[:Fs]
Qrs = iqtree_params[:Qrs]
freqs = iqtree_params[:freqs]

display(iqtree_params)

if isempty(rs)
	rs = [1.0]
	props = [1.0]
elseif any(x->x==0, rs)
#elseif rs[1] == 0.0
	#rs = filter(!=(0), rs) |> x -> vcat(x, fill(0, count(==(0), rs))) # move ele==0 to the end of rs
	inv_info = Dict(:is_do_inv=>true, :inv_prop=>props[1])
	rs = iqtree_params[:rs][2:end]
	props = iqtree_params[:props][2:end]
end


######################################################################
# obtain bls
bls = Vector()


######################################################################
function do_phylo_log_lk(q_pis::Vector{Any}, param::Vector{Float64}, liks_ori::Array{Float64,2}, cherry::Dict=Dict(), flouri=Dict(), is_pmsf::Bool=false, is_flouri::Bool=false; rs::Vector{Float64}=[1.0], props::Vector{Float64}=[1.0], is_do_inv::Bool=false, inv_prop::Float64=0.0, Qrs::Vector{Float64}, freqs::Vector{Float64}, index::Int)
    bl = param

	# gamma rate
	#rs = [1]

	# rate's weight
	#props = fill(1/length(rs), length(rs))
	props = props ./ sum(props)

	# for inv
	if(is_do_inv)
		gap_cols = all(liks_ori[1:nb.tip,:] .== 1.0, dims=2)
		res = liks_ori[1:nb.tip, :][ .! gap_cols[:], : ]
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
	if startswith(sub_model, "LG4"); weights = props; end #LG4M

	# other related params
	log_lk = 0
	lks = zeros(length(weights))
	counter = 0

	@inbounds for r_i in 1:length(rs)
		for q_i in 1:length(q_pis)
			#LG4M or LG4X
			if startswith(sub_model, "LG4")
				if r_i != q_i; continue; end
			end
			counter += 1

			Q = q_pis[q_i].Q
			Pi = q_pis[q_i].Pi
			U = q_pis[q_i].U; Lambda = q_pis[q_i].Lambda; U_inv = q_pis[q_i].U_inv
			r = Qrs[q_i] * rs[r_i] # r defined

			#liks = similar(liks_ori)
			#[ liks[i] = liks_ori[i] for i in 1:length(liks_ori) ]
			liks = copy(liks_ori)
			j = 0

			comp = zeros(Float64, nb.tip+nb.node)

			@inbounds @simd for anc in ((nb.node + nb.tip):-1:(1 + nb.tip))
				children = all_children[anc] # used to be "anc-nb.tip"
				m = zeros(Float64, nl, length(children))

				desc_lik_pattern = [ liks[descendant,:] for descendant in descendants[anc] ]
			 	@inbounds @simd for i in 1:length(children)
					j += 1
					if any(is_complex_matrix, [U, Lambda, U_inv])
						#U = convert_complex_to_float_in_matrix(U)
						#Lambda = convert_complex_to_float_in_matrix(Lambda)
						#U_inv = convert_complex_to_float_in_matrix(U)
						@inbounds m[:,i] = exp(Q*bl[j]*r)*liks[children[i], :]
					else
						@inbounds m[:,i] = eigendecom_pt(U,Lambda,U_inv,bl[j]*r)*liks[children[i], :]
					end
				end
				m_prod = prod(m, dims=2)
				@inbounds comp[anc] = (anc == (1 + nb.tip)) ? dot(Pi, m_prod) : sum(m_prod)
				@inbounds liks[anc, :] = m_prod / comp[anc]
			end
			# very rare numerical error
			if any(x->x<0, comp)
				println(stderr, "Error" * "\t" * string(index) * "\t")
				println(stderr, findall(x->x<0, comp))
			end
			comp[comp.<0] .= 1e-10
			lks[counter] = sum(log.(comp[nb.tip+1:end]))
		end
	end

	if(is_do_inv)
		if(is_inv_site)
			#log_lk = -log(sum(weights .* exp.(lks)))
			log_lk = -logsumexp(log.(weights) .+ lks)
		else
			lks[isnan.(lks)] .= 0
			inv_weights = weights * (1-inv_prop)
			log_lk = -logsumexp(log.(inv_weights) .+ lks)
		end
	else
		log_lk = -logsumexp(log.(weights) .+ lks)
	end

	#return (ismissing(log_lk) ? Inf : log_lk)
	return ([log_lk, flouri])
end


######################################################################
function get_liks_pattern(pattern::Vector)
	pattern2 = Vector()
	for x in pattern
		liks = zeros(Float64, nb.tip + nb.node, nl)
		for i in 1:nb.tip
			if x[1][i] == 999 # gap
				liks[i,:] = fill(1, nl)
			else
				liks[i, x[1][i]] = 1.0
			end
		end
		push!(pattern2, (liks, x[2]))
	end
	pattern2 = sort(pattern2, by=x->x[2], rev=true)
end


function calculate_lk_cherry(cherry_nodes::Vector{Int}, param::Vector{Float64})
	cherry_liks = Dict(i => Dict{Any,Any}() for i in cherry_nodes)
	#cherry_liks = Dict(i => Dict{Vector{Float64},Array{Float64}}() for i in cherry_nodes)
	for cherry_node in cherry_nodes
		cherry_liks[cherry_node] = Dict(Folds.map((x)->generate_cherry_liks(x[1],x[2],cherry_node,param), Iterators.product(1:nl, 1:nl)))
	end
	return(cherry_liks)
end


function generate_cherry_liks(i,j,cherry_node,param)
	v1 = zeros(Float64, 4)
	v2 = zeros(Float64, 4)
	v1[i] = 1.0; v2[j] = 1.0
	m = (nb.node+nb.tip - cherry_node + 1) * 2
	m_prod = (eigendecom_pt(param[m-1]) * v1) .* (eigendecom_pt(param[m]) * v2)
	return( (v1,v2)=>m_prod )
end


#@inline function eigendecom_pt(U::SMatrix{Float64}, Lambda::Vector{Float64}, U_inv::SMatrix{Float64}, scale::Float64=1.0)
@inline function eigendecom_pt(U::Matrix{Float64}, Lambda::Vector{Float64}, U_inv::Matrix{Float64}, scale::Float64=1.0)
	nrow = size(U)[1]
	P1 = similar(U)
	P2 = similar(U)
	#P1 = Matrix{Float64}(undef, nrow, nrow)
	#P2 = Matrix{Float64}(undef, nrow, nrow)
	mul!(P1, U, Diagonal(exp.(Lambda * scale)))
	#mul!(P2, U*Diagonal(exp.(Lambda * scale)), U_inv)
	mul!(P2, P1, U_inv)
	return(P2)
	#println(typeof(U)); println(typeof(Lambda))
	@fastmath U * Diagonal(exp.(Lambda * scale)) * U_inv
end


######################################################################
function delta_log_lk(bls::Vector{Float64}, indices::Vector{Int}, nums, i, pattern2, cherry, flouri)
	bls_copy = copy(bls)
	bls_copy[indices] .+= nums
	rv = do_phylo_log_lk((is_pmsf ? q_pis_sites[i] : q_pis), bls_copy, pattern2[i][1], cherry, flouri, is_pmsf, false;
		rs = rs, props = props, is_do_inv=inv_info[:is_do_inv], inv_prop=inv_info[:inv_prop],
		Qrs = Qrs, freqs = freqs, index = i)[1] * pattern2[i][2]
end


#function hessian_SKT2004(bls::Vector{Float64}, pattern2::Vector, bl_order::Dict{String, Dict{Any, Any}})
function hessian_SKT2004(bls::Vector{Float64}, pattern2::Vector)
	flouri=Dict()
	#is_pmsf=false
	cherry = Dict("cherry_nodes"=>[], "cherry_liks"=>Dict())
	#h = zeros(Float64, nb.branch, nb.branch)
	h = Dict(i => zeros(Float64, nb.branch, nb.branch) for i in 1:length(pattern2))

	# pre-assign argu
	delta_log_lk2 = (bls, indices, nums, k) -> delta_log_lk(bls, indices, nums, k, pattern2, cherry, flouri)

	Threads.@threads for k in 1:length(pattern2)
		#println(join(["site:\t", k], "\t"))
		thetas = zeros(Float64, nb.branch)
		for i in 1:nb.branch
			step = (1 + bls[i]) ./ 1e6
			thetas[i] = (delta_log_lk2(bls,[i],(step),k) - delta_log_lk2(bls,[i],(-step),k))/(2*step)
		end
		h[k] += thetas .* thetas'
	end
	h = -sum(values(h)) # negative Hessian
	println()

	return(h)
end


function hessian_fd(bls::Vector{Float64}, sum_phylo_log_lk2::Function)
	h = ones(nb.branch, nb.branch)
	function f(bls0, indices, nums)
		a = deepcopy(bls0)
		a[indices] .+= nums
		loglk = sum_phylo_log_lk2(a)
	end
	
	f0 = f(bls,[1],(0))

	Threads.@threads for i in 1:nb.branch
		for j in 1:nb.branch
			if h[i,j] != 1.0
				h[i,j] = h[j,i]
				continue
			end
			bls0 = deepcopy(bls) #bls is mutable
			△bls = max.((bls0[i], bls0[j])./1e6, 1e-6)
			#△bls = (1, 1) ./ 1e6
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
	h .= -h

	return(h)
end


######################################################################
@inline function get_new_order(bl_order::Dict{String, Dict{Any, Any}})
	sorted_keys = sort(collect(keys(bl_order["mcmctree2ape"])))
	new_order = [bl_order["mcmctree2ape"][k] for k in sorted_keys]
end


#function sum_phylo_log_lk(param, pattern2, q_pis, q_pis_sites)
function sum_phylo_log_lk(param::Vector{Float64}, pattern2::Vector{Any})
	# pattern2: v[1]=liks, v[2]=numsum_phylo_log_lk
	cherry = Dict("cherry_nodes"=>[], "cherry_liks"=>Dict())

	flouri = ConcurrentDict{Int, ConcurrentDict}()
	for int_node in keys(descendants)
		#flouri[int_node]=ConcurrentDict{Vector{Vector},Tuple}()
		flouri[int_node]=ConcurrentDict{Any,Any}()
	end
	#Folds.sum(Folds.map(i->do_phylo_log_lk(param,pattern2[i][1],cherry,flouri,Q[i],is_pmsf,true)[1]*pattern2[i][2], 1:length(pattern2)))
	if(true)
		Folds.sum( Folds.map(i -> do_phylo_log_lk((is_pmsf ? q_pis_sites[i] : q_pis), param, pattern2[i][1], cherry, flouri, is_pmsf, false; 
			rs = rs, props = props, is_do_inv=inv_info[:is_do_inv], inv_prop=inv_info[:inv_prop],
			Qrs = Qrs, freqs = freqs, index = i
		)[1]*pattern2[i][2], 1:length(pattern2)) )
	else
		sum( map(i -> do_phylo_log_lk((is_pmsf ? q_pis_sites[i] : q_pis), param, pattern2[i][1], cherry, flouri, is_pmsf, false;
			rs = rs, props = props, is_do_inv=inv_info[:is_do_inv], inv_prop=inv_info[:inv_prop],
			Qrs = Qrs, freqs = freqs, index = i
		)[1]*pattern2[i][2], 1:length(pattern2)) )
	end
end


function calculate_gradient(bls::Vector{Float64}, sum_phylo_log_lk2::Function)
	delta_loglks = Vector{Float64}(undef, length(bls))
	bls = map!(x->x<=1e-6 ? 0 : x, bls, bls)
	for i in 1:length(bls)
		zs = zeros(length(bls))
		step = 2*(1e-3+bls[i]) / 1e6
		zs[i] = step
		delta_loglk = (sum_phylo_log_lk2(bls + zs) - sum_phylo_log_lk2(bls - zs)) / (2*step)
		#delta_loglk = (sum_phylo_log_lk2(bls + zs) - sum_phylo_log_lk2(bls)) / step
		delta_loglks[i] = delta_loglk
	end
	return(delta_loglks)
end


######################################################################
function compare_true_vs_approx_lnL(bls, g, h, lnL0, sum_phylo_log_lk2, transform_method)
	cov = -inv( Hermitian(h+Diagonal(ones(dim(h)[1]))*1e-6) )
	println(diag(cov))

	for i in 1:100
		#new_bls = bls .+ sqrt.(diag(cov)) .* max.(randn(length(bls)), 1e-6)
		#new_bls = bls .+ sqrt.(diag(cov)) .* max.(randn(length(bls)), 1e-8)
		new_bls = max.(vec(rand(MvNormal(bls,cov),1)), 1e-6)
		lnL = -sum_phylo_log_lk2(new_bls)

		if transform_method != nothing
			bls2 = sqrt.(bls)
			new_bls2, gu, hu = transform_bl(new_bls, g, h)
			delta_bls = new_bls2 - bls2
			lnL_approx = lnL0 + (delta_bls' * gu) + (delta_bls' * hu * delta_bls/2)
		else
			delta_bls = new_bls-bls
			lnL_approx = lnL0 + (delta_bls' * g) + (delta_bls' * h * delta_bls/2)
		end

		println(join([lnL, lnL_approx, lnL-lnL_approx], "\t"))
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


######################################################################
function main(rs, props, Fs, Qrs, freqs, pattern)
#function main()
	pattern2 = get_liks_pattern(pattern)
	bls = fill(0.1, nb.tip+nb.node-1)

	bls = read_bls_from_branchout_matrix(branchout_matrix)
	bl_order = get_bl_order(branchout_matrix)
	# reordered according MCMCtree's order
	new_bl_order = get_new_order(bl_order)
	println(now())

	# pre-assign pattern2
	sum_phylo_log_lk2 = (bls::Vector{Float64}) -> sum_phylo_log_lk(bls, pattern2)
	lnL0 = -sum_phylo_log_lk2(bls)
	println(join(["best lnL", lnL0, "\n"], "\t"))
	println(std_outfh, join(["best lnL", lnL0, "\n"], "\t"))

	println(now())

	# Gradient
	println("Gradients")
	g = zeros(length(bls))
	#g = FiniteDiff.finite_difference_gradient(sum_phylo_log_lk2, bls)
	#g = round.(g, digits=6)
	#g_mcmctree = g[new_bl_order]
	#println(.-g_mcmctree)

	g = calculate_gradient(bls, sum_phylo_log_lk2)
	g = round.(g, digits=6)
	g = .-g
	g_mcmctree = g[new_bl_order]
	println(g_mcmctree)
	println()
	println(now())

	# Hessian
	if hessian_outfile != nothing
		out_fh = open(hessian_outfile, "w")
		#println(out_fh, lnL0)
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

		if hessian_type == "SKT2004"
			h = hessian_SKT2004(bls, pattern2)
		elseif occursin(r"^finite_difference|fd$", hessian_type)
			h = hessian_fd(bls, sum_phylo_log_lk2)
			#h = -FiniteDiff.finite_difference_hessian(sum_phylo_log_lk2, bls)
		end

		sorted_keys = sort(collect(keys(bl_order["mcmctree2ape"])))
		new_order = [bl_order["mcmctree2ape"][k] for k in sorted_keys]
		h_mcmctree = h[new_order, new_order]
		h_mcmctree_out = format_number.(h_mcmctree)

		open(hessian_outfile, "a") do out_fh
			writedlm(out_fh, h_mcmctree_out, ' ')
		end
		println()
		
		#compare_true_vs_approx_lnL(bls, g, h, lnL0, sum_phylo_log_lk2, transform_method)
	end

	close(std_outfh)
	println(now())

	exit()


	println()

function hessian_fd2(bls)
	function f(bls0, indices, nums)
		a = deepcopy(bls0); a[indices] .+= nums; loglk = sum_phylo_log_lk2(a)
		#bls0[indices] .+= nums; loglk = sum_phylo_log_lk2(bls0)
	end

	for i in 1:nb.branch, j in 1:nb.branch
		bls0 = deepcopy(bls) #bls is mutable
		#bls[1]=ori_bls[1]+2h; println(bls[1]); f1 = sum_phylo_log_lk(bls)
		#bls[1]=ori_bls[1]-2h; println(bls[1]); f2 = sum_phylo_log_lk(bls)
		#v1 = (f(bls,i,2h)+f(bls,i,-2h)-2*lnL0)/(4h^2)
		△bls = (bls0[i], bls0[j]) ./ 1e6
		(bl1, bl2) = △bls
		v1 = (f(bls0,[i,j],(bl1,bl2)) + f(bls0,[i,j],(-bl1,-bl2)) - f(bls0,[i,j],(-bl1,bl2)) - f(bls0,[i,j],(bl1,-bl2)))/(4*bl1*bl2)
		#println(f(bls0,[i,j],(bl1,bl2)))
		#v1 = (f(bls0,[i,j],(bl1,bl2)) + f(bls0,[i,j],(-2bl1,-2bl2)) - f(bls0,[i,j],(0,2bl2)) - f(bls0,[i,j],(2bl1,-2bl2)))/(4*bl1*bl2)

		#v2 = (-f(bls,i,4h) + 16*f(bls,i,2h) + lnL0 - 16*lnL0 + lnL0 - 16*lnL0 - f(bls,i,-4h) + 16*f(bls,i,-2h))/48h^2
		#v2 = (-f(bls,(i,j),(2bl1,2bl2)) + 16f(bls,(i,j),(bl1,bl2)) + f(bls,(i,j),(2bl1,-2bl2)) - 16f(bls,(i,j),(bl1,-bl2)) + f(bls,(i,j),(-2bl1,2bl2)) - 16f(bls,(i,j),(-bl1,bl2)) - f(bls,(i,j),(-2bl1,-2bl2)) + 16f(bls,(i,j),(-bl1,-bl2)))/(48*bl1*bl2)
		print(i,"-",j,":",round.((v1),digits=3),"\t")
	end
end

	println()
	#println(i, "\t", bls[i], "\t", "order(n^2, n^4)", "\t", v1, "\t", v2)
	#g = FiniteDiff.finite_difference_gradient(sum_phylo_log_lk, bls); println(g)
	#h = hessian((x) -> sum_phylo_log_lk(x), bls); display(h); exit()
	#h = FiniteDiff.finite_difference_hessian(sum_phylo_log_lk, bls); display(h); exit()
	println(now())

	exit()
end


######################################################################
# Get "all_children"
#infile = "julia/all_children"
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
#infile = "julia/pattern"
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
#infile = "julia/basics"
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


struct Nb
    node::Int
	tip::Int
	branch::Int
end
nb = Nb(basics["nb.node"], basics["nb.tip"], basics["nb.node"]+basics["nb.tip"]-1)


# Get cherry nodes
#infile = "julia/cherry"
infile = joinpath(basics_indir, "cherry")
cherry_nodes = Vector{Int}()
for line in eachline(infile)
    push!(cherry_nodes, parse(Int, line))
end


# Get descendants
#infile = "julia/descendants"
infile = joinpath(basics_indir, "descendants")
#descendants = Dict{Int,Dict}()
descendants = Dict()
for line in eachline(infile)
    line_arr = split(line, '\t')
	int_node = parse(Int, line_arr[1])
    descendants[int_node] = [parse(Int, x) for x in line_arr[2:end]]
end


######################################################################
# Q
q_pis_sites = Vector()

if(is_pmsf)
	pmsf_pis = read_pmsf_file(pmsf_file)
	# Qs q_pis_sites
	q_pis_sites = generate_Qs(sub_model, pmsf_pis)
	#display(q_pis_sites[3][1].Q); exit()
else
	is_mix_freq = false

	# solve mixture like C2, but not work when EX2+F

	if mix_freq_model != nothing || mix_freq_model == ""
		Fs = isempty(Fs) ? Fs : [ first(Fs) ]
		q_pis, is_mix_freq = generate_Q(mix_freq_model) 
		[ push!(Fs, fs) for fs in map(x->x.Pi, q_pis) ]
		q_pis, _ = generate_Q(sub_model)
	else
		q_pis, _ = generate_Q(sub_model) 
	end
	println(std_outfh, join(["q_pis size:", size(q_pis)], "\t") * "\n")
	#Fs[1] = fill(0.05, 20)
	#q_pis = generate_Q_w_F(q_pis, Fs)
	
	if (! isempty(Fs)) && (! isempty(Fs[1]))
		q_pis = generate_Q_w_F(q_pis, Fs, is_mix_freq)
	end
end

nl = type == "AA" ? 20 : 4



######################################################################
if isempty(Qrs)
	Qrs = isempty(q_pis_sites) ? fill(1.0, length(q_pis)) : fill(1.0, length(q_pis_sites[1]))
end
if isempty(freqs)
    freqs = fill(1/length(Qrs), length(Qrs))
end

#println(length(q_pis)); exit()


######################################################################
main(rs, props, Fs, Qrs, freqs, pattern)
#main()


