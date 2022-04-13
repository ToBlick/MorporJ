### A Pluto.jl notebook ###
# v0.19.0

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 5128e10e-b5a2-11ec-2371-d3c4d4e24b53
begin
	import Pkg
	Pkg.add("HDF5")
	Pkg.add(url="https://github.com/ToBlick/MorporJ")
	using PlutoUI, HDF5, Plots, Statistics, LinearAlgebra, MorporJ
end;

# ╔═╡ be1f1362-34af-4902-88e4-032d3cfd22b9
md""" After reading in the snapshots, the gBar algorithm finds 5 atoms. 
"""

# ╔═╡ 789541fe-e001-44ab-bb70-e50036a209d4
md"""
	Select barycentric weights:\
	λ₁ = $(@bind λ₁ Slider(0:0.01:1; default=0.33, show_value=true))\
	λ₂ = $(@bind λ₂ Slider(0:0.01:1; default=0.33, show_value=true))\
	λ₃ = $(@bind λ₃ Slider(0:0.01:1; default=0.34, show_value=true))\
	λ₄ = $(@bind λ₄ Slider(0:0.01:1; default=0.0, show_value=true))\
	λ₅ = $(@bind λ₅ Slider(0:0.01:1; default=0.0, show_value=true))\
	"""

# ╔═╡ 25b561c3-a03f-4aaa-a0b8-5dc45657113d
begin
	md""" User input projected to Σₙ: λ =
	"""
end

# ╔═╡ 4c62cd98-aad7-4215-a879-21606463ab1a
pal = palette(:viridis, 10)[1:7];

# ╔═╡ 99dda5be-3c26-4119-b832-e42f2559ca12
begin
	### import snapshots and get icdfs
	# parameters
	N = 1002
	M = 2*N
	
	tol_icdf = 1e-12
	tol_cdf = 1e-12

	β = [2,3,4,5,6]
	μ = [1,2,3,6,12,25]
	t = collect(3:2:51)
	paramnames = ("β", "μ", "t")
	
	params = (1.0*β, 1.0*μ, 1.0*t)
	A = zeros(length(β), length(μ), length(t), N)
	ci = CartesianIndices(size(A)[1:(end-1)])
	
	# get snaps
	fid = h5open("../snapshots/data_homogeneous_mu_exp.h5", "r");

	for i in eachindex(β)
	    for j in eachindex(μ)
	            g = read(fid, "snap_mu_$(μ[j])")
	            dset = g["snap_exp_$(β[i])"]
	            for l in eachindex(t)
	                A[i,j,l,1] = 0    # add zero at position -2ε
	                A[i,j,l,2] = 1    # add 1 at position -ε
	                A[i,j,l,3:end] .= dset[t[l],:]
	            end
	    end
	end
	# Snapshot matrix: N × nₜₚ
	S = reshape(A,:,N)';
	nₜₚ = size(S)[2]
	close(fid);

	xgrid = collect(range(0, 1, length = N));
	pgrid = collect(range(0, 1, length = M));
	Δx = xgrid[2] - xgrid[1]; Δp = pgrid[2] - pgrid[1];
	
	# normalize
	S₁, mass = MorporJ.get_S₁(S, Δx)
	# icdfs
	s = MorporJ.get_s(S₁, Δx, xgrid, pgrid)
end;

# ╔═╡ 2a1f31cd-2879-4396-a062-a12a8d27fdb6
begin
	a, basis_index = MorporJ.get_initial_atoms(s)
	Λ, ΔW, ΔW_max, ΔW_avg, a, basis_index = MorporJ.greedy_algo!(s, a, Δp, basis_index; rtol=-1e3, tol=1e-6, max_iter_greedy=4, max_iter_solver=Int(1e7), eps_abs_solver=1e-5, eps_rel_solver=5e-4)
end;

# ╔═╡ f5b01747-23c0-422f-8681-ac4be31d6b3c
begin
	print("Parameter values chosen by the greedy algorithm: \n")
	print(paramnames)
	print("\n")
	for a in basis_index
	    (i,j,l) = Tuple(ci[a])
	    print(params[1][i], " ", params[2][j], " ", params[3][l])
	    print("\n")
	end
end

# ╔═╡ f4d5a587-9c24-4125-87dc-c54b3666b372
md"""Select parameters for snapshot to be approximated:\
β = $(@bind β_target Slider(β; default=2, show_value=true))\
μ = $(@bind μ_target Slider(μ; default=1, show_value=true))\
t = $(@bind t_target Slider(t; default=31, show_value=true))\
Select number of atoms (<=5):\
n = $(@bind n Slider(2:5; default=3, show_value=true))
"""

# ╔═╡ d3db0ff2-f40c-460e-9076-c2537c581fcf
begin
	λ = MorporJ.proj_to_simplex([λ₁, λ₂, λ₃, λ₄, λ₅][1:n])
end

# ╔═╡ e24b0c02-6ffd-4e79-8754-4616ac850e84
begin
	icdf_target = zeros(M)
	for i in 1:n
		icdf_target .+= λ[i] .* a[i]
	end
	cdf_target = zeros(N)
	s_target = zeros(N)
	MorporJ.inv_F!(cdf_target, icdf_target, pgrid, xgrid, tol_icdf)
    MorporJ.cdf_to_pdf!(s_target, cdf_target, Δx, 1)
    s_target ./= s_target[2]
end;

# ╔═╡ 82fdb311-369c-49b5-9c77-21b2ee4a29a6
begin
	i_target = 1
	# get the corresponding index
	for k in 1:nₜₚ
		(i,j,l) = Tuple(ci[k])
		if β_target == β[i] && μ_target == μ[j] && t_target == t[l]
			i_target = k
			break
		end
	end
end


# ╔═╡ fc395d85-0831-4c79-a392-9169d2ff090f
begin
	plot(xgrid, S[:,basis_index[1:n]], 
		xlabel="p", ylabel="x", 
		linewidth=[2 2 2 2 2], ylim=(0,1), 
	    label=[ "s(a₁)" "s(a₂)" "s(a₃)" "s(a₄)" "s(a₅)"],
	    legend = :topright, 
		palette = pal[1:n], 
		linealpha=[0.5 0.5 0.5 0.5 0.5])
	plot!(xgrid, S[:,i_target],
		linewidth=2, 
		linealpha=1, 
		label = "s(β,μ,t)",
		color = pal[end])
	plot!(xgrid, s_target,
		linewidth=2, 
		linealpha=1, 
		label = "B({λᵢ}ᵢ)",
		color = :red)
end

# ╔═╡ 2ec48ea2-f5ca-4f35-b83a-d7bccc0ca89a
ΔL1 = Δx * sum(abs.(s_target - S[:,i_target]));

# ╔═╡ 46ccf4bb-57ad-499e-a40b-c6e50bfb9e76
begin
	md""" L1 error: $(ΔL1)
	"""
end

# ╔═╡ Cell order:
# ╟─5128e10e-b5a2-11ec-2371-d3c4d4e24b53
# ╟─be1f1362-34af-4902-88e4-032d3cfd22b9
# ╟─2a1f31cd-2879-4396-a062-a12a8d27fdb6
# ╟─f5b01747-23c0-422f-8681-ac4be31d6b3c
# ╟─f4d5a587-9c24-4125-87dc-c54b3666b372
# ╟─789541fe-e001-44ab-bb70-e50036a209d4
# ╟─25b561c3-a03f-4aaa-a0b8-5dc45657113d
# ╟─d3db0ff2-f40c-460e-9076-c2537c581fcf
# ╟─fc395d85-0831-4c79-a392-9169d2ff090f
# ╟─46ccf4bb-57ad-499e-a40b-c6e50bfb9e76
# ╟─2ec48ea2-f5ca-4f35-b83a-d7bccc0ca89a
# ╟─e24b0c02-6ffd-4e79-8754-4616ac850e84
# ╟─4c62cd98-aad7-4215-a879-21606463ab1a
# ╟─82fdb311-369c-49b5-9c77-21b2ee4a29a6
# ╟─99dda5be-3c26-4119-b832-e42f2559ca12
