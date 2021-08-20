### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 37bf7c73-acda-42ee-acce-e10131f83e12
import Pkg; Pkg.add(url="https://github.com/ToBlick/MorporJ")

# ╔═╡ 1e8a18a4-c2ba-482b-9624-fdafc49521c9
using LinearAlgebra, Plots, MorporJ, PlutoUI

# ╔═╡ c3dd7a97-aea6-4f58-80a4-175803b55167
begin
	n=2
	x = collect(0:0.001:20)
	N = length(x)
	Δx = x[2]-x[1]
	#f₁ = 0.5*exp.(-(x.-3.5).^2 ./ 1.2) + 0.5*exp.(-(x.-5.5).^2 ./0.7)
	f₁ = zeros(N)
	for i in eachindex(x)
		x[i] < 10 ? f₁[i] = maximum([sin(1.2*(x[i])/π), 0]) : nothing
	end
	f₁ ./= sum(f₁*Δx)
	f₂ = 0.3*exp.(-(x.-16).^2 ./0.9) + 0.7*exp.(-(x.-18).^2 ./0.4)
	f₂ ./= sum(f₂*Δx)
	g = 0.4*exp.(-(x.-10.5).^2 ./1.6) + 0.6*exp.(-(x.-13).^2 ./1.3)
	g ./= sum(g*Δx)
	f = [f₁, f₂]
end;

# ╔═╡ fcb46c03-0802-4ee5-bb6f-d09f82991df4
begin
	p = collect(0:0.0001:1)
	M = length(p)
	F = [zeros(N) for i in 1:n]
	F⁻¹ = [zeros(M) for i in 1:n]
	B = zeros(N)
	b = zeros(N)
end;

# ╔═╡ 5a2e20c6-8e14-4984-a607-553d52b47c23
md"""
$B_{L^2}(Λ; f_1, \dots, f_n) = {\arg\min}_{g \in \mathcal P \cap L^2(\Omega)} \sum_{i=1}^n \lambda_i \Vert g - f_i \Vert^2_{L^2(\Omega)} = \sum_{i=1}^n \lambda_i f_i$
"""

# ╔═╡ 3bbd9c53-f44e-4948-9017-c4a2e6cc8fd5
@bind λ_L2 Slider(0:0.01:1; default=0.0, show_value=true)

# ╔═╡ fe47bb3e-12b7-4eed-86a4-ba1b01d67eac
plot(x, [f, (1-λ_L2)*f[1]+λ_L2*f[2]], ylimits=[0,1], label=["f₁" "f₂" "λ₁f₁ + λ₂f₂"], linewidth=[1 1 3], legendfontsize=10)

# ╔═╡ 695ffb1b-7bfe-4f05-a5cf-b7b32c784bab
md"""
$B_{W^2}(Λ; f_1, \dots, f_n) = {\arg\min}_{g \in \mathcal P(\Omega)} \sum_{i=1}^n \lambda_i W_2(g, f_i)^2$
"""

# ╔═╡ 7a3a25b1-572b-40c9-a522-00cfa92af6b0
@bind λ Slider(0:0.01:1; default=0.0, show_value=true)

# ╔═╡ 9744811c-fdfb-4697-b6d9-3232d795438e
begin
	for i in 1:n
		MorporJ.cdf!(F[i], f[i], Δx)
		MorporJ.icdf!(F⁻¹[i], F[i], x, p)
	end
	B⁻¹ = (1-λ)*F⁻¹[1] + λ*F⁻¹[2]
	MorporJ.iicdf!(B, b, B⁻¹, x, p, 1e-12)
	for i in eachindex(b)
		b[i] < 0 ? b[i] = 0 : nothing
	end
end;

# ╔═╡ e7c0fe22-4baf-4530-a63a-9ad211f26a8d
md"Fit $(@bind fit CheckBox())"

# ╔═╡ 2643241e-c45a-4bd7-96a6-39ab56352db2
intfit = Int(fit);

# ╔═╡ 96420763-da4e-47de-9fbd-8ba6e24d655f
plot(x, [f, b, g], ylimits=[0,(λ+1-λ)], label=["f₁" "f₂" "(λ₁id + λ₂T)#f₁" "g"], linewidth=[1 1 3 2*intfit], legendfontsize=10)

# ╔═╡ 0c676bce-3964-43b8-96c6-b9e54b5db9f5
md"""
$Λ^*(g; f_1, \dots, f_n) \in {\arg\min}_{Λ \in Σ_n} \; \text{Loss} \big (g, B_{W^2}(Λ; f_1, \dots, f_n) \big)$
$g \approx B_{W^2}(Λ^*(g); f_1, \dots, f_n)$
"""

# ╔═╡ 4e6e0002-5703-4b03-8e92-c6f9a43b7100
md"""
- How does one find well-suited $f_1, \dots, f_n$?
- If $g = g(\mu)$, what are the properties of $\mu \mapsto Λ^* \big (g(\mu) \big )$?
- How many does the (expected, maximal) approximation error change with $n$?
"""

# ╔═╡ Cell order:
# ╟─37bf7c73-acda-42ee-acce-e10131f83e12
# ╟─1e8a18a4-c2ba-482b-9624-fdafc49521c9
# ╠═c3dd7a97-aea6-4f58-80a4-175803b55167
# ╟─fcb46c03-0802-4ee5-bb6f-d09f82991df4
# ╟─9744811c-fdfb-4697-b6d9-3232d795438e
# ╟─5a2e20c6-8e14-4984-a607-553d52b47c23
# ╟─fe47bb3e-12b7-4eed-86a4-ba1b01d67eac
# ╟─3bbd9c53-f44e-4948-9017-c4a2e6cc8fd5
# ╟─695ffb1b-7bfe-4f05-a5cf-b7b32c784bab
# ╟─96420763-da4e-47de-9fbd-8ba6e24d655f
# ╟─7a3a25b1-572b-40c9-a522-00cfa92af6b0
# ╟─e7c0fe22-4baf-4530-a63a-9ad211f26a8d
# ╟─2643241e-c45a-4bd7-96a6-39ab56352db2
# ╟─0c676bce-3964-43b8-96c6-b9e54b5db9f5
# ╟─4e6e0002-5703-4b03-8e92-c6f9a43b7100
