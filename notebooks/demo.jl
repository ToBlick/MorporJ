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

# ╔═╡ 5a2e20c6-8e14-4984-a607-553d52b47c23
md"""
$B_{L^2}(Λ; f_1, \dots, f_n) = {\arg\min}_{g \in \mathcal P \cap L^2(\Omega)} \sum_{i=1}^n \lambda_i \Vert g - f_i \Vert^2_{L^2(\Omega)} = \sum_{i=1}^n \lambda_i f_i$
"""

# ╔═╡ 3bbd9c53-f44e-4948-9017-c4a2e6cc8fd5
md""" λ = $(@bind λ_L2 Slider(0:0.01:1; default=0.0, show_value=true))"""

# ╔═╡ 695ffb1b-7bfe-4f05-a5cf-b7b32c784bab
md"""
$B_{W^2}(Λ; f_1, \dots, f_n) = {\arg\min}_{g \in \mathcal P(\Omega)} \sum_{i=1}^n \lambda_i W_2(g, f_i)^2$
"""

# ╔═╡ 7a3a25b1-572b-40c9-a522-00cfa92af6b0
md""" λ = $(@bind λ Slider(0:0.01:1; default=0.0, show_value=true))"""

# ╔═╡ e7c0fe22-4baf-4530-a63a-9ad211f26a8d
md"Fit $(@bind fit CheckBox())"

# ╔═╡ 0c676bce-3964-43b8-96c6-b9e54b5db9f5
md"""
$Λ^*(g; f_1, \dots, f_n) \in {\arg\min}_{Λ \in Σ_n} \; \text{Loss} \big (g, B_{W^2}(Λ; f_1, \dots, f_n) \big)$
$g \approx B_{W^2}(Λ^*(g); f_1, \dots, f_n)$
"""

# ╔═╡ 9ce7e85d-6d39-4cad-8146-456e97d01e0e
md"""λ₁ = $(@bind λ₁ Slider(0:0.01:1; default=0.0, show_value=true))"""

# ╔═╡ ad9dfabf-4dbb-4b32-93cf-94782fd294c6
md"""λ₂ = $(@bind λ₂ Slider(0:0.01:1; default=0.0, show_value=true))"""

# ╔═╡ 4e6e0002-5703-4b03-8e92-c6f9a43b7100
md"""
- How does one find well-suited $f_1, \dots, f_n$?
- If $g = g(\mu)$, what are the properties of $\mu \mapsto Λ^* \big (g(\mu) \big )$?
- How does the (expected, maximal) approximation error change with $n$?
"""

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
	g = 0.4*exp.(-(x.-10.5).^2 ./1.6) + 0.6*exp.(-(x.-12.5).^2 ./1.3)
	g ./= sum(g*Δx)
	f = [f₁, f₂]
end;

# ╔═╡ fe47bb3e-12b7-4eed-86a4-ba1b01d67eac
plot(x, [f, (1-λ_L2)*f[1]+λ_L2*f[2]], ylimits=[0,1], label=["f₁" "f₂" "λ₁f₁ + λ₂f₂"], linewidth=[1 1 3], legendfontsize=10)

# ╔═╡ fcb46c03-0802-4ee5-bb6f-d09f82991df4
begin
	p = collect(0:0.0001:1)
	M = length(p)
	F = [zeros(N) for i in 1:n]
	F⁻¹ = [zeros(M) for i in 1:n]
	B = zeros(N)
	b = zeros(N)
end;

# ╔═╡ 9744811c-fdfb-4697-b6d9-3232d795438e
begin
	for i in 1:n
		MorporJ.cdf!(F[i], f[i], Δx)
		MorporJ.icdf!(F⁻¹[i], F[i], x, p)
	end
	B⁻¹ = (1-λ)*F⁻¹[1] + λ*F⁻¹[2]
	MorporJ.iicdf!(B, B⁻¹, x, p, 1e-12)
	MorporJ.cdf_to_pdf!(b, B, Δx; order=1)
	for i in eachindex(b)
		b[i] < 0 ? b[i] = 0 : nothing
	end
end;

# ╔═╡ 2643241e-c45a-4bd7-96a6-39ab56352db2
intfit = Int(fit);

# ╔═╡ 96420763-da4e-47de-9fbd-8ba6e24d655f
plot(x, [f, b, g], ylimits=[0,(λ+1-λ)], label=["f₁" "f₂" "(λ₁id + λ₂T)#f₁" "g"], linewidth=[1 1 3 2*intfit], legendfontsize=10)

# ╔═╡ 41b82709-19d5-489c-8486-13fde19eece2
begin
	f₃ = zeros(N)
	for i in eachindex(x)
		9 < x[i] < 12 ? f₃[i] = 1 : nothing
	end
	f₃ ./= sum(f₃*Δx)
	fₑ = [ f₁, f₂, f₃ ]
	nₑ = 3
	Fₑ = [zeros(N) for i in 1:nₑ]
	F⁻¹ₑ = [zeros(M) for i in 1:nₑ]
	Bₑ = zeros(N)
	bₑ = zeros(N)
end;

# ╔═╡ 19f912c0-5ec4-4d80-8977-7c7a12aa90d2
plot(x, [fₑ, bₑ], ylimits=[0,1+0*(λ₁+λ₂)], label=["f₁" "f₂" "f₃" "B(Λ)" "g"], linewidth=[1 1 1 3], legend = :topleft, legendfontsize=10)

# ╔═╡ 068869f1-08ed-4eaf-88d5-e595ae2028ca
begin
	for i in 1:nₑ
		MorporJ.cdf!(Fₑ[i], fₑ[i], Δx)
		MorporJ.icdf!(F⁻¹ₑ[i], Fₑ[i], x, p)
	end
	λ₃ = 1-λ₁-λ₂
	B⁻¹ₑ = λ₁*F⁻¹ₑ[1] + λ₂*F⁻¹ₑ[2] + λ₃*F⁻¹ₑ[3]
	MorporJ.iicdf!(Bₑ, B⁻¹ₑ, x, p, 1e-12)
	MorporJ.cdf_to_pdf!(bₑ, Bₑ, Δx; order=1)
	for i in eachindex(bₑ)
		bₑ[i] < 0 ? bₑ[i] = 0 : nothing
	end
end;

# ╔═╡ 44105627-7662-4ada-9a81-d9f232f41c7e
md"""λ₃ = $(round(λ₃,digits=2))"""

# ╔═╡ Cell order:
# ╟─5a2e20c6-8e14-4984-a607-553d52b47c23
# ╠═fe47bb3e-12b7-4eed-86a4-ba1b01d67eac
# ╟─3bbd9c53-f44e-4948-9017-c4a2e6cc8fd5
# ╟─695ffb1b-7bfe-4f05-a5cf-b7b32c784bab
# ╠═96420763-da4e-47de-9fbd-8ba6e24d655f
# ╟─7a3a25b1-572b-40c9-a522-00cfa92af6b0
# ╟─e7c0fe22-4baf-4530-a63a-9ad211f26a8d
# ╟─0c676bce-3964-43b8-96c6-b9e54b5db9f5
# ╠═19f912c0-5ec4-4d80-8977-7c7a12aa90d2
# ╟─9ce7e85d-6d39-4cad-8146-456e97d01e0e
# ╟─ad9dfabf-4dbb-4b32-93cf-94782fd294c6
# ╟─44105627-7662-4ada-9a81-d9f232f41c7e
# ╟─4e6e0002-5703-4b03-8e92-c6f9a43b7100
# ╟─37bf7c73-acda-42ee-acce-e10131f83e12
# ╟─1e8a18a4-c2ba-482b-9624-fdafc49521c9
# ╟─c3dd7a97-aea6-4f58-80a4-175803b55167
# ╟─fcb46c03-0802-4ee5-bb6f-d09f82991df4
# ╟─9744811c-fdfb-4697-b6d9-3232d795438e
# ╟─2643241e-c45a-4bd7-96a6-39ab56352db2
# ╟─41b82709-19d5-489c-8486-13fde19eece2
# ╟─068869f1-08ed-4eaf-88d5-e595ae2028ca
