"""
barycenter_fit(t, a, Δp)

Given a k-vector of target icdfs t and a n-vector of atoms a, finds a vector of weights 

Λᵢʲ = [λ₁, …, λₙ]ʲ : 1 ≤ j ≤ k

such that 

‖ tⱼ - ∑ᵢ aᵢ λᵢʲ ‖₂ = ∑ⱼ W₂²(tⱼ, B(λ₁ʲ, …, λₙʲ; a₁, …, aₙ))  is minimal ∀j.

Returns the weight vector Λ and the fitting W₂ errors ΔW
"""

function barycenter_fit(t, a, Δp, max_iter=1e7, eps_abs=1e-9, eps_rel=1e-6)
    
    n = length(a)
    k = length(t)
    Λ = [zeros(n) for _ in 1:k]
    ΔW = zeros(k)

    M = zeros(n, n)

    # Assemble mass matrix
    for i in 1:n
        for j in i:n
            M[i,j] = dot(a[i], a[j])
            M[j,i] = M[i,j]
        end
    end
    # linear term
    c = [dot(a[i], t[1]) for i in eachindex(a)] 

    #set up the problem
    x = Variable(n)

    x₀ = [1/norm(a[i] - t[1]) for i in 1:n]
    x₀ ./= sum(x₀)

    # if the targets are somewhat ordered, last optimal x may also be a good starting point
    set_value!(x, x₀)
    #set_value!(x, [1/n for _ in 1:n]) # initial guess is λᵢ ≡ 1/n ∀ i

    expression = quadform(x, M) - 2*dot(x,c)
    problem = minimize( expression , [x >= 0, sum(x) == 1])

    solver = () -> COSMO.Optimizer( verbose=false,
                                    eps_abs = eps_abs,
                                    eps_rel = eps_rel,
                                    max_iter = max_iter)

    # check if target is in dictionary
    skip = false
    for i in 1:n
        if t[1] == a[i]
            Λ[1][i] = 1
            ΔW[1] = 0
            skip = true
        end
    end
    if !skip
        solve!(problem, solver; warmstart=true)
    
        # optimal barycentric weights
        Λ[1] .= evaluate(x)
        # W₂ distance between optimal barycenter and target
        ΔW[1] = Δp * ( evaluate( expression ) + dot(t[1],t[1]) )
    end

    # fit the rest of the targets
    for j in 2:k
        # new linear term
        for i in eachindex(c)
            c[i] = dot(a[i], t[j])
        end

        skip = false
        for i in 1:n
            if t[j] == a[i]
                Λ[j][i] = 1
                ΔW[j] = 0
                skip = true
            end
        end
        if !skip

            x₀ .= [1/norm(a[i] - t[j]) for i in 1:n]
            x₀ ./= sum(x₀)

            # if the targets are somewhat ordered, last optimal x may also be a good starting point
            set_value!(x, x₀)
            solve!(problem, solver; warmstart=true)

            # if not converged, set all lambda equal
            if !(Int(problem.status) in [1])
                # set_value!(x, [1/n for _ in 1:n])
            end
            
            Λ[j] .= evaluate(x)
            ΔW[j] = Δp * (evaluate(expression) + dot(t[j],t[j]) )
        end
    end

    return Λ, sqrt.(abs.(ΔW))
end