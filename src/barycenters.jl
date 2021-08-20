"""
barycenter_fit(t, a, Δx)

Given a k-vector of target icdfs t and a n-vector of atoms a, finds a vector of weights 

Λᵢʲ = [λ₁, …, λₙ]ʲ : 1 ≤ j ≤ k

such that 

‖ tⱼ - ∑ᵢ aᵢ λᵢʲ ‖₂ = ∑ⱼ W₂²(tⱼ, B(λ₁ʲ, …, λₙʲ; a₁, …, aₙ))  is minimal ∀j.

Returns the weight vector Λ and the fitting W₂ errors Δ
"""

function barycenter_fit(t, a, Δx)
    
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
    # t_ = [dot(t[1],t[1])]

    #set up the problem
    x = Variable(n)
    set_value!(x, [1/n for _ in 1:n])

    expression = quadform(x, M) - 2*dot(x,c) # + t_[1]
    problem = minimize( expression , [x >= 0, sum(x) == 1])

    solver = () -> COSMO.Optimizer(verbose=false, max_iter = 100000)
    solve!(problem, solver; warmstart=true)
    
    # optimal barycentric weights
    Λ[1] .= evaluate(x)
    # W₂ distance between optimal barycenter and target
    ΔW[1] = Δx * evaluate( expression )

    # fit the rest of the targets
    for j in 2:k
        # new linear term
        for i in eachindex(c)
            c[i] = dot(a[i], t[j])
        end
        # t_ .= dot(t[j],t[j])

        # if the targets are somewhat ordered, last optimal x is a good starting point
        solve!(problem, solver; warmstart=true)

        if !(Int(problem.status) in [1])
            set_value!(x, [1/n for _ in 1:n])
        end
        
        Λ[j] .= evaluate(x)
        ΔW[j] = Δx * (evaluate(expression) + dot(t[j],t[j]) )
    end

    return Λ, ΔW
end