"""
barycenter_fit(t, a, Δp)

Given a k-vector of target icdfs t and a n-vector of atoms a, finds a vector of weights 

Λᵢʲ = [λ₁, …, λₙ]ʲ : 1 ≤ j ≤ k

such that 

‖ tⱼ - ∑ᵢ aᵢ λᵢʲ ‖₂ = ∑ⱼ W₂²(tⱼ, B(λ₁ʲ, …, λₙʲ; a₁, …, aₙ))  is minimal ∀j.

Returns the weight vector Λ and the fitting W₂ errors ΔW
"""

function barycenter_fit(t, a, Δp; max_iter=50000)
    
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
    c = [dot(a[i], t[1]) for i in eachindex(a)] # initial guess is λᵢ ≡ 1/n ∀ i

    #set up the problem
    x = Variable(n)
    set_value!(x, [1/n for _ in 1:n]) # initial guess 

    expression = quadform(x, M) - 2*dot(x,c)
    problem = minimize( expression , [x >= 0, sum(x) == 1])

    solver = () -> COSMO.Optimizer( verbose=false,
                                    eps_abs = 1e-6,
                                    eps_rel = 1e-5,
                                    max_iter = max_iter)
    solve!(problem, solver; warmstart=true)
    
    # optimal barycentric weights
    Λ[1] .= evaluate(x)
    # W₂ distance between optimal barycenter and target
    ΔW[1] = Δp * ( evaluate( expression ) + dot(t[1],t[1]) )

    # fit the rest of the targets
    for j in 2:k
        # new linear term
        for i in eachindex(c)
            c[i] = dot(a[i], t[j])
        end

        # if the targets are somewhat ordered, last optimal x may also be a good starting point
        set_value!(x, [1/n for _ in 1:n])
        solve!(problem, solver; warmstart=true)

        if !(Int(problem.status) in [1])
            set_value!(x, [1/n for _ in 1:n])
        end
        
        Λ[j] .= evaluate(x)
        ΔW[j] = Δp * (evaluate(expression) + dot(t[j],t[j]) )
    end

    return Λ, sqrt.(abs.(ΔW))
end

function barycenter_fit_pgrid_projected(t, a, p_t, p_a, tol=1e-9)
    
    n = length(a)
    k = length(t)
    Λ = [zeros(n) for _ in 1:k]
    ΔW = zeros(k)

    M = zeros(n, n)

    at = copy(a)
    p_at = copy(p_a)
    push!(at, t[1])
    push!(p_at, p_t[1])

    at_comm, p_at_comm = iicdf_common_grid(at, p_at, tol)
    m = length(p_at_comm)
    # Assemble mass matrix
    for i in 1:n
        for j in i:n
            M[i,j] = 0
            for l in 1:(m-1)
                M[i,j] += (p_at_comm[l+1]-p_at_comm[l])*0.5*(
                            at_comm[i+1][l]*at_comm[j+1][l] + at_comm[i+1][l+1]*at_comm[j+1][l+1])
            end
            M[j,i] = M[i,j]
        end
    end
    # linear term
    c = zeros(n)
    for i in 1:n
        c[i] = 0
        for l in 1:(m-1)
            c[i] += (p_at_comm[l+1]-p_at_comm[l])*0.5*(
                        at_comm[i+1][l]*at_comm[1][l] + at_comm[i+1][l+1]*at_comm[1][l+1])
        end
    end
    # constant term
    tt = 0
    for l in 1:(m-1)
        tt += (p_at_comm[l+1]-p_at_comm[l])*0.5*(
                    at_comm[1][l]*at_comm[1][l] + at_comm[1][l+1]*at_comm[1][l+1])
    end

    #set up the problem
    x = Variable(n)
    set_value!(x, [1/n for _ in 1:n])

    expression = quadform(x, M) - 2*dot(x,c)
    problem = minimize( expression , [x >= 0, sum(x) == 1])

    solver = () -> COSMO.Optimizer(verbose=false, max_iter = 100000)
    solve!(problem, solver; warmstart=true)
    
    # optimal barycentric weights
    Λ[1] .= evaluate(x)
    # W₂ distance between optimal barycenter and target
    ΔW[1] = evaluate( expression ) + tt

    # fit the rest of the targets
    for j in 2:k

        pop!(at)
        pop!(p_at)
        push!(at, t[j])
        push!(p_at, p_t[j])

        at_comm, p_at_comm = iicdf_common_grid(at, p_at, tol)

        m = length(p_at_comm)
        # Assemble mass matrix
        for i in 1:n
            for j in i:n
                M[i,j] = 0
                for l in 1:(m-1)
                    M[i,j] += (p_at_comm[l+1]-p_at_comm[l])*0.5*(
                                at_comm[i+1][l]*at_comm[j+1][l] + at_comm[i+1][l+1]*at_comm[j+1][l+1])
                end
                M[j,i] = M[i,j]
            end
        end
        # linear term
        for i in 1:n
            c[i] = 0
            for l in 1:(m-1)
                c[i] += (p_at_comm[l+1]-p_at_comm[l])*0.5*(
                            at_comm[i+1][l]*at_comm[1][l] + at_comm[i+1][l+1]*at_comm[1][l+1])
            end
        end
        # constant term
        tt = 0
        for l in 1:(m-1)
            tt += (p_at_comm[l+1]-p_at_comm[l])*0.5*(
                        at_comm[1][l]*at_comm[1][l] + at_comm[1][l+1]*at_comm[1][l+1])
        end

        # if the targets are somewhat ordered, last optimal x is a good starting point
        solve!(problem, solver; warmstart=true)

        if !(Int(problem.status) in [1])
            set_value!(x, [1/n for _ in 1:n])
        end
        
        Λ[j] .= evaluate(x)
        ΔW[j] = evaluate(expression) + tt
    end

    return Λ, ΔW
end