"""
icdf_exact(F, xgrid, tol=1e-12)

Given a discrete cdf F, returns its discrete icdf F⁻¹

F⁻¹ : [0,1] ∋ p ↦ inf { x: F(x) ≥ p } ∈ [0,L]

F⁻¹ will be defined on a non-uniform grid which is returned aswell
"""
function icdf_exact(F, xgrid, tol=1e-12)
    pgrid = []
    F⁻¹ = []
    for i in eachindex(F)
        if i == 1 || (F[i]-F[i-1]>tol)
            push!(pgrid, F[i])
            push!(F⁻¹, xgrid[i])
        end
    end
    return F⁻¹, pgrid
end

function iicdf_exact(F⁻¹, pgrid, xgrid, tol=1e-12)
    F = zero(xgrid)
    M = length(pgrid)
    j = 1 # index in the pgrid
    for i in eachindex(xgrid) # index in the xgrid
        x = xgrid[i]
        while F⁻¹[j] <= x && j < M
            j += 1
        end
        if j == 1 # F⁻¹[1] > x : icdf does not start at x=0
            F[i] = pgrid[1]
        elseif j == M && F⁻¹[j] <= x # F⁻¹[M] =: l <= x : icdf has reached the last point of support
            F[i] = pgrid[M]
        elseif pgrid[j] - pgrid[j-1] < tol
            F[i] = pgrid[j-1] 
        else
            # F(x) =  p₋ + (x - x₋) * (p₊ - p₋) / (x₊ - x₋)
            F[i] = pgrid[j-1] + (x - F⁻¹[j-1]) * (pgrid[j] - pgrid[j-1]) / (F⁻¹[j] - F⁻¹[j-1])
        end
    end
    return F
end

function interpolate_to_grid(F⁻¹, pgrid, new_pgrid, tol=1e-12)

    m = length(pgrid)
    M = length(new_pgrid)
    F⁻¹_new = zero(new_pgrid)

    j = 1
    for i in 1:M
        p = new_pgrid[i]

        while pgrid[j] < p && j < m
            j += 1
        end
        # now pⱼ >= p
        if j == 1
            if length(pgrid) == 1 # dirac delta
                F⁻¹_new[i] = F⁻¹[j]
                continue
            else
                p₋ = pgrid[j]
                p₊ = pgrid[j+1]
                F⁻¹₋ = F⁻¹[j]
                F⁻¹₊ = F⁻¹[j+1]
            end
        else
            p₋ = pgrid[j-1]
            p₊ = pgrid[j]
            F⁻¹₋ = F⁻¹[j-1]
            F⁻¹₊ = F⁻¹[j]
        end
        if p₊ - p₋ < tol
            F⁻¹_new[i] = F⁻¹₋
        else
            F⁻¹_new[i] = F⁻¹₋ + (p - p₋) * (F⁻¹₊ - F⁻¹₋) / (p₊ - p₋)
        end
    end
    return F⁻¹_new    
end

function iicdf_common_grid(F⁻¹_v, pgrid_v, tol=1e-12)

    n = length(F⁻¹_v)

    M = 0
    for i in 1:n
        M += length(pgrid_v[i])
    end
    
    ΣF⁻¹_v = [ [] for _ in 1:n]
    Σpgrid = []
    i_vec = [0 for _ in 1:n]
    p₀ = 2

    for j in 1:M
        p₀₋ = p₀ # last p₀
        p₀ = 2
        iₐ = 0 # active index
        for i in 1:n # find next smallest pgrid value
            if i_vec[i] == length(pgrid_v[i])
                continue
            elseif pgrid_v[i][i_vec[i]+1] < p₀
                p₀ = pgrid_v[i][i_vec[i]+1]
                iₐ = i
            end
        end
        if iₐ == 0 # finished
            break
        end

        # update interpolation point for active function
        i_vec[iₐ] += 1

        if p₀ - p₀₋ < tol && j != 1 # if new p₀ is identical with a previous one, no need to duplicate
            continue
        end

        # add p-grid value to the grid union
        push!(Σpgrid, p₀)
        #Σpgrid[j] = p₀
        for i in 1:n
            # calculate interpolant value
            if i_vec[i] == 0 # this icdf is not yet different from zero
                push!(ΣF⁻¹_v[i], 0.0)
                #ΣF⁻¹_v[i][j] = 0
            elseif i_vec[i] == 1 # first icdf value different from zero
                if length(pgrid_v[i]) == 1 # dirac delta
                    push!(ΣF⁻¹_v[i], F⁻¹_v[i][i_vec[i]])
                    #ΣF⁻¹_v[i][j] = F⁻¹_v[i][i_vec[i]]
                else
                    p₋ = pgrid_v[i][i_vec[i]]
                    p₊ = pgrid_v[i][i_vec[i]+1]
                    F⁻¹₋ = F⁻¹_v[i][i_vec[i]]
                    F⁻¹₊ = F⁻¹_v[i][i_vec[i]+1]
                    F⁻¹_i = F⁻¹₋ + (p₀ - p₋) * (F⁻¹₊ - F⁻¹₋) / (p₊ - p₋)
                    push!(ΣF⁻¹_v[i], F⁻¹_i)
                    #ΣF⁻¹_v[i][j] = F⁻¹_i
                end
            else # linear interpolation
                p₋ = pgrid_v[i][i_vec[i]-1]
                p₊ = pgrid_v[i][i_vec[i]]
                F⁻¹₋ = F⁻¹_v[i][i_vec[i]-1]
                F⁻¹₊ = F⁻¹_v[i][i_vec[i]]
                F⁻¹_i = F⁻¹₋ + (p₀ - p₋) * (F⁻¹₊ - F⁻¹₋) / (p₊ - p₋)
                push!(ΣF⁻¹_v[i], F⁻¹_i)   
                #ΣF⁻¹_v[i][j] = F⁻¹_i        
            end
        end
    end
    return ΣF⁻¹_v, Σpgrid
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