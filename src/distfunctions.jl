"""
cdf!(f, F)

Given a discrete pdf f, writes its cdf into F. The value in F[i] includes the integral of the ith cell.

Limited to length(f) == length(xgrid)-1 == length(F)
"""
function cdf!(F, f, Δx, tol = 1e-12)
    F[1] = f[1] * Δx
    for i in 2:length(f)
        F[i] = F[i-1] + f[i] * Δx
        if 1 - F[i] < tol
            F[i] = 1
        end
    end
end

"""
icdf!(F⁻¹, F, xgrid, pgrid)

Given a discrete cdf F, writes its discrete icdf into F⁻¹

F⁻¹ : [0,1] ∋ p ↦ inf { x: F(x) ≥ p } ∈ [0,L]

Limited to length(F) == length(F⁻¹), uniform grids
"""
function icdf!(F⁻¹, F, xgrid, pgrid, tol = 1e-9)
    N = length(F)
    M = length(F⁻¹)
    i = 1
    for j in 1:M
        p = pgrid[j]
        while F[i] < p && i < N
            i += 1
        end
        if i == 1 # p₀ := F[1] ≥ p[j] ⇒ cdf does not start at zero
            # return smallest x
            F⁻¹[j] = xgrid[1]
            continue
        # F[i] ≥ p[j] from here on
        elseif  F[i] < p && i == N
            # F[end] < p[j] ⇒ cdf does not reach 1
            # return largest x
            F⁻¹[j] = xgrid[end]
            continue
        else
            #quadratic interpolation
            if i == N #|| (1 - F[i]) < tol # the p[j] = 1 case: F[i], F[i+1], … will be 1.
                p₋ = F[i-2]
                p₊ = F[i-1]
                p₊₊ = F[i]
                x₋ = xgrid[i-2]
                x₊ = xgrid[i-1]
                x₊₊ = xgrid[i]
            #= elseif (1 - F[i+1]) < tol # avoid interpolationg over the kink in the cdf at p = 1
                p₋ = F[i-2]
                p₊ = F[i-1]
                p₊₊ = F[i]
                x₋ = xgrid[i-2]
                x₊ = xgrid[i-1]
                x₊₊ = xgrid[i] =#
            else
                p₋ = F[i-1]
                p₊ = F[i]
                p₊₊ = F[i+1]
                x₋ = xgrid[i-1]
                x₊ = xgrid[i]
                x₊₊ = xgrid[i+1]
            end

            if abs(p₊-p₋) < tol
                F⁻¹[j] = x₋
                continue
            elseif abs(p₊₊-p₊) < tol
                if i > 2
                    p₋ = F[i-2]
                    p₊ = F[i-1]
                    p₊₊ = F[i]
                    x₋ = xgrid[i-2]
                    x₊ = xgrid[i-1]
                    x₊₊ = xgrid[i]
                else
                    F⁻¹[j] = x₋
                    continue
                end
            end

            # interpolant of p ↦ F⁻¹(p)
            F⁻¹[j] = ( x₋ * (p-p₊)/(p₋-p₊) * (p-p₊₊)/(p₋-p₊₊) 
                    + x₊ * (p-p₋)/(p₊-p₋) * (p-p₊₊)/(p₊-p₊₊) 
                    + x₊₊ * (p-p₋)/(p₊₊-p₋) * (p-p₊)/(p₊₊-p₊) )

            #=
            # linear interpolation
            x₋ = xgrid[i-1]
            x₊ = xgrid[i] 
            p₋ = F[i-1]
            p₊ = F[i]

            if isapprox(p₊, p₋)
                F⁻¹[j] =  x₋ # see the inf in the definition
            else
                F⁻¹[j] =  x₋ + (p - p₋) * (x₊ - x₋) / (p₊ - p₋)
            end
            =#
            if F⁻¹[j] < xgrid[1]
                F⁻¹[j] = xgrid[1]
            elseif F⁻¹[j] > xgrid[end]
                F⁻¹[j] = xgrid[end]
            end
        end
    end
end

"""
iicdf!(F, F⁻¹, xgrid, pgrid)

Given a discrete icdf F⁻¹, writes its discrete cdf into F and its discrete pdf into f

F⁻¹ : [0,1] ∋ p ↦ inf { x: F(x) ≥ p } ∈ [0,L]

Limited to length(F) == length(F⁻¹), uniform grids
"""
function iicdf!(F, f, F⁻¹, xgrid, pgrid, tol=1e-9)
    N = length(F)
    M = length(F⁻¹)
    j = 1
    for i in 1:N
        x = xgrid[i]
        while F⁻¹[j] <= x && j < M
            j += 1
        end
        # now either F⁻¹[j] > x[i] or there exists no j: F⁻¹[M] <= x[i]
        if j == 1 # x₀ := F⁻¹[1] = F⁻¹(0) > x[i] ⇒ icdf does not start at x[1]
            # return smallest p
            F[i] = pgrid[1]
            f[i] = 0
            continue
        elseif j == M && F⁻¹[j] <= x
            # F⁻¹[end] = F⁻¹(1) <= x[i] ⇒ icdf does not reach x[end] ⇒ fill up with ones
            F[i] = pgrid[end]
            if i == 1
                f[i] = F[i] / (xgrid[2] - xgrid[1])
            else
                f[i] = (F[i] - F[i-1]) / (xgrid[i] - xgrid[i-1])
            end
            continue
        else
            #quadratic interpolation
            if i == 1 && j <= (M-2) # kink at x = 0 in the cdf
                p₋ = pgrid[j]
                p₊ = pgrid[j+1]
                p₊₊ = pgrid[j+2]
                x₋ = F⁻¹[j]
                x₊ = F⁻¹[j+1]
                x₊₊ = F⁻¹[j+2]
                #= F[i] = p₋ + (x - x₋) * (p₊ - p₋) / (x₊ - x₋) 
                f[i] = (p₊ - p₋) / (x₊ - x₋) 
                if F[i] < pgrid[1]
                    F[i] = pgrid[1]
                elseif F[i] > pgrid[end]
                    F[i] = pgrid[end]
                end
                continue =#
            elseif j == M
                p₋ = pgrid[j-2]
                p₊ = pgrid[j-1]
                p₊₊ = pgrid[j]
                x₋ = F⁻¹[j-2]
                x₊ = F⁻¹[j-1]
                x₊₊ = F⁻¹[j]
            #=elseif j == (M-1)
                p₋ = pgrid[j-2]
                p₊ = pgrid[j-1]
                p₊₊ = pgrid[j]abs
                x₋ = F⁻¹[j-2]
                x₊ = F⁻¹[j-1]
                x₊₊ = F⁻¹[j]=#
            else
                p₋ = pgrid[j-1]
                p₊ = pgrid[j]
                p₊₊ = pgrid[j+1]
                x₋ = F⁻¹[j-1]
                x₊ = F⁻¹[j]
                x₊₊ = F⁻¹[j+1]
            end

            if abs(x₊-x₋) < tol
                F[i] = p₊₊
                if i == 1
                    f[i] = F[i] / (xgrid[2] - xgrid[1])
                else
                    f[i] = (F[i] - F[i-1]) / (xgrid[i] - xgrid[i-1])
                end
                continue
            elseif abs(x₊₊-x₊) < tol
                if j > 2
                    p₋ = pgrid[j-2]
                    p₊ = pgrid[j-1]
                    p₊₊ = pgrid[j]
                    x₋ = F⁻¹[j-2]
                    x₊ = F⁻¹[j-1]
                    x₊₊ = F⁻¹[j]
                else
                    F[i] = p₊₊
                    if i == 1
                        f[i] = F[i] / (xgrid[2] - xgrid[1])
                    else
                        f[i] = (F[i] - F[i-1]) / (xgrid[i] - xgrid[i-1])
                    end
                    continue
                end
            end

            # interpolant of x ↦ F(x)
            F[i] = ( p₋ * (x-x₊)/(x₋-x₊) * (x-x₊₊)/(x₋-x₊₊) 
                    + p₊ * (x-x₋)/(x₊-x₋) * (x-x₊₊)/(x₊-x₊₊) 
                    + p₊₊ * (x-x₋)/(x₊₊-x₋) * (x-x₊)/(x₊₊-x₊) )

            f[i] = ( p₋ * 1/(x₋-x₊) * (x-x₊₊)/(x₋-x₊₊)
                    + p₋ * (x-x₊)/(x₋-x₊) * 1/(x₋-x₊₊) 
                    + p₊ * 1/(x₊-x₋) * (x-x₊₊)/(x₊-x₊₊)
                    + p₊ * (x-x₋)/(x₊-x₋) * 1/(x₊-x₊₊)
                    + p₊₊ * 1/(x₊₊-x₋) * (x-x₊)/(x₊₊-x₊)
                    + p₊₊ * (x-x₋)/(x₊₊-x₋) * 1/(x₊₊-x₊) )


            #= 
            I_icdf = y -> ( x₋*(y-p₊)/(p₋-p₊)*(y-p₊₊)/(p₋-p₊₊) + x₊*(y-p₋)/(p₊-p₋)*(y-p₊₊)/(p₊-p₊₊) 
                            + x₊₊*(y-p₋)/(p₊₊-p₋)*(y-p₊)/(p₊₊-p₊) )

            if j == N
                upper = p₊₊
            else
                upper = p₊
            end
            lower = p₋
            F[i] = 0.5*(upper + lower)
            ϵ = I_icdf(F[i]) - x
            for k in 1:100
                if ϵ > 0 # I_icdf(F[i]) > x
                    upper = F[i]
                else
                    lower = F[i]
                end
                F[i] = 0.5*(F[i] + upper)
                ϵ = I_icdf(F[i]) - x
                if abs(ϵ) < tol
                    break
                end
            end 
            =#

            # fix values outside of allowed domain
            if F[i] < pgrid[1]
                F[i] = pgrid[1]
            elseif F[i] > pgrid[end]
                F[i] = pgrid[end]
            end

            #=
            if i == 1
                f[1] = F[1] / (xgrid[2] - xgrid[1])
            else
                ∂xI_icdf_i = ( x₋*(F[i]-p₊)/(p₋-p₊)/(p₋-p₊₊) + x₋/(p₋-p₊)*(F[i]-p₊₊)/(p₋-p₊₊)
                           + x₊*(F[i]-p₋)/(p₊-p₋)/(p₊-p₊₊) + x₊/(p₊-p₋)*(F[i]-p₊₊)/(p₊-p₊₊) 
                           + x₊₊*(F[i]-p₋)/(p₊₊-p₋)/(p₊₊-p₊) + x₊₊/(p₊₊-p₋)*(F[i]-p₊)/(p₊₊-p₊) )
                f[i] = 1 / ∂xI_icdf_i
            end
            =#
            #= 
            if isapprox(x₊, x₋)
                F[i] = p₋
                if i == 1
                    f[i] = F[i] / (xgrid[2] - xgrid[1])
                else
                    f[i] = (F[i] - F[i-1]) / (xgrid[i] - xgrid[i-1])
                end
            else
                F[i] = p₋ + (x - x₋) * (p₊ - p₋) / (x₊ - x₋) 
                f[i] = (p₊ - p₋) / (x₊ - x₋) 
            end
            =#
        end
    end
end

function get_pdf_from_cdf(a_cdf, xgrid; order=2)
    N = length(a_cdf)
    a_pdf_rr = zeros(N)
    if order == 1
        for i in eachindex(a_cdf)
            if i == 1
                a_pdf_rr[i] = (a_cdf[i] - 0) / (xgrid[i+1] - xgrid[i])
            elseif i == N
                a_pdf_rr[i] = (a_cdf[i] - a_cdf[i-1]) / (xgrid[i] - xgrid[i-1])
            else
                a_pdf_rr[i] = (a_cdf[i] - a_cdf[i-1]) / (xgrid[i] - xgrid[i-1])
            end
        end
    elseif order == 2
        for i in eachindex(a_cdf)
            if i == 1
                a_pdf_rr[i] = (-3/2*a_cdf[1] + 2*a_cdf[1+1] - 1/2*a_cdf[1+2]) / (xgrid[i+1] - xgrid[i]) #2nd order forward
            elseif i == N
                a_pdf_rr[i] = (3/2*a_cdf[i-2] - 2*a_cdf[i-1] + 1/2*a_cdf[i]) / (xgrid[i] - xgrid[i-1]) #2nd order backward
            else
                a_pdf_rr[i] = (a_cdf[i] - a_cdf[i-1]) / (xgrid[i+1] - xgrid[i-1])
            end
        end
    else
        error("not implemented")
    end
    return a_pdf_rr
end

"""
cdf_to_pdf!(f, F, xgrid)

Given a discrete cdf F, writes its discrete pdf into f
"""
function cdf_to_pdf!(f, F, Δx)
    N = length(F)
    f[1] = (-3/2*F[1] + 2*F[1+1] - 1/2*F[1+2]) / Δx #2nd order forward
    #f[1] = (F[1] - 0) / Δx #1st order backwards   
    #for i in 1:2
    #    f[i] = (-25/12*F[i] + 4*F[i+1] - 3*F[i+2] + 4/3*F[i+3] - 1/4*F[i+4] ) / Δx #4th order forward
    #end

    for i in 3:(N-2)    
    #for i in 2:(N-1)
        #f[i] = (F[i+1] - F[i]) / Δx #1st order fwd 
        #f[i] = (1/12*F[i-2] - 2/3*F[i-1] + 2/3*F[i+1] - 1/12*F[i+2]) / Δx #4th order centered
        f[i] = 1/2*(F[i+1] - F[i-1]) / Δx #2nd order centered
    end
    #f[N] = (1 - F[N]) / Δx #1st order fwd 
    f[N] = (3/2*F[N] - 2*F[N-1] + 1/2*F[N-2]) / Δx #2nd order backward     
    #for i in (N-1):N
    #   f[i] = (25/12*F[i] - 4*F[i-1] + 3*F[i-2] - 4/3*F[i-3] + 1/4*F[i-4] ) / Δx #4th order backward
    #end
end

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

function add_iicdf(F⁻¹_v, pgrid_v, tol=1e-12)

    n = length(F⁻¹_v)

    M = 0
    for i in 1:n
        M += length(pgrid_v[i])
    end
    
    ΣF⁻¹_v = [ zeros(M) for _ in 1:n]
    Σpgrid = zeros(M)
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

#        if p₀ - p₀₋ < tol && j != 1 # if new p₀ is identical with a previous one, no need to duplicate
#            continue
#        end

        # add p-grid value to the grid union
        #push!(Σpgrid, p₀)
        Σpgrid[j] = p₀
        for i in 1:n
            # calculate interpolant value
            if i_vec[i] == 0 # this icdf is not yet different from zero
                #push!(ΣF⁻¹_v[i], 0.0)
                ΣF⁻¹_v[i][j] = 0
            elseif i_vec[i] == 1 # first icdf value different from zero
                if length(pgrid_v[i]) == 1 # dirac delta
                    #push!(ΣF⁻¹_v[i], F⁻¹_v[i][i_vec[i]])
                    ΣF⁻¹_v[i][j] = F⁻¹_v[i][i_vec[i]]
                else
                    p₋ = pgrid_v[i][i_vec[i]]
                    p₊ = pgrid_v[i][i_vec[i]+1]
                    F⁻¹₋ = F⁻¹_v[i][i_vec[i]]
                    F⁻¹₊ = F⁻¹_v[i][i_vec[i]+1]
                    F⁻¹_i = F⁻¹₋ + (p₀ - p₋) * (F⁻¹₊ - F⁻¹₋) / (p₊ - p₋)
                    #push!(ΣF⁻¹_v[i], F⁻¹_i)
                    ΣF⁻¹_v[i][j] = F⁻¹_i
                end
            else # linear interpolation
                p₋ = pgrid_v[i][i_vec[i]-1]
                p₊ = pgrid_v[i][i_vec[i]]
                F⁻¹₋ = F⁻¹_v[i][i_vec[i]-1]
                F⁻¹₊ = F⁻¹_v[i][i_vec[i]]
                F⁻¹_i = F⁻¹₋ + (p₀ - p₋) * (F⁻¹₊ - F⁻¹₋) / (p₊ - p₋)
                #push!(ΣF⁻¹_v[i], F⁻¹_i)   
                ΣF⁻¹_v[i][j] = F⁻¹_i        
            end
        end
    end
    return ΣF⁻¹_v, Σpgrid
end