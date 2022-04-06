"""
icdf!(F⁻¹, F, xgrid, pgrid)

Given a discrete cdf F, writes its discrete icdf into F⁻¹

F⁻¹ : [0,1] ∋ p ↦ inf { x: F(x) ≥ p } ∈ [0,L]

Limited to length(F) == length(F⁻¹), uniform grids
"""
function icdf!(F⁻¹, F, xgrid, pgrid, tol = 1e-6)
    N = length(F)
    M = length(F⁻¹)
    i = 1
    for j in 1:M
        p = pgrid[j]
        while F[i] < p && i < N
            i += 1
        end
        if i == 1 # p₀ := F[1] ≥ p[j] ⇒ cdf does not start at zero or p₀ = 0 = F[1]
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
            if i == N || (1 - F[i+1]) < tol # Do not interpolate over the kink.
                if i < 3 # dirac delta case
                    F⁻¹[j] = xgrid[1]
                    continue
                end
                p₋ = F[i-2]
                p₊ = F[i-1]
                p₊₊ = F[i]
                x₋ = xgrid[i-2]
                x₊ = xgrid[i-1]
                x₊₊ = xgrid[i]
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
            #=
            F⁻¹[j] = ( x₋ * (p-p₊)/(p₋-p₊) * (p-p₊₊)/(p₋-p₊₊) 
                    + x₊ * (p-p₋)/(p₊-p₋) * (p-p₊₊)/(p₊-p₊₊) 
                    + x₊₊ * (p-p₋)/(p₊₊-p₋) * (p-p₊)/(p₊₊-p₊) )
            =#

            F⁻¹[j] = x₋ + (p-p₋)/(p₊-p₋) * (x₊-x₋)

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

Given a discrete icdf F⁻¹, writes its discrete cdf into F

F⁻¹ : [0,1] ∋ p ↦ inf { x: F(x) ≥ p } ∈ [0,L]

Limited to length(F) == length(F⁻¹), uniform grids
"""
function iicdf!(F, F⁻¹, xgrid, pgrid, tol=1e-6)
    N = length(F)
    M = length(F⁻¹)
    j = 1
    for i in 1:N
        x = xgrid[i]
        while F⁻¹[j] < x && j < M
            j += 1
        end
        # now either F⁻¹[j] >= x[i] or there exists no j: F⁻¹[M] <= x[i]
        if j == 1 # x₀ := F⁻¹[1] = F⁻¹(0) >= x[i] ⇒ icdf does not start at x[1] or 
            # return smallest p
            F[i] = pgrid[1]
            continue
        elseif j == M && F⁻¹[j] < x
            # F⁻¹[end] = F⁻¹(1) <= x[i] ⇒ icdf does not reach x[end] ⇒ fill up with ones
            F[i] = pgrid[end]
            continue
        else
            #quadratic interpolation
            #if i == 1 && j <= (M-2) # kink at x = 0 in the cdf
            #    p₋ = pgrid[j]
            #    p₊ = pgrid[j+1]
            #    p₊₊ = pgrid[j+2]
            #    x₋ = F⁻¹[j]
            #    x₊ = F⁻¹[j+1]
            #    x₊₊ = F⁻¹[j+2]
            if j == M
                p₋ = pgrid[j-2]
                p₊ = pgrid[j-1]
                p₊₊ = pgrid[j]
                x₋ = F⁻¹[j-2]
                x₊ = F⁻¹[j-1]
                x₊₊ = F⁻¹[j]
            else
                p₋ = pgrid[j-1]
                p₊ = pgrid[j]
                p₊₊ = pgrid[j+1]
                x₋ = F⁻¹[j-1]
                x₊ = F⁻¹[j]
                x₊₊ = F⁻¹[j+1]
            end

            if abs(x₊-x₋) < tol
                F[i] = p₊
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
                    F[i] = p₊
                    continue
                end
            end

            # interpolant of x ↦ F(x)
            #=
            F[i] = ( p₋ * (x-x₊)/(x₋-x₊) * (x-x₊₊)/(x₋-x₊₊) 
                    + p₊ * (x-x₋)/(x₊-x₋) * (x-x₊₊)/(x₊-x₊₊) 
                    + p₊₊ * (x-x₋)/(x₊₊-x₋) * (x-x₊)/(x₊₊-x₊) )
            =#

            F[i] = p₋ + (x-x₋)/(x₊-x₋) * (p₊-p₋) 

            # fix values outside of allowed domain
            if F[i] < pgrid[1]
                F[i] = pgrid[1]
            elseif F[i] > pgrid[end]
                F[i] = pgrid[end]
            end
        end
    end
end