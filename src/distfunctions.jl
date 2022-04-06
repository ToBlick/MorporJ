"""
cdf!(f, F)

Given a discrete pdf f, writes its cdf into F. The value in F[i] includes the integral of the ith cell.

Limited to length(f) == length(xgrid)-1 == length(F)
"""
function cdf!(F, f, Δx, tol = 1e-9)
    F[1] = f[1] * Δx
    for i in 2:length(f)
        F[i] = F[i-1] + f[i] * Δx
        if 1 - F[i] < tol
            F[i] = 1
        end
    end
end

"""
iicdf!(F⁻¹, F, xgrid, pgrid)

Given a discrete F⁻¹, writes its discrete inverse into F⁻¹

F⁻¹ : [a,b] ∋ p ↦ inf { x: F(x) ≥ p } ∈ [c,d]

Limited to length(F) == length(F⁻¹), uniform grids
"""
function inv_F!(F⁻¹, F, xgrid, pgrid, tol = 1e-12)
    N = length(F)
    M = length(F⁻¹)
    i = 1
    for j in 1:M
        p = pgrid[j]
        while F[i] < p && i < N
            i += 1
        end
        if i == 1 
            # p₀ := F[1] ≥ p[j] ⇒ F does not start at zero or p₀ = 0 = F[1]
            # return smallest x
            F⁻¹[j] = xgrid[1]
            continue
        # F[i] ≥ p[j] from here on
        elseif  F[i] < p && i == N
            # F[end] < p[j] ⇒ F does not reach 1
            # return largest x
            F⁻¹[j] = xgrid[end]
            continue
        elseif i == N || (1 - F[i+1]) < tol # Do not interpolate over the kink.
            if i < 3 # dirac delta case
                F⁻¹[j] = xgrid[1]
                continue
            end
            p₋ = F[i-2]
            p₊ = F[i-1]
            x₋ = xgrid[i-2]
            x₊ = xgrid[i-1]
        else
            #linear interpolation
            p₋ = F[i-1]
            p₊ = F[i]
            x₋ = xgrid[i-1]
            x₊ = xgrid[i]
        end

        # catch singular case
        if abs(p₊-p₋) < tol
            F⁻¹[j] = x₋
            continue
        end

        F⁻¹[j] = x₋ + (p-p₋)/(p₊-p₋) * (x₊-x₋)

        if F⁻¹[j] < xgrid[1]
            F⁻¹[j] = xgrid[1]
        elseif F⁻¹[j] > xgrid[end]
            F⁻¹[j] = xgrid[end]
        end
    end
end

"""
cdf_to_pdf!(f, F, Δx)

Given a discrete cdf F, writes its discrete pdf into f
"""
function cdf_to_pdf!(f, F, Δx, order=1)
    N = length(F)
    if order == 1
        f[1] = (F[1] - 0) / Δx #1st order backwards 
        for i in 2:(N)
            f[i] = (F[i] - F[i-1]) / Δx #1st order backwd. 
        end
    elseif order == 2
        f[1] = (-3/2*F[1] + 2*F[1+1] - 1/2*F[1+2]) / Δx
        #f[1] = (-F[1] + F[1+1]) / Δx
        for i in 2:(N-1)    
            f[i] = 1/2*(F[i+1] - F[i-1]) / Δx #2nd order centered
        end
        f[N] = (3/2*F[N] - 2*F[N-1] + 1/2*F[N-2]) / Δx
        #f[N] = (F[N] - F[N-1]) / Δx
    elseif order == 4
        for i in 1:2
            f[i] = (-25/12*F[i] + 4*F[i+1] - 3*F[i+2] + 4/3*F[i+3] - 1/4*F[i+4] ) / Δx #4th order forward
        end
        for i in 3:(N-2)    
            f[i] = (1/12*F[i-2] - 2/3*F[i-1] + 2/3*F[i+1] - 1/12*F[i+2]) / Δx #4th order centered
        end
        for i in (N-1):N
            f[i] = (25/12*F[i] - 4*F[i-1] + 3*F[i-2] - 4/3*F[i-3] + 1/4*F[i-4] ) / Δx #4th order backward
        end
    else
        error("implemented orders: 1, 2, 4")
    end
end
