function get_S₁(S, Δx)
    S₁ = zero(S)
    N, nₜₚ = size(S)
    mass = zeros(nₜₚ)
    for i in 1:size(S,2)
        mass[i] = sum(S[:,i]) * Δx
        S₁[:,i] .= S[:,i] ./ mass[i]
    end
    return S₁, mass
end

function get_s(S₁, Δx, xgrid, pgrid, tol_cdf=1e-12, tol_icdf=1e-6)
    N, nₜₚ = size(S₁)
    M = length(pgrid)
    s_cdf = zeros(N)
    s = [zeros(M) for _ in 1:nₜₚ]

    for i in eachindex(s)
        MorporJ.cdf!(s_cdf, S₁[:,i], Δx, tol_cdf)
        MorporJ.inv_F!(s[i], s_cdf, xgrid, pgrid, tol_icdf)
    end

    return s
end

function test_interpolations(λ, mass, a, params, massₜ, paramsₜ, Sₜ, sₜ, ciₜ, Δx, xgrid, pgrid, tol_icdf)
    M = length(pgrid)
    nₜₚ = length(λ)
    N, Nₜₚ = size(Sₜ)
    n = length(a)

    # build the interpolations
    itp, mass_itp = MorporJ.get_interpolates(params, λ, mass);

    Sᵣ = zero(Sₜ)
    Sₑ = zero(Sₜ)

    mᵣ = zeros(Nₜₚ)

    _s = zeros(M)
    _s_cdf = zeros(N)
    _s_pdf = zeros(N)

    for k in 1:Nₜₚ
        # recover the parameter values
        indices = Tuple(ciₜ[k])
    
        _params = Tuple( paramsₜ[i][indices[i]] for i in 1:length(paramsₜ) )
        _itp_params = ()
        for i in 1:length(_params)
            if length(paramsₜ[i]) != 1
                _itp_params = (_itp_params..., _params[i])
            end
        end
    
        # interpolate mass and weights
        _m = mass_itp(_itp_params...)
        _λ = [ itp[i](_itp_params...) for i in 1:n  ]
        _λ = MorporJ.proj_to_simplex(_λ)
    
        # barycenter icdf
        _s .= 0
        for i in 1:n
            _s += _λ[i] .* a[i]
        end
    
        # calculate the density of the reconstruction
        MorporJ.inv_F!(_s_cdf, _s, pgrid, xgrid, tol_icdf)
        MorporJ.cdf_to_pdf!(_s_pdf, _s_cdf, Δx, 1)
        Sᵣ[:,k] .= _s_pdf .* _m
    
        # calculate the density of the exact reconstruction
        MorporJ.inv_F!(_s_cdf, sₜ[k], pgrid, xgrid, tol_icdf)
        MorporJ.cdf_to_pdf!(_s_pdf, _s_cdf, Δx, 1)
        Sₑ[:,k] .= _s_pdf .* massₜ[k]
    
        mᵣ[k] = _m
    end

    return Sᵣ, mᵣ, Sₑ
end