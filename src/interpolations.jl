function get_interpolates(params, Λ, mass)
    
    itp_params = ()

    for nd in params
        if length(nd) != 1
            itp_params = (itp_params..., nd)
        end
    end

    n = length(Λ[1])
    lambda_vec_vec = [ [Λ[i][j] for i in eachindex(Λ)] for j in 1:n]

    paramlengths = Tuple(length(i) for i in itp_params)
    mass_r = reshape(mass, paramlengths);

    lambda_mat_vec = [zeros(paramlengths) for _ in 1:n]
    for i in 1:n
        lambda_mat_vec[i] = reshape(lambda_vec_vec[i], paramlengths)
    end

    itp = [interpolate(itp_params, lambda_mat_vec[i], Gridded(Linear())) for i in 1:n];
    mass_itp = interpolate(itp_params, mass_r, Gridded(Linear()));

    return itp, mass_itp
end


function proj_to_simplex(λ)
    u = sort(λ)

    ρ = 1
    for i in eachindex(u)
        if u[i] + 1/i * (1 - sum(u[1:i])) > 0
            ρ = i
        end
    end
    δ = 1/ρ * (1 - sum(u[1:ρ]))

    return [ maximum([0.0, λ[i] + δ]) for i in eachindex(λ) ]
end