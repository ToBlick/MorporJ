function greedy_algo!(s, a, Δx, basis_index; tol=1e-6, rtol=1e-1, max_iter_greedy=25, max_iter_solver=1000000)
    
    n = length(a) # number of atoms
    k = length(s) # size of training data
    Λ = [ ]
    ΔW = [ ]
    ΔW_max = []
    ΔW_avg = []
    
    Λ_int = [zeros(n) for _ in 1:k]
    ΔW_int = zeros(k)
    
    iter=1
    
    while iter <= max_iter_greedy
        Λ_int, ΔW_int = barycenter_fit(s, a, Δx, max_iter=max_iter_solver)
        
        Δ = maximum(ΔW_int)
        push!(ΔW_max, Δ)
        push!(Λ, Λ_int)
        push!(ΔW, ΔW_int)
        push!(ΔW_avg, sum(ΔW_int)/k )

        @printf "------------------------------------- \n"
        @printf "Iteration %.0f \n" iter
        @printf "Maximum W2 error: %.2e \n" Δ
        @printf "Average W2 error: %.2e \n" ΔW_avg[end]
        @printf "Atoms used:"
        for idx in basis_index
            @printf " %.0f" idx
        end
        @printf " \n"
        @printf "------------------------------------- \n"

        
        # check relative tolerance
        if iter !=1 && (ΔW_avg[end-1] - Δ)/ΔW_avg[end-1] < rtol
            relative_error = ( ΔW_avg[end-1] - Δ)/ΔW_avg[end-1]
            print("Relative tolerance not met: Δᵣₑ = $relative_error. Removing last dictionary atom. \n")
            pop!(basis_index) # remove last basis addition that didnt substantially improve the fit
            pop!(a)
            pop!(Λ)
            pop!(ΔW)
            pop!(ΔW_max)
            pop!(ΔW_avg)
            break
        # check absolute tolerance
        elseif ΔW_avg[end] < tol
            print("Absolute tolerance met.")
            break
        # check maximum iterations
        elseif iter == max_iter_greedy
            print("Maximum iterations reached: $iter \n")
            break
        # add next atom
        else
            i_star=argmax(ΔW_int)
            print("Adding snapshot number: $i_star \n")
            push!(a, s[i_star])
            push!(basis_index, i_star)
        end
        iter+=1
    end
    
    return Λ, ΔW, ΔW_max, ΔW_avg, a, basis_index
end


function get_initial_atoms(s)
    nₜₚ = length(s)

    D = zeros(nₜₚ,nₜₚ)
    for i in 1:nₜₚ
        for j in i:nₜₚ
            D[i,j] = norm(s[i]-s[j],2)
        end
    end
    maxD = Tuple(argmax(D))

    return [s[maxD[2]], s[maxD[1]]], [maxD[2], maxD[1]]
end