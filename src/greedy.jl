function greedy_algo!(s, a, Δx, basis_index; tol=1e-6, rtol=1e-1, max_iter_greedy=25, max_iter_solver=1000000)
    
    n = length(a) # number of atoms
    k = length(s) # size of training data
    Λ = [ [zeros(n) for _ in 1:k] ]
    ΔW = [ zeros(k) ]
    ΔW_max = []
    
    Λ_int = [zeros(n) for _ in 1:k]
    ΔW_int = zeros(k)
    
    iter=1
    
    while iter <= max_iter_greedy
        Λ_int, ΔW_int = barycenter_fit(s, a, Δx, max_iter=max_iter_solver)
        
        Δ = maximum(ΔW_int)
        push!(ΔW_max, Δ)
        push!(Λ, Λ_int)
        push!(ΔW, ΔW_int)

        print("Current W2 error: $Δ with atoms:")
        for idx in basis_index
            print(" $idx")
        end
        print(".\n")

        
        # check relative tolerance
        if iter !=1 && (ΔW_max[end-1] - Δ)/ΔW_max[end-1] < rtol
            relative_error = ( ΔW_max[end-1] - Δ)/ΔW_max[end-1]
            print("Relative tolerance not met - $relative_error. Removing last dictionary atom. \n")
            pop!(basis_index) # remove last basis addition that didnt substantially improve the fit
            pop!(a)
            pop!(Λ)
            pop!(ΔW)
            pop!(ΔW_max)
            break
        # check absolute tolerance
        elseif Δ < tol
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
    
    return Λ, ΔW, ΔW_max, a, basis_index
end