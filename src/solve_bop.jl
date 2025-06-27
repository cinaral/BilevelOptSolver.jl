"""
0. Create BilevelOptProb (symbolics in setup_bop)
1. Initialize valid v = [x₁, x₂, λ, s] i.e. solves the follower's problem (initialize_z!())
2. while True:
    2.1. Find all i such that x ∈ Hᵢ (find_feas_index_sets())
    2.2. Construct Hᵢ sets to get BOPᵢ
    2.3. Check if x, λ can solve BOPᵢ for all i: x ∈ Hᵢ (setup_Λ_feas_solver()):
        If True:
            2.3.1. (Optional*) Check if x, λ is a minimum
            2.3.2. Success
        Otherwise: Solve for x and λ for the violating i (setup_leader_nlp()), go to 2.1.

*After we find a feasible Λ, we could check if (x, λ, Λ) is indeed a minimum of (BOPᵢ) by writing the sufficient conditions or using an NLP solver again.

```
Verbosity:
    0: off
    1: minimum: warnings and errors
    2: basic: iterations
    3: extended: number of index sets, infeasible index sets
    4: full: show all index sets
    5: function trace
    6: full: v trace
```
"""
function solve_bop(bop; x_init=zeros(bop.n1 + bop.n2), tol=1e-6, max_iter=100, verbosity=0, n_J_max=20, is_using_HSL=false, seed=0, max_restart_count=10, is_check_v_agreem=false, v_agreement_tol=1e-3, is_using_PATH_to_init=false, conv_tol=1e-3)

    if is_using_HSL && !haskey(ENV, "HSL_PATH")
        is_using_HSL = false
        if verbosity > 0
            print("HSL_PATH not found: Setting is_using_HSL = false! If you would like to use HSL, please obtain a license and download HSL_jll.jl (https://licences.stfc.ac.uk/products/Software/HSL/LibHSL), and set HSL_PATH environment variable to the extracted location.\n")
        end
    end

    if bop.n_θ > 300 && !haskey(ENV, "PATH_LICENSE_STRING")
        if verbosity > 0
            print("PATH_LICENSE_STRING not found and problem size is too large: Please obtain a license (https://pages.cs.wisc.edu/~ferris/path/julia/LICENSE), and set PATH_LICENSE_STRING environment variable.\n")
        end
    end

    nx = bop.n1 + bop.n2

    iter_count = 0
    is_converged = false
    is_sol_valid = false
    is_v_valid = false

    x::Vector{Float64} = zeros(nx)
    v = zeros(bop.nv)
    Λ = zeros(bop.mΛ)

    v_temp = copy(v)
    Λ_temp = copy(Λ)


    v[bop.v_inds["x"]] .= x_init[1:nx]

    Ghs_l = copy(bop.Ghs_l₀)
    Ghs_u = copy(bop.Ghs_u₀)
    θ_l = copy(bop.θ_l₀)
    θ_u = copy(bop.θ_u₀)
    θ_init = zeros(bop.n_θ)

    prev_iter_v = zeros(bop.nv)
    prev_iter_v .= v

    is_Λ_feas_arr = fill(false, n_J_max)
    restart_count = -1

    if isnothing(seed)
        seed = 0
    end
    fol_feas_set_tol = copy(tol)

    #v_arr = [Vector{Float64}(undef, bop.nᵥ) for _ = 1:n_J_max]

    #Λ_all_arr = [Vector{Float64}(undef, bop.m₁ + bop.mₕ + 2 * bop.nᵥ) for _ = 1:n_J_max]
    #is_v_valid_arr = fill(false, n_J_max)


    #θ_init = zeros(bop.n_θ)

    #function terminate()
    #    break
    #end

    while !is_converged
        if iter_count >= max_iter
            if verbosity > 0
                print("Max iterations reached!\n")
            end
            break
        end
        iter_count += 1


        if verbosity > 1
            print("--Iteration $iter_count\n")
        end

        if !is_v_valid
            if verbosity > 2
                print("v isn't valid, initializing v\n")
            end
            init_z_success = initialize_z!(v, bop; is_using_HSL, verbosity, is_using_PATH=is_using_PATH_to_init, tol)

            restart_count += 1
            if restart_count >= max_restart_count
                if verbosity > 0
                    print("Reached maximum restart count, terminating!\n")
                end
                break
            end

            if !init_z_success
                if verbosity > 1
                    print("Failed to initialize z! Randomly initializing and restarting...\n")
                end

                v[bop.v_inds["x"]] .= 10^(restart_count) * (0.5 .- rand(MersenneTwister(seed + restart_count), bop.n1 + bop.n2))
            end
        end

        # By this point, v must at least satisfy the follower's problem
        if verbosity > 5
            print("Computing feasible sets for the follower at v: ")
            display(v)
        end

        if !is_follower_KKT_satisfied(bop, v)
            if verbosity > 1
                print("wtf somehow follower KKT isn't satisfied\n")
            end
            continue
        end
        follow_feas_Js = compute_follow_feas_ind_sets(bop, v; tol=fol_feas_set_tol)

        if length(follow_feas_Js) < 1
            if verbosity > 2
                print("Could not compute follower's feasible sets despite satisfying follower's KKT! Maybe it's a tolerance issue?\n")
            end
            restart_count += 1

            if restart_count > max_restart_count
                if verbosity > 0
                    print("Reached max restarts!\n")
                end
                break
            end
            fol_feas_set_tol = min(10^(restart_count) * fol_feas_set_tol, 1e-2)
            if verbosity > 1
                print("Relaxing follower feasible set tolerance to $(round(fol_feas_set_tol, sigdigits=2)) and restarting...\n")
            end
            continue
        end
        #else
        #    if verbosity > 0
        #        print("Despite claimed initialize_z! success, somehow follower's KKT are not satisfied!\n")
        #    end
        #    break
        #end

        n_J = length(follow_feas_Js)

        if verbosity > 1
            if n_J > 1
                print("Multiple feasible sets ($n_J) detected!\n")
            end
        end

        # we expect a maximum n_J_max number of manifolds
        if n_J > n_J_max
            if verbosity > 0
                print("Too many feasible sets (max: $(n_J_max)): $(n_J)!\n")
            end
            n_J = n_J_max
        end

        is_Λ_feas_arr .= false

        # check if there exists Λ for all i: x ∈ Hᵢ
        # this step could be done in parallel
        for i in 1:n_J
            Ji = follow_feas_Js[i]
            Ji_bounds = convert_J_to_bounds(Ji, bop)
            # does there exist Λ for v? feasibility problem so is_minimizing=false
            is_Λ_feas = update_v!(v_temp, Λ_temp, Ghs_l, Ghs_u, θ_l, θ_u, θ_init, bop, Ji_bounds; is_minimizing=false, is_using_HSL, tol)

            # debug
            is_Λ_feas_2 = bop.check_Λ_lp_feas(v, Ji_bounds.z_l, Ji_bounds.z_u)

            if (is_Λ_feas && !is_Λ_feas_2) || (!is_Λ_feas && is_Λ_feas_2)
                #@info is_Λ_feas
                #@info is_Λ_feas_2
                #Main.@infiltrate
            end

            if is_Λ_feas
                is_Λ_feas_arr[i] = true
                if verbosity > 2
                    print("i=$i: Feasible J1/2/3/4: $(Ji[1])/$(Ji[2])/$(Ji[3])/$(Ji[4])\n")
                end
            else
                if verbosity > 1
                    print("i=$i: INFEASIBLE J1/2/3/4: $(Ji[1])/$(Ji[2])/$(Ji[3])/$(Ji[4])\n")
                end
            end
        end

        is_all_Λ_feas = all(is_Λ_feas_arr[1:n_J])

        if is_all_Λ_feas
            if verbosity > 4
                print("All $n_J follower solutions have feasible Λ.\n")
            end
        else
            continue
        end

        # by this point all follower solutions have feasible v, Λ
        is_v_valid = false
        is_v_sol = false
        #is_v_valid = true

        # We must check all BOPᵢ solutions to ensure the solution is valid, but first we need a solution
        v_temp .= v
        Λ_temp .= Λ

        for i in 1:n_J
            Ji = follow_feas_Js[i]
            Ji_bounds = convert_J_to_bounds(Ji, bop)

            if !is_v_valid || is_check_v_agreem
                is_v_sol = update_v!(v_temp, Λ_temp, Ghs_l, Ghs_u, θ_l, θ_u, θ_init, bop, Ji_bounds; is_minimizing=true, is_using_HSL, tol)

                if !is_v_sol
                    if verbosity > 4
                        print("Could not solve BOPᵢ i=$i\n")
                    end
                    continue
                end
                is_v_valid = true
                v .= v_temp
                Λ .= Λ_temp
                # we just need one solution
                break
            end
        end

        if !is_v_valid
            if verbosity > 1
                print("v is not valid (could not solve any BOPᵢ)!\n")
            end
            continue
        end

        # now we check if the solution is valid. for this we need to find if the solution v is valid on all solutions, or if all solutions agree
        is_sol_valid = true

        for i in 1:n_J
            Ji = follow_feas_Js[i]
            Ji_bounds = convert_J_to_bounds(Ji, bop)

            if !is_check_v_agreem
                is_v_Λ_valid_KKT = check_v_Λ_BOPᵢ_KKT(v, Λ, bop, Ghs_l, Ghs_u; verbosity, tol)

                if !is_v_Λ_valid_KKT
                    if verbosity > 4
                        print("v and Λ did not satisfy KKT of BOPᵢ i=$i\n")
                    end
                    is_sol_valid = false
                    break
                end
                continue
            end

            is_v_temp_sol = update_v!(v_temp, Λ_temp, Ghs_l, Ghs_u, θ_l, θ_u, θ_init, bop, Ji_bounds; is_minimizing=true, is_using_HSL, tol)

            if !is_v_temp_sol
                if verbosity > 1
                    print("Could not solve BOPᵢ i=$i so we cannot check for v agreement!\n")
                end
                is_sol_valid = false
                break
            end

            norm_x_err = LinearAlgebra.norm(v_temp[bop.v_inds["x"]] - v[bop.v_inds["x"]]) # only checking x error
            if (norm_x_err > v_agreement_tol)
                if verbosity > 1
                    print("BOPᵢ solutions disagree! norm x err: $norm_x_err\n")
                end
                is_sol_valid = false
                break
            end
        end
        # by this point is_sol_valid=true: v's are sols and agree, Λs are feasible, all we have to check now is if the solution has converged
        dv = v[bop.v_inds["x"]] - prev_iter_v[bop.v_inds["x"]]
        prev_iter_v .= v
        #@info LinearAlgebra.norm(dv)
        if LinearAlgebra.norm(dv) < conv_tol
            is_converged = true
            break
        else
            restart_count = 0
            if verbosity > 3
                print("Found new v that's a valid solution ")
            end
            if verbosity > 5
                display(v)
            elseif verbosity > 3
                print("\n")
            end
        end
    end

    if is_converged
        if verbosity > 1
            print("Converged in $iter_count iterations\n")
        end
    end

    x .= v[bop.v_inds["x"]]

    # final sanity check
    if is_sol_valid && !(all(bop.G(x) .≥ 0 - tol) && all(bop.g(x) .≥ 0 - tol))# && is_follower_KKT_satisfied(bop, v))
        if verbosity > 0
            print("Something went wrong, this allegedly valid solution was primal infeasible!\n")
        end
        is_sol_valid = false
    end

    (; x, is_converged, is_sol_valid, iter_count)
end

function check_v_Λ_BOPᵢ_KKT(v, Λ, bop, Ghs_l, Ghs_u; verbosity=0, tol=1e-6)
    ∇ᵥGhs_rows, ∇ᵥGhs_cols, ∇ᵥGhs_shape, ∇ᵥᵥL_rows, ∇ᵥᵥL_cols, ∇ᵥᵥL_shape = bop.info_BOPᵢ()
    Ghs = zeros(bop.mΛ)
    ∇ᵥF = zeros(bop.nv)
    ∇ᵥGhs = sparse(∇ᵥGhs_rows, ∇ᵥGhs_cols, zeros(length(∇ᵥGhs_rows)), ∇ᵥGhs_shape[1], ∇ᵥGhs_shape[2])
    ∇ᵥᵥL = sparse(∇ᵥᵥL_rows, ∇ᵥᵥL_cols, zeros(length(∇ᵥᵥL_rows)), ∇ᵥᵥL_shape[1], ∇ᵥᵥL_shape[2])
    bop.eval_BOPᵢ!(Ghs, ∇ᵥF, ∇ᵥGhs.nzval, ∇ᵥᵥL.nzval, v, Λ)

    is_stationary = all(isapprox.(∇ᵥF - ∇ᵥGhs' * Λ, 0; atol=2 * tol))
    is_complement = all(isapprox.(Λ .* Ghs, 0; atol=2 * tol)) # complementarity
    is_primal_feas = all(Ghs .≥ Ghs_l .- tol) && all(Ghs .≤ Ghs_u .+ tol)

    is_sol_valid = true
    if !(is_stationary && is_complement && is_primal_feas)
        is_sol_valid = false
    end
    is_sol_valid
end

function is_follower_KKT_satisfied(bop, v; tol=1e-6)
    x = @view v[bop.v_inds["x"]]
    λ = @view v[bop.v_inds["λ"]]

    ∇ₓf = zeros(bop.n1 + bop.n2)
    bop.deriv_funs.∇ₓf!(∇ₓf, x)
    ∇ₓg_vals = zeros(length(bop.deriv_funs.∇ₓg_rows))
    bop.deriv_funs.∇ₓg_vals!(∇ₓg_vals, x)
    ∇ₓg = sparse(bop.deriv_funs.∇ₓg_rows, bop.deriv_funs.∇ₓg_cols, ∇ₓg_vals)
    x2_inds = bop.n1+1:bop.n1+bop.n2
    if bop.m2 > 1
        is_stationary = all(isapprox.(∇ₓf[x2_inds] - ∇ₓg[:, x2_inds]' * λ, 0; atol=2 * tol))
    else
        is_stationary = all(isapprox.(∇ₓf[x2_inds], 0; atol=2 * tol))
    end
    is_primal_feas = all(bop.g(x) .≥ 0 - tol) # primal feas
    is_dual_feas = all(λ .≥ 0 - tol) # dual feas
    is_complement = all(isapprox.(λ .* bop.g(x), 0; atol=2 * tol)) # complementarity

    return is_stationary && is_primal_feas && is_dual_feas && is_complement
end

function initialize_z!(v, bop; verbosity=0, is_using_PATH=false, is_using_HSL=false, tol=1e-6)
    # if BOPᵢ wasn't solved the low level solution may be invalid, and we have to call the follower nlp
    x1 = @view v[bop.v_inds["x1"]]
    x2 = @view v[bop.v_inds["x2"]]
    λ = zeros(bop.m2)

    if is_using_PATH
        θ_out, status, _ = bop.solve_follower_KKT_mcp(x1; tol)

        x2 = θ_out[1:bop.n2]
        init_z_success = status == PATHSolver.MCP_Solved
    else
        x2, λ, solvestat, _ = bop.solve_follower_nlp(x1; x2_init=x2, is_using_HSL, tol)
        init_z_success = solvestat == 0 || solvestat == 1 # accept Solve_Succeeded and Solved_To_Acceptable_Level
    end

    # if the feasible region of the follower is empty for x₁, this finds a bilevel feasible x to move to
    if !init_z_success
        if verbosity > 1
            print("Resetting x to a bilevel feasible point. \n")
        end

        x_feas, _, solvestat, _ = bop.find_bilevel_feas_pt(x_init=[x1; x2]; tol)
        bilevel_feas_success = solvestat == 0 || solvestat == 1

        if bilevel_feas_success
            x1 .= x_feas[bop.v_inds["x1"]]
            x2 .= x_feas[bop.v_inds["x2"]]

            if is_using_PATH
                θ_out, status, _ = bop.solve_follower_KKT_mcp(x1; tol)
                x2 = θ_out[1:bop.n2]
                init_z_success = status == PATHSolver.MCP_Solved
            else
                x2, λ, solvestat, _ = bop.solve_follower_nlp(x1; x2_init=x2, is_using_HSL, tol)
                init_z_success = solvestat == 0 || solvestat == 1 # accept Solve_Succeeded and Solved_To_Acceptable_Level
            end
        else
            if verbosity > 0
                print("Failed to find a bilevel feasible point to re-attempt, the problem may be bilevel infeasible! Check your x_init, G(x), and g(x).\n")
            end
        end
    end

    if init_z_success
        v[bop.v_inds["x1"]] .= x1
        v[bop.v_inds["x2"]] .= x2
        v[bop.v_inds["λ"]] .= λ
        v[bop.v_inds["s"]] .= bop.g([x1; x2])
    end

    return init_z_success
end

"""
is_minimizing = false: BOPᵢ KKT is solved as an MCP
is_minimizing = true: BOPᵢ is solved as an NLP
is_using_HSL = true: Use HSL LP solver back-end when solving the NLP
"""
function update_v!(v, Λ, Ghs_l, Ghs_u, θ_l, θ_u, θ_init, bop, Ji_bounds; is_minimizing=false, is_using_HSL=false, tol=1e-6)
    is_v_valid::Bool = false
    if is_minimizing
        Ghs_l[bop.Ghs_inds["h"]] .= Ji_bounds.h_l
        Ghs_u[bop.Ghs_inds["h"]] .= Ji_bounds.h_u
        Ghs_l[bop.Ghs_inds["s"]] .= Ji_bounds.z_l[bop.z_inds["s"]]
        Ghs_u[bop.Ghs_inds["s"]] .= Ji_bounds.z_u[bop.z_inds["s"]]
        
        v_out, Λ_out, solvestat = bop.solve_BOPᵢ_nlp(g_l=Ghs_l, g_u=Ghs_u, x_init=v; is_using_HSL, tol)
        is_v_valid = solvestat == 0 || solvestat == 1
    else
        θ_init[bop.θ_inds["v"]] .= v
        θ_l[bop.θ_inds["z"]] .= Ji_bounds.z_l
        θ_u[bop.θ_inds["z"]] .= Ji_bounds.z_u
        θ_l[bop.θ_inds["rh"]] .= Ji_bounds.h_l
        θ_u[bop.θ_inds["rh"]] .= Ji_bounds.h_u

        θ_out, status, _ = bop.solve_BOPᵢ_KKT_mcp(x_l=θ_l, x_u=θ_u, x_init=θ_init; tol)
        v_out = θ_out[bop.θ_inds["v"]]
        Λ_out = θ_out[bop.θ_inds["Λ"]]
        is_v_valid = status == PATHSolver.MCP_Solved
    end

    v .= v_out # even if it's not solved, we update v so we can try to initialize z again
    Λ .= Λ_out
    is_v_valid
end

"""
Check viable ways to satisfy h ⟂ z (aka Λₕ)

This function constructs, j∈{1,…,mₕ}:
        { j: hⱼ > 0, l = z    } case 1  : hⱼ inactive (positive), zⱼ at lb
        { j: hⱼ = 0, l = z    } case 1/2: ambiguous
K[j] =  { j: hⱼ = 0, l < z < u} case 2  : hⱼ active (l ≠ u)
        { j: hⱼ = 0,     z = u} case 2/3: ambiguous 
        { j: hⱼ < 0,     z = u} case 3  : hⱼ inactive (negative), zⱼ at ub
        { j: hⱼ,     l = z = u} case 4  : hⱼ free, zⱼ fixed

Then we sort out ambiguities by enumerating and collect them in J[1], J[2], J[3] and J[4]
"""
function check_is_sol_valid(bop, v; tol=1e-3)
    Gh = zeros(bop.m1 + bop.mh)
    bop.Gh!(Gh, v)
    h = @view Gh[bop.Ghs_inds["h"]]
    z = @view v[bop.v_inds["z"]]
    z_l = bop.z_l₀
    z_u = bop.z_u₀
    K = Dict{Int,Vector{Int}}()

    # note which constraints are active
    for j in 1:bop.mh
        Kj = Int[]

        if isapprox(z_l[j], z_u[j]; atol=2 * tol)
            push!(Kj, 4) # case 4
        elseif tol ≤ h[j] && z[j] < z_l[j] + tol
            push!(Kj, 1) # case 1
        elseif -tol < h[j] < tol && z[j] < z_l[j] + tol
            push!(Kj, 1) # case 1 OR  
            push!(Kj, 2) # case 2
        elseif -tol < h[j] < tol && z_l[j] + tol ≤ z[j] ≤ z_u[j] - tol
            push!(Kj, 2) # case 2
        elseif -tol < h[j] < tol && z_u[j] - tol < z[j]
            push!(Kj, 2) # case 2 OR 
            push!(Kj, 3) # case 3
        elseif h[j] ≤ tol && z_u[j] - tol < z[j]
            push!(Kj, 2) # case 3
        end
        K[j] = Kj
    end
    is_sol_valid = !any(isempty.(Kj for Kj in values(K)))
    is_sol_valid, K
end

function compute_follow_feas_ind_sets(bop, v; tol=1e-3, verbosity=0)
    Js = Vector{Dict{Int,Set{Int}}}()
    is_sol_valid, K = check_is_sol_valid(bop, v; tol)

    if !is_sol_valid
        if verbosity > 0
            print("Not a valid solution!\n")
        end
        return Js
    end

    # enumerate the ambiguous indexes
    ambigu_inds = [i for i in 1:bop.mh if length(K[i]) > 1]
    single_inds = setdiff(1:bop.mh, ambigu_inds)
    It = Iterators.product([K[i] for i in ambigu_inds]...)

    for assignment in It
        J = Dict(k => Set{Int}() for k in 1:4) # number of cases = 4
        for (j, ej) in enumerate(assignment)
            push!(J[ej], ambigu_inds[j])
        end
        for j in single_inds
            push!(J[K[j][1]], j)
        end
        push!(Js, J)
    end
    Js
end

function convert_J_to_bounds(J, bop)
    h_l = zeros(bop.mh)
    h_u = zeros(bop.mh)
    z_l = zeros(bop.mh)
    z_u = zeros(bop.mh)
    z_l₀ = bop.z_l₀
    z_u₀ = bop.z_u₀

    for j in J[1] # hⱼ inactive (positive), zⱼ at lb
        h_l[j] = 0.0
        h_u[j] = Inf
        z_l[j] = z_l₀[j]
        z_u[j] = z_l₀[j]
    end
    for j in J[2] # hⱼ active (l ≠ u)
        h_l[j] = 0.0
        h_u[j] = 0.0
        z_l[j] = z_l₀[j]
        z_u[j] = z_u₀[j]
    end
    for j in J[3] # hⱼ inactive (negative), zⱼ at ub
        h_l[j] = -Inf
        h_u[j] = 0.0
        z_l[j] = z_u₀[j]
        z_u[j] = z_u₀[j]
    end
    for j in J[4] # hⱼ free, zⱼ fixed (l = u)
        h_l[j] = -Inf
        h_u[j] = Inf
        z_l[j] = z_l₀[j]
        z_u[j] = z_u₀[j]
    end

    (; h_l, h_u, z_l, z_u)
end
