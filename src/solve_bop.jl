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
function solve_bop(bop; x_init=zeros(bop.n₁ + bop.n₂), tol=1e-6, max_iter=200, verbosity=0, n_J_max=20, is_using_PATH=false, is_using_HSL=false)
    #x_init::Vector{Float64}=zeros(bop.n₁ + bop.n₂); tol::Float64=1e-3; max_iter::Int64=200; verbosity::Int64=0; n_J_max::Int64=20; is_using_PATH::Bool=false; is_using_HSL::Bool=true

    if is_using_HSL && !haskey(ENV, "HSL_PATH")
        is_using_HSL = false
        if verbosity > 0
            print("HSL_PATH not found: Setting is_using_HSL = false! If you would like to use HSL, please obtain a license and download HSL_jll.jl (https://licences.stfc.ac.uk/products/Software/HSL/LibHSL), and set HSL_PATH environment variable to the extracted location.\n")
        end
    end

    if is_using_PATH && bop.n_θ > 300 && !haskey(ENV, "PATH_LICENSE_STRING")
        is_using_PATH = false
        if verbosity > 0
            print("PATH_LICENSE_STRING not found and problem size is too large: Setting is_using_PATH = false! If you would like to use the PATH Solver for this problem, please obtain a license (https://pages.cs.wisc.edu/~ferris/path/julia/LICENSE), and set PATH_LICENSE_STRING environment variable.\n")
        end
    end

    nₓ = bop.n₁ + bop.n₂

    iter_count = 0
    is_converged = false
    is_sol_valid = false
    is_success = false
    rolling_v_idx = 1

    x::Vector{Float64} = zeros(nₓ)
    #λ = zeros(bop.m₂)
    #s = zeros(bop.m₂)
    v = zeros(bop.nᵥ)

    v[bop.v_inds["x"]] .= x_init[1:nₓ]
    v_l = copy(bop.v_l₀)
    v_u = copy(bop.v_u₀)
    Gh_l = copy(bop.Gh_l₀)
    Gh_u = copy(bop.Gh_u₀)
    θ_l = copy(bop.θ_l₀)
    θ_u = copy(bop.θ_u₀)
    θ_init = zeros(bop.n_θ)

    prev_iter_v = zeros(bop.nᵥ)
    prev_iter_v .= v

    is_BOPᵢ_solved = false

    v_arr = [Vector{Float64}(undef, bop.nᵥ) for _ = 1:n_J_max]
    is_Λ_feas_arr = fill(false, n_J_max)
    is_BOPᵢ_solved_arr = fill(false, n_J_max)

    #θ_init = zeros(bop.n_θ)

    #function terminate()
    #    break
    #end

    while !is_converged
        iter_count += 1

        if iter_count > max_iter
            if verbosity > 0
                print("Max iterations reached!\n")
            end
            break
        end

        if verbosity > 1
            print("--Iteration $iter_count\n")
        end

        if !is_BOPᵢ_solved
            init_z_success = initialize_z!(v, bop; is_using_HSL, verbosity, is_using_PATH)

            if !init_z_success
                if verbosity > 0
                    print("Failed to initialize z!\n")
                end
                break
            end
        end

        # By this point, v must satisfy the follower's problem
        if verbosity > 5
            print("Computing feasible sets for the follower at v: ")
            display(v)
        end

        follow_feas_Js = compute_follow_feas_ind_sets(bop, v)
        if length(follow_feas_Js) < 1
            if verbosity > 0
                # this should never happen, somehow follower KKT isn't satisfied:
                @info "is follower KKT satisfied? $(is_follower_KKT_satisfied(bop, v))"
                print("Invalid v: Could not compute follower feasible sets!\n")
            end
        end

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

        for i in 1:n_J
            v_arr[i] .= zeros(bop.nᵥ)
        end
        is_Λ_feas_arr[1:n_J] .= false
        is_BOPᵢ_solved_arr[1:n_J] .= false
        #is_Js_outdated = false

        for i in 1:n_J
            Ji = follow_feas_Js[i]
            Ji_bounds = convert_J_to_bounds(Ji, bop)

            is_Λ_feas = bop.check_Λ_lp_feas(v, Ji_bounds.z_l, Ji_bounds.z_u, Ji_bounds.h_l, Ji_bounds.h_u) # does there exist Λ_all=[Λ; Λ_v_l; Λ_v_u] for v

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

            # if we can't confirm there exists a feasible Λ at this J, we need to update v:
            if !is_Λ_feas
                is_BOPᵢ_solved = update_v!(v, bop, Ji_bounds, v_l, v_u, Gh_l, Gh_u, θ_l, θ_u, θ_init; is_using_PATH, is_using_HSL)

                if is_BOPᵢ_solved
                    if verbosity > 2
                        print("i=$i: Λ not feasible, but we could solve BOPᵢ.\n")
                    end
                    is_BOPᵢ_solved_arr[i] = true
                    v_arr[i] .= v
                    #is_Js_outdated = true
                    break
                end
            end
        end

        # By this point, all manifolds must be feasible
        if all(is_Λ_feas_arr[1:n_J])
            if verbosity > 4
                print("There exists Λ for all ($n_J) follower feasible index sets.\n")
            end

            # We must compute the BOPᵢ solutions to ensure the solution is valid
            # TODO: or... check some other way
            for i in 1:n_J
                Ji = follow_feas_Js[i]
                Ji_bounds = convert_J_to_bounds(Ji, bop)
                if !is_BOPᵢ_solved_arr[i]
                    is_BOPᵢ_solved = update_v!(v, bop, Ji_bounds, v_l, v_u, Gh_l, Gh_u, θ_l, θ_u, θ_init; is_using_PATH, is_using_HSL)

                    if is_BOPᵢ_solved
                        if verbosity > 2
                            print("i=$i: BOPᵢ solved.\n")
                        end
                        is_BOPᵢ_solved_arr[i] = true
                        v_arr[i] .= v
                    end
                end
            end
        end

        if !all(is_BOPᵢ_solved_arr[1:n_J])
            if !any(is_BOPᵢ_solved_arr[1:n_J])
                # oh no, we have no valid v, maybe try again?
                if verbosity > 1
                    is_sol_valid = false
                    print("None BOPᵢ could be solved (out of $(n_J)). v is not valid!\n")
                end
                continue
            end

            for (i, _) in enumerate(is_BOPᵢ_solved_arr[1:n_J])
                # we have some potentially valid solutions to choose from
                if is_BOPᵢ_solved_arr[i]
                    v = v_arr[i]
                    prev_iter_v .= v_arr[i]
                    if verbosity > 3
                        print("v updated: ")
                    end
                    if verbosity > 5
                        display(v)
                    elseif verbosity > 3
                        print("\n")
                    end
                    break
                end
            end
        end

        # by this point all BOPᵢ must be solved, but we don't know if their solutions agree
        is_sol_valid = true

        if verbosity > 1
            print("All ($n_J) BOPᵢ are solved.\n")
        end

        if n_J > 1
            for i in 2:n_J
                norm_x_err = LinearAlgebra.norm(v_arr[i][1:nₓ] - v_arr[1][1:nₓ]) # only checking x error
                if (norm_x_err > tol)
                    if verbosity > 1
                        print("BOPᵢ solutions disagree! norm x err: $norm_x_err\n")
                    end
                    is_sol_valid = false
                    break
                end
            end
        end

        # TODO: selecting new v 
        # we choose a new v arbitrarily, if there are multiple choices, we try to select another index every iteration to avoid getting stuck, but it also breaks some solutions. we should check if the solution is valid here
        check_count = 0
        is_v_valid, _ = check_is_sol_valid(bop, v_arr[rolling_v_idx])
        while !is_v_valid
            rolling_v_idx += 1
            if rolling_v_idx > n_J
                rolling_v_idx = 1
            end
            is_v_valid, _ = check_is_sol_valid(bop, v_arr[rolling_v_idx])

            check_count += 1
            if check_count >= n_J
                if verbosity > 0
                    print("Failed to find a valid solution!\n")
                end
                break
            end
        end
        v .= v_arr[rolling_v_idx]

        # check if v is stable between iterations if it's a valid solution
        dv = v - prev_iter_v
        prev_iter_v .= v

        if LinearAlgebra.norm(dv) < tol && is_sol_valid
            is_converged = true
            break
        else
            if verbosity > 3
                print("Found new v: ")
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
            print("Success in $iter_count iterations\n")
        end
    end

    x .= v[bop.v_inds["x"]]

    # final feasibility check
    if !(all(bop.G(x) .≥ 0 - tol) && all(bop.g(x) .≥ 0 - tol))
        if verbosity > 0
            print("Something went wrong, this solution is not valid!\n")
        end
        is_sol_valid = false
    end

    is_success = is_converged && is_sol_valid

    (; x, is_success, iter_count)
end

function is_follower_KKT_satisfied(bop, v; tol=1e-6)
    x = @view v[bop.v_inds["x"]]
    λ = @view v[bop.v_inds["λ"]]

    ∇ₓf = zeros(bop.n₁ + bop.n₂)
    bop.deriv_funs.∇ₓf!(∇ₓf, x)
    ∇ₓg_vals = zeros(length(bop.deriv_funs.∇ₓg_rows))
    bop.deriv_funs.∇ₓg_vals!(∇ₓg_vals, x)
    ∇ₓg = sparse(bop.deriv_funs.∇ₓg_rows, bop.deriv_funs.∇ₓg_cols, ∇ₓg_vals)
    x₂_inds = bop.n₁+1:bop.n₁+bop.n₂
    is_stationary = all(isapprox.(∇ₓf[x₂_inds] - ∇ₓg[:, x₂_inds]' * λ, 0; atol=2 * tol))
    is_primal_feas = all(bop.g(x) .≥ 0 - tol) # primal feas
    is_dual_feas = all(λ .≥ 0 - tol) # dual feas
    is_complement = all(isapprox.(λ .* bop.g(x), 0; atol=2 * tol)) # complementarity

    return is_stationary && is_primal_feas && is_dual_feas && is_complement
end


function initialize_z!(v, bop; verbosity=0, is_using_PATH=false, is_using_HSL=false)
    # if BOPᵢ wasn't solved the low level solution may be invalid, and we have to call the follower nlp
    x₁ = @view v[bop.v_inds["x₁"]]
    x₂ = @view v[bop.v_inds["x₂"]]
    λ = zeros(bop.m₂)

    if is_using_PATH
        θ_out, status, _ = bop.solve_follower_KKT_mcp(x₁)
        x₂ = θ_out[1:bop.n₂]
        init_z_success = status == PATHSolver.MCP_Solved
    else
        x₂, λ, solvestat, _ = bop.solve_follower_nlp(x₁; x₂_init=x₂, is_using_HSL)
        init_z_success = solvestat == 0 || solvestat == 1 # accept Solve_Succeeded and Solved_To_Acceptable_Level
    end

    # if the feasible region of the follower is empty for x₁, this finds a bilevel feasible x to move to
    if !init_z_success
        if verbosity > 1
            print("Resetting x to a bilevel feasible point. \n")
        end

        x_feas, _, solvestat, _ = bop.find_bilevel_feas_pt(x_init=[x₁; x₂])
        bilevel_feas_success = solvestat == 0 || solvestat == 1

        if bilevel_feas_success
            x₁ .= x_feas[bop.v_inds["x₁"]]
            x₂ .= x_feas[bop.v_inds["x₂"]]

            if is_using_PATH
                θ_out, status, _ = bop.solve_follower_KKT_mcp(x₁)
                x₂ = θ_out[1:bop.n₂]
                init_z_success = status == PATHSolver.MCP_Solved
            else
                x₂, λ, solvestat, _ = bop.solve_follower_nlp(x₁; x₂_init=x₂, is_using_HSL)
                init_z_success = solvestat == 0 || solvestat == 1 # accept Solve_Succeeded and Solved_To_Acceptable_Level
            end
        else
            if verbosity > 0
                print("Failed to find a bilevel feasible point to re-attempt, the problem may be bilevel infeasible! Check your x_init, G(x), and g(x).\n")
            end
        end
    end

    if init_z_success
        v[bop.v_inds["x₁"]] .= x₁
        v[bop.v_inds["x₂"]] .= x₂
        v[bop.v_inds["λ"]] .= λ
        v[bop.v_inds["s"]] .= bop.g([x₁; x₂])
    end

    return init_z_success
end

function update_v!(v, bop, Ji_bounds, v_l, v_u, Gh_l, Gh_u, θ_l, θ_u, θ_init; is_using_PATH=false, is_using_HSL=false)
    is_BOPᵢ_solved::Bool = false
    if is_using_PATH
        θ_init[bop.θ_inds["v"]] .= v
        θ_l[bop.θ_inds["z"]] .= Ji_bounds.z_l
        θ_u[bop.θ_inds["z"]] .= Ji_bounds.z_u
        θ_l[bop.θ_inds["rₕ"]] .= Ji_bounds.h_l
        θ_u[bop.θ_inds["rₕ"]] .= Ji_bounds.h_u
        θ_out, status, _ = bop.solve_BOPᵢ_KKT_mcp(x_l=θ_l, x_u=θ_u, x_init=θ_init)
        v_out = θ_out[bop.θ_inds["v"]]
        is_BOPᵢ_solved = status == PATHSolver.MCP_Solved
    else
        v_l[bop.v_inds["z"]] .= Ji_bounds.z_l
        v_u[bop.v_inds["z"]] .= Ji_bounds.z_u
        Gh_l[bop.Gh_inds["h"]] .= Ji_bounds.h_l
        Gh_u[bop.Gh_inds["h"]] .= Ji_bounds.h_u
        v_out, _, solvestat = bop.solve_BOPᵢ_nlp(x_l=v_l, x_u=v_u, g_l=Gh_l, g_u=Gh_u, x_init=v, is_using_HSL=is_using_HSL)
        is_BOPᵢ_solved = solvestat == 0 || solvestat == 1
    end

    v .= v_out # even if it's not solved, we update v so we can try to initialize z again
    is_BOPᵢ_solved
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
    Gh = zeros(bop.m₁ + bop.mₕ)
    bop.Gh!(Gh, v)
    h = @view Gh[bop.Gh_inds["h"]]
    z = @view v[bop.v_inds["z"]]
    z_l = @view bop.v_l₀[bop.v_inds["z"]]
    z_u = @view bop.v_u₀[bop.v_inds["z"]]
    K = Dict{Int,Vector{Int}}()

    # note which constraints are active
    for j in 1:bop.mₕ
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
    ambigu_inds = [i for i in 1:bop.mₕ if length(K[i]) > 1]
    single_inds = setdiff(1:bop.mₕ, ambigu_inds)
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
    h_l = zeros(bop.mₕ)
    h_u = zeros(bop.mₕ)
    z_l = zeros(bop.mₕ)
    z_u = zeros(bop.mₕ)
    z_l₀ = bop.v_l₀[bop.v_inds["z"]]
    z_u₀ = bop.v_u₀[bop.v_inds["z"]]

    for j in J[1] # hⱼ inactive (positive), Λₕⱼ at lb
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
    for j in J[3] # hⱼ inactive (negative), Λₕⱼ at ub
        h_l[j] = -Inf
        h_u[j] = 0.0
        z_l[j] = z_u₀[j]
        z_u[j] = z_u₀[j]
    end
    for j in J[4] # hⱼ free, Λₕⱼ fixed
        h_l[j] = -Inf
        h_u[j] = Inf
        z_l[j] = z_l₀[j]
        z_u[j] = z_u₀[j]
    end

    (; h_l, h_u, z_l, z_u)
end
