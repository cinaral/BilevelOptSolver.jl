"""
```
Verbosity:
    0: silent
    1: minimum: warnings and errors
    2: basic: iterations
    3: extended: number of index sets, infeasible index sets
    4: full: show all index sets
    5: function trace
    6: full: v trace
```
"""
function solve_bop(bop; x_init=zeros(bop.nx), param=zeros(bop.np), tol=1e-6, fol_feas_set_tol_max=1e0, x_agree_tol=1e-4, max_iter=50, verbosity=0, init_solver="IPOPT", solver="IPOPT", is_checking_x_agree=false, conv_tol=1e-4, norm_dv_len=10, conv_dv_len=2, is_checking_min=false, max_no_sols=2, max_inits=2, is_always_hp=false, is_forcing_inactive_inds=false, is_require_all_solved=false)

    if is_forcing_inactive_inds && is_require_all_solved
        if verbosity > 0
            print("is_forcing_inactive_inds and is_require_all_solved cannot be on at the same time, setting is_require_all_solved=false!\n")
        end
        is_require_all_solved = false
    end

    # buffers
    v = zeros(bop.n + bop.np)
    v[bop.inds.v["x"]] .= x_init[1:bop.nx]
    v[bop.inds.v["p"]] .= param
    Λ = zeros(bop.m)

    v_correspond_Js = copy(v)

    is_necc_fol::Vector{Bool} = []
    is_sufc_fol::Vector{Bool} = []
    #is_necc_SBOPi::Vector{Bool} = []
    #is_sufc_SBOPi::Vector{Bool} = []
    i_arr::Vector{Int64} = []
    v_arr::Vector{Vector{Float64}} = []
    Λ_arr::Vector{Vector{Float64}} = []
    J2_seen_arr::Vector{Set{Int64}} = []
    v_seen_arr::Vector{Vector{Float64}} = []
    Λ_seen_arr::Vector{Vector{Float64}} = []
    is_solved_seen_arr::Vector{Bool} = []
    ## restart stuff
    fol_feas_set_tol = 1e3 * copy(tol)

    # convergence stuff 
    is_sol_valid = false
    is_prev_v_set = false
    is_norm_dv_full = false
    no_sol_counter = 0
    init_counter = 0
    norm_dv_arr = zeros(norm_dv_len)
    chron_norm_dv_arr = copy(norm_dv_arr)
    norm_dv_cur_idx = 1
    prev_iter_v = zeros(bop.n)

    # these get used a lot, so we only set them up once
    if solver == "IPOPT"
        solve_SBOPi_IPOPT! = setup_solve_SBOPi_IPOPT(bop)
    elseif solver == "PATH"
        solve_SBOPi_PATH! = setup_solve_SBOPi_PATH(bop)
    end

    iter_count = 0
    is_converged = false
    is_initialized = false
    status = "uninitialized"

    while !is_converged
        if iter_count >= max_iter
            if verbosity > 0
                print("Max iterations ($iter_count) reached!\n")
            end
            status = "max_iter"
            break
        end
        iter_count += 1
        if verbosity > 1
            print("--Iteration $iter_count\n")
        end

        #### initialize
        if !is_initialized
            init_counter += 1
            if init_counter > max_inits
                status = "max_init"
                if verbosity > 0
                    print("Max initializations ($init_counter) reached!\n")
                end
                break
            end

            is_initialized = initialize_z!(v, bop; verbosity, init_solver, tol, is_always_hp)

            norm_dv_arr = zeros(norm_dv_len)
            chron_norm_dv_arr = copy(norm_dv_arr)
            norm_dv_cur_idx = 1
            prev_iter_v = zeros(bop.n)

            if !is_initialized
                if verbosity > 0
                    print("Iteration $iter_count: Failed to (re)initialize!\n")
                end
                status = "init_fail"
                break
            end
        end

        #### compute feasible sets
        # by this point, v must at least satisfy the follower's problem
        if verbosity > 5
            print("Computing feasible sets for the follower at v: ")
            display(v)
        end

        #v[bop.inds.v["λ"]] .= 0.
        follow_feas_Js = compute_follow_feas_ind_sets(bop, v; tol=fol_feas_set_tol, verbosity, is_forcing_inactive_inds)
        if length(follow_feas_Js) == 0
            if verbosity > 1
                print("Failed to compute follower feasible index sets: Not a valid solution! Maybe it's a tolerance issue?\n")
            end
            tol_restart_count = 1

            while length(follow_feas_Js) < 1
                fol_feas_set_tol = 10^(tol_restart_count) * tol

                if fol_feas_set_tol > fol_feas_set_tol_max
                    if verbosity > 1
                        print("Reached max follow set tol relaxations!\n")
                    end
                    fol_feas_set_tol = fol_feas_set_tol_max
                    break
                end
                if verbosity > 1
                    print("Relaxing follower feasible set tolerance $(round(fol_feas_set_tol,sigdigits=2)) and trying again...\n")
                end
                follow_feas_Js = compute_follow_feas_ind_sets(bop, v; tol=fol_feas_set_tol, is_forcing_inactive_inds)
                tol_restart_count += 1
            end
        end

        n_J = length(follow_feas_Js)
        if n_J == 0
            status = "max_tol_relax"
            break
        end

        if verbosity > 2
            print("x=$(v[bop.inds.v["x"]])\n")
        end

        if n_J > 1 && verbosity > 1
            print("Multiple feasible sets ($n_J) detected at!\n")
        end

        # solve for all feasible sets
        v_correspond_Js .= v
        has_v_changed = false

        empty!(is_necc_fol)
        empty!(is_sufc_fol)
        empty!(v_arr)
        empty!(Λ_arr)
        empty!(i_arr)

        for (i, J) in enumerate(follow_feas_Js)
            hl, hu, zl, zu = convert_J_to_h_z_bounds(J, bop)

            i_seen = findfirst(isequal(J[2]), J2_seen_arr)
            if !isnothing(i_seen)
                v .= v_seen_arr[i_seen]
                Λ .= Λ_seen_arr[i_seen]
                is_solved = is_solved_seen_arr[i_seen]
                if verbosity > 3
                    print("Already attempted: ")
                end
            else
                if solver == "IPOPT"
                    is_solved = solve_SBOPi_IPOPT!(v, Λ, hl, hu, zl, zu; tol, max_iter)
                elseif solver == "PATH"
                    is_solved = solve_SBOPi_PATH!(v, Λ, hl, hu, zl, zu; tol, max_iter)
                end
                push!(J2_seen_arr, copy(J[2]))
                push!(v_seen_arr, copy(v))
                push!(Λ_seen_arr, copy(Λ))
                push!(is_solved_seen_arr, is_solved)
            end

            if is_solved
                dv = v_correspond_Js - v
                norm_dv = LinearAlgebra.norm(dv)

                if norm_dv > conv_tol
                    #if verbosity > 3
                    #    print("New v=$v\n")
                    #end
                    has_v_changed = true
                end
                is_fol_nec, is_fol_suf = check_follower_sol(v, bop; tol=1e2 * tol)
                #is_SBOPi_nec, is_SBOPi_suf = check_SBOPi_sol(v, Λ, bop, hl, hu, zl, zu; tol=1e2 * tol) # is_SBOPi_suf is always false
                push!(is_necc_fol, is_fol_nec)
                push!(is_sufc_fol, is_fol_suf)
                #push!(is_necc_SBOPi, is_SBOPi_nec)
                #push!(is_sufc_SBOPi, is_SBOPi_suf)
                push!(v_arr, copy(v))
                push!(Λ_arr, copy(Λ))
                push!(i_arr, i)

                if verbosity > 2
                    print("SBOP$i: Solved (follower nc: $is_fol_nec sc: $is_fol_suf) J2: $(J[2]) x=$(v[bop.inds.v["x"]])\n")
                end

            else
                if verbosity > 1
                    print("SBOP$i: FAILED J1: $(J[1]) J2: $(J[2])\n")
                end
            end
        end

        #### check if we need to re-initialize
        # if solved none, we have to check if we have to re-initialize
        if length(is_necc_fol) == 0
            is_sol_valid = false
            no_sol_counter += 1
            if verbosity > 1
                print("Failed to solve any SBOPi!\n")
            end
            if no_sol_counter > max_no_sols
                status = "max_no_sols"
                break
            end

            is_necessary, _ = check_follower_sol(v, bop; tol=1e2 * tol)
            if !is_necessary
                is_initialized = false
                if verbosity > 0
                    print("Iteration $iter_count: Invalid follower solution, we must initialize again!\n")
                end
                is_prev_v_set = false
                norm_dv_arr .= 0.0
                norm_dv_cur_idx = 1
                status = "uninitialized"
            end
            continue
        end

        # if none of the follower's are actually solved, we still need to reinitialize
        if !any(is_necc_fol)
            is_initialized = false
            if verbosity > 0
                print("Iteration $iter_count: No valid follower solution, we must initialize again!\n")
            end
            is_prev_v_set = false
            norm_dv_arr .= 0.0
            norm_dv_cur_idx = 1
            status = "uninitialized"
            continue
        end
        no_sol_counter = 0

        # update v based on cost among the solutions
        Fs = map(v -> bop.F(v[bop.inds.v["x"]]), v_arr)
        fs = map(v -> bop.f(v[bop.inds.v["x"]]), v_arr)
        F_ = Inf
        chosen_i = 0
        for (i, is_fol_nec) in enumerate(is_necc_fol)
            if !is_fol_nec
                continue
            end
            if is_checking_min && any(is_sufc_fol) # if checking min, skip non-sufficient solutions unless there's none
                if !is_sufc_fol[i]
                    continue
                end
            end
            if verbosity > 4
                print("Considering SBOP$(i_arr[i]): F(x): $(Fs[i]) f(x): $(fs[i])...\n")
            end
            if Fs[i] < F_  # choose the smallest available F value
                chosen_i = i_arr[i]
                F_ = Fs[i]
                v .= v_arr[i]
                Λ .= Λ_arr[i]
            end
        end

        if verbosity > 1
            print("Chose SBOP$chosen_i x=$(v[bop.inds.v["x"]])\n")
        end

        #all(is_necc_SBOPi) # this will usually be false
        is_sol_valid = false
        is_all_solved = length(v_arr) == n_J
        if is_all_solved
            if verbosity > 0
                print("All SBOPi ($n_J) solved.\n")
            end
        end

        # if we solved all, there's a good chance we found a solution. now we can check optimality conditions
        if is_require_all_solved
            if is_all_solved
                if is_checking_min
                    if all(is_sufc_fol)
                        is_sol_valid = true
                    end
                else # ignore sufc
                    is_sol_valid = true
                end
            end
        else # out of all solved solutions...
            if is_checking_min
                if all(is_sufc_fol)
                    is_sol_valid = true
                end
            else # ignore sufc
                if all(is_necc_fol)
                    is_sol_valid = true
                end
            end
        end

        # optional check
        if is_checking_x_agree && is_sol_valid && length(v_arr) > 1
            for (i, vv) in enumerate(v_arr)
                if i < n_J
                    norm_x_err = LinearAlgebra.norm(v_arr[i+1][bop.inds.v["x"]] - vv[bop.inds.v["x"]]) # only checking x error
                    if (norm_x_err > x_agree_tol)
                        if verbosity > 1
                            print("SBOPi solutions disagree! norm x err: $norm_x_err\n")
                        end
                        is_sol_valid = false
                        break
                    end
                end
            end
        end

        # we check if the solution has converged even if the solution is not valid
        if !is_prev_v_set
            prev_iter_v .= v[1:bop.n]
            is_prev_v_set = true
        else
            dv = v[bop.inds.v["x"]] - prev_iter_v[bop.inds.v["x"]]
            prev_iter_v .= v[1:bop.n]

            norm_dv = LinearAlgebra.norm(dv)
            norm_dv_arr[norm_dv_cur_idx] = norm_dv

            if !is_norm_dv_full && norm_dv_cur_idx == norm_dv_len
                is_norm_dv_full = true
            end

            # in chronological order
            if !is_norm_dv_full
                chron_norm_dv_arr[end-norm_dv_cur_idx+1:end] .= norm_dv_arr[1:norm_dv_cur_idx]
            else
                if norm_dv_cur_idx < norm_dv_len
                    chron_norm_dv_arr .= [norm_dv_arr[norm_dv_cur_idx+1:end]; norm_dv_arr[1:norm_dv_cur_idx]]
                else
                    chron_norm_dv_arr .= norm_dv_arr
                end
            end

            # we check convergence by tracking norm_dv_len dv's, if all is less than conv tol we're done 
            if (is_norm_dv_full || norm_dv_cur_idx ≥ conv_dv_len) && all(chron_norm_dv_arr[end-conv_dv_len+1:end] .< conv_tol)
                if verbosity > 3
                    print("Converged in $iter_count iterations! (Last $conv_dv_len norm(dv) is less than conv tol $conv_tol)\n")
                end
                is_converged = true
                status = "converged"
                break
            elseif is_norm_dv_full && all(diff(chron_norm_dv_arr) .≤ -tol)
                # also it's worth checking if dv is slowly
                if verbosity > 3
                    print("last $norm_dv_len norm(dv) is monotonously decreasing without meeting conv tol, terminating\n")
                end
                status = "monoton_decreasing"
                break
            end

            if verbosity > 1
                print("norm(dv) $norm_dv\n")
            end

            norm_dv_cur_idx += 1
            if norm_dv_cur_idx > norm_dv_len
                norm_dv_cur_idx = 1
            end
        end

        if !has_v_changed
            # are there even any new v's to choose from?
            status = "v_unchanged"
            break
        end
    end

    x = @view v[bop.inds.v["x"]]
    λ = @view v[bop.inds.v["λ"]]
    # final sanity check
    #is_fol_necessary, is_fol_sufficient = check_follower_sol(v, bop; tol)

    if is_sol_valid
        if !(all(bop.G(x) .≥ 0 - tol) && all(bop.g(x) .≥ 0 - tol))
            if verbosity > 0
                print("Something went VERY wrong!\n")
            end
            status = "very_wrong"
            is_sol_valid = false
        end
    end
    if verbosity > 1
        print("\n")
    end

    (; is_sol_valid, x, λ, iter_count, status)
end

function initialize_z!(v, bop; verbosity=0, init_solver="IPOPT", tol=1e-6, max_iter=100, is_always_hp=false)
    # if BOPᵢ wasn't solved the low level solution may be invalid, and we have to call the follower nlp
    if verbosity > 1
        print("Initializing...\n")
    end
    x = zeros(bop.nx)
    x1 = @view v[bop.inds.v["x1"]]
    x2 = @view v[bop.inds.v["x2"]]
    λ = zeros(bop.m2)
    success = false

    if !is_always_hp
        (; x, λ, success) = solve_follower_nlp(bop, x1; x2_init=x2, solver=init_solver, tol, param=v[bop.inds.v["p"]])
    end
    # if failure it may be that the feasible region of the follower is empty for x₁
    if !success
        if verbosity > 0
            print("Resetting x to a bilevel feasible point using high-point relaxation...\n")
        end
        x, λ, is_x_feasible = solve_high_point_nlp(bop; x_init=x, tol=1e-6, max_iter=1000)
        x1 = @view x[bop.inds.x["x1"]]
        x2 = @view x[bop.inds.x["x2"]]

        if is_x_feasible # we try again
            (; x, λ, success) = solve_follower_nlp(bop, x1; x2_init=x2, solver=init_solver, tol, param=v[bop.inds.v["p"]])
        else
            if verbosity > 0
                print("Failed resetting x to a bilevel feasible point, the problem may be bilevel infeasible! Check your x_init, G(x), and g(x).\n")
            end
        end
    end

    if success
        v[bop.inds.v["x"]] .= x[bop.inds.v["x"]]
        v[bop.inds.v["λ"]] .= λ
    end

    return success
end

function solve_follower_nlp(bop, x1; x2_init=zeros(bop.n2), param=zeros(bop.np), solver="IPOPT", tol=1e-6, max_iter=1000)
    x = zeros(bop.nx + bop.np)
    x[bop.inds.x["x1"]] .= x1
    x[bop.inds.x["p"]] .= param
    λ = zeros(bop.m2)
    success = false

    if solver == "IPOPT"
        x2_l = fill(-Inf, bop.n2)
        x2_u = fill(Inf, bop.n2)
        gl = fill(0.0, bop.m2)
        gu = fill(Inf, bop.m2)

        function eval_f(x2::Vector{Float64})
            x[bop.inds.x["x2"]] .= x2
            bop.f(x)
        end
        function eval_g(g::Vector{Float64}, x2::Vector{Float64})
            x[bop.inds.x["x2"]] .= x2
            bop.g!(g, x)
        end
        function eval_∇ₓ₂f(∇ₓ₂f::Vector{Float64}, x2::Vector{Float64})
            x[bop.inds.x["x2"]] .= x2
            bop.∇ₓ₂f!(∇ₓ₂f, x)
        end
        function eval_∇ₓ₂g_vals(∇ₓ₂g_vals::Vector{Float64}, x2::Vector{Float64})
            x[bop.inds.x["x2"]] .= x2
            bop.∇ₓ₂g_vals!(∇ₓ₂g_vals, x)
        end
        function eval_∇²ₓ₂L2_vals(∇²ₓ₂L2_vals::Vector{Float64}, x2::Vector{Float64}, λ::Vector{Float64}, obj_factor::Float64)
            x[bop.inds.x["x2"]] .= x2
            bop.∇²ₓ₂L2_vals!(∇²ₓ₂L2_vals, x, λ, obj_factor)
        end
        solve = setup_nlp_solve_IPOPT(bop.n2, bop.m2, x2_l, x2_u, gl, gu, eval_f, eval_g, eval_∇ₓ₂f, bop.∇ₓ₂g_rows, bop.∇ₓ₂g_cols, eval_∇ₓ₂g_vals, bop.∇²ₓ₂L2_rows, bop.∇²ₓ₂L2_cols, eval_∇²ₓ₂L2_vals)

        x_out, λ_out, solvestat, _ = solve(; x_init=x2_init, tol, max_iter, print_level=0)

        x[bop.inds.x["x2"]] .= x_out
        λ .= λ_out
        success = solvestat == 0

    elseif solver == "PATH"
        v = zeros(bop.n)
        v[bop.inds.v["x"]] .= x
        z_init = zeros(bop.nz)
        z_init[bop.inds.z["x2"]] = x2_init

        function eval_F!(h, z::Vector{Float64})
            v[bop.inds.v["z"]] .= z
            bop.h!(h, v, 1.0)
        end
        function eval_J_vals!(∇h, z::Vector{Float64})
            v[bop.inds.v["z"]] .= z
            bop.∇h_vals!(∇h, v, 1.0)
        end
        solve = setup_mcp_solve_PATH(bop.nz, bop.zl₀, bop.zu₀, eval_F!, bop.∇h_rows, bop.∇h_cols, eval_J_vals!)

        z_out, status, _ = solve(; x_init=z_init, tol, max_iter, is_silent=true)

        x[bop.inds.x["x2"]] .= z_out[bop.inds.z["x2"]]
        λ .= z_out[bop.inds.z["λ"]]
        success = status == PATHSolver.MCP_Solved
    end

    (; x, λ, success)
end

function solve_high_point_nlp(bop; x_init=zeros(bop.nx + bop.np), tol=1e-6, max_iter=1000)
    xl = fill(-Inf, bop.n1 + bop.n2 + bop.np)
    xu = fill(Inf, bop.n1 + bop.n2 + bop.np)
    Gg_l = fill(0.0, bop.m1 + bop.m2)
    Gg_u = fill(Inf, bop.m1 + bop.m2)

    solve_high_point_nlp = setup_nlp_solve_IPOPT(bop.nx + bop.np, bop.m1 + bop.m2, xl, xu, Gg_l, Gg_u, bop.F, bop.Gg!, bop.∇ₓF!, bop.∇ₓGg_rows, bop.∇ₓGg_cols, bop.∇ₓGg_vals!, bop.∇²ₓL3_rows, bop.∇²ₓL3_cols, bop.∇²ₓL3_vals!)

    x, λ, solvestat, _ = solve_high_point_nlp(; x_init, tol, max_iter, print_level=0)

    success = solvestat == 0 || solvestat == 1 # this is acceptable too since we only care about feasibility
    (; x, λ, success)
end

function setup_solve_SBOPi_IPOPT(bop)
    vl = fill(-Inf, bop.n + bop.np)
    vu = fill(Inf, bop.n + bop.np)
    Γl = fill(0.0, bop.m)
    Γu₀ = fill(Inf, bop.m) # doesn't change
    solve_SBOPi_IPOPT = setup_nlp_solve_IPOPT(bop.n + bop.np, bop.m, vl, vu, Γl, Γu₀, bop.Fv, bop.Γ!, bop.∇ᵥF!, bop.∇ᵥΓ_rows, bop.∇ᵥΓ_cols, bop.∇ᵥΓ_vals!, bop.∇²ᵥL1_rows, bop.∇²ᵥL1_cols, bop.∇²ᵥL1_vals!)

    function solve_SBOPi!(v, Λ, hl, hu, zl, zu; tol=1e-6, max_iter=1000)
        Γl[bop.inds.Γ["hl"]] .= hl
        Γl[bop.inds.Γ["zl"]] .= zl
        Γl[bop.inds.Γ["hu"]] .= -hu
        Γl[bop.inds.Γ["zu"]] .= -zu
        v_out, Λ_out, solvestat, _ = solve_SBOPi_IPOPT(; gl=Γl, x_init=v, λ_init=Λ, tol, max_iter, print_level=0)
        v[1:bop.n] .= v_out[1:bop.n]
        Λ .= Λ_out
        success = solvestat == 0
        success
    end
    solve_SBOPi!
end

function setup_solve_SBOPi_PATH(bop)
    θ_init = zeros(bop.nθ)
    θl = copy(bop.θl₀)
    θu = copy(bop.θu₀)
    solve_SBOPi_PATH = setup_mcp_solve_PATH(bop.nθ, θl, θu, bop.Φ!, bop.∇Φ_rows, bop.∇Φ_cols, bop.∇Φ_vals!)

    function solve_SBOPi!(v, Λ, hl, hu, zl, zu; tol=1e-6, max_iter=5000)
        θ_init[bop.inds.θ["v"]] .= v
        θ_init[bop.inds.θ["Λ"]] .= Λ
        θl[bop.inds.θ["rzl"]] .= zl
        θl[bop.inds.θ["rhl"]] .= hl
        θl[bop.inds.θ["rzu"]] .= -zu
        θl[bop.inds.θ["rhu"]] .= -hu

        θ_out, status, _ = solve_SBOPi_PATH(; xl=θl, x_init=θ_init, tol, max_iter, is_silent=true)
        v .= θ_out[bop.inds.θ["v"]]
        Λ .= θ_out[bop.inds.θ["Λ"]]
        success = status == PATHSolver.MCP_Solved
        success
    end
    solve_SBOPi!
end

function check_follower_sol(v, bop; tol=1e-5)
    x = @view v[bop.inds.v["x"]]
    λ = @view v[bop.inds.v["λ"]]

    check_nlp_sol(x, λ, bop.n2, bop.m2, zeros(bop.m2), bop.g!, bop.∇ₓ₂f!, bop.∇ₓ₂g_size, bop.∇ₓ₂g_rows, bop.∇ₓ₂g_cols, bop.∇ₓ₂g_vals!, bop.∇²ₓ₂L2_size, bop.∇²ₓ₂L2_rows, bop.∇²ₓ₂L2_cols, bop.∇²ₓ₂L2_vals!; tol)
end

function check_SBOPi_sol(v, Λ, bop, hl, hu, zl, zu; tol=1e-5)
    Γl = fill(0.0, bop.m)
    Γl[bop.inds.Γ["hl"]] .= hl
    Γl[bop.inds.Γ["zl"]] .= zl
    Γl[bop.inds.Γ["hu"]] .= -hu
    Γl[bop.inds.Γ["zu"]] .= -zu

    check_nlp_sol(v, Λ, bop.n, bop.m, Γl, bop.Γ!, bop.∇ᵥF!, bop.∇ᵥΓ_size, bop.∇ᵥΓ_rows, bop.∇ᵥΓ_cols, bop.∇ᵥΓ_vals!, bop.∇²ᵥL1_size, bop.∇²ᵥL1_rows, bop.∇²ᵥL1_cols, bop.∇²ᵥL1_vals!; tol)
end

"""
n_actual refers to part of v that corresponds to x1 and x2, and without the λ part which always violates SC 
"""
function check_nlp_sol(x, λ, n, m, gl, g!, ∇ₓf!, ∇ₓg_size, ∇ₓg_rows, ∇ₓg_cols, ∇ₓg_vals!, ∇²ₓL_size, ∇²ₓL_rows, ∇²ₓL_cols, ∇²ₓL_vals!; tol=1e-5)
    #@assert(n == length(x))
    @assert(m == length(λ))
    @assert(m == length(gl))

    g = zeros(m)
    ∇ₓf = zeros(n)
    ∇ₓg = sparse(∇ₓg_rows, ∇ₓg_cols, zeros(length(∇ₓg_rows)), ∇ₓg_size[1], ∇ₓg_size[2])
    ∇²ₓL = sparse(∇²ₓL_rows, ∇²ₓL_cols, zeros(length(∇²ₓL_rows)), ∇²ₓL_size[1], ∇²ₓL_size[2])
    g!(g, x)
    ∇ₓf!(∇ₓf, x)
    ∇ₓg_vals!(∇ₓg.nzval, x)
    ∇²ₓL_vals!(∇²ₓL.nzval, x, λ, 1.0)

    # necessary conditions
    is_stationary = all(isapprox.(∇ₓf - ∇ₓg' * λ, 0; atol=tol))
    is_complement = all(isapprox.(g .* λ, 0; atol=tol)) # complementarity
    is_primal_feas = all(g .≥ gl .- tol)
    is_dual_feas = all(λ .≥ 0 - tol) # dual feas
    #is_necessary = is_stationary && is_complement && is_primal_feas && is_dual_feas



    # strongly active constraints
    active_js = []
    for i in 1:m
        if isapprox(g[i], 0; atol=tol) && λ[i] ≥ tol
            push!(active_js, i)
        end
    end

    min_eig = -Inf

    if isempty(active_js) # unconstrained problem
        if n > 1
            min_eig = eigmin(Matrix(∇²ₓL))
        else
            min_eig = ∇²ₓL[1]
        end
    else
        C = Matrix(∇ₓg[active_js, :])
        _, _, V = svd(C; full=true)

        r = rank(C)
        # 2025-08-04 TODO: what about when nullspace is empty?
        if n - r > 0
            Z = V[:, n-r+1:n]
            min_eig = eigmin(Z' * ∇²ₓL * Z)
        end
    end

    is_sufficient = false
    if min_eig ≥ tol
        is_sufficient = true # strict minimum
    elseif isapprox(min_eig, 0; atol=tol)
        is_sufficient = true # not strict minimum
    end

    #Main.@infiltrate
    is_necessary = is_stationary && is_complement && is_primal_feas && is_dual_feas
    is_sufficient = is_necessary && is_sufficient

    #Main.@infiltrate

    (; is_necessary, is_sufficient)
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
function check_mcp_sol(bop, n, h, z, zl, zu; tol=1e-3, is_forcing_inactive_inds=false)
    @assert(n == length(h))
    @assert(n == length(z))
    @assert(n == length(zl))
    @assert(n == length(zu))
    K = Dict{Int,Vector{Int}}()

    # note which constraints are active
    for j in 1:n
        Kj = Int[]

        if isapprox(zl[j], zu[j]; atol=2 * tol)
            push!(Kj, 4) # case 4
        elseif tol ≤ h[j] && z[j] < zl[j] + tol
            push!(Kj, 1) # case 1
        elseif -tol < h[j] < tol && z[j] < zl[j] + tol
            push!(Kj, 1) # case 1 OR  
            push!(Kj, 2) # case 2
        elseif -tol < h[j] < tol && zl[j] + tol ≤ z[j] ≤ zu[j] - tol
            # without this, we fail to find most sols
            if is_forcing_inactive_inds
                if j > bop.n2 # should be cleaner
                    #Main.@infiltrate
                    push!(Kj, 1) # case 1 # IF ACTIVE, also consider inactive...
                end
            end
            push!(Kj, 2) # case 2
        elseif -tol < h[j] < tol && zu[j] - tol < z[j]
            push!(Kj, 2) # case 2 OR 
            push!(Kj, 3) # case 3
        elseif h[j] ≤ tol && zu[j] - tol < z[j]
            push!(Kj, 2) # case 3
        end

        K[j] = Kj
    end
    is_valid = !any(isempty.(Kj for Kj in values(K)))
    (; is_valid, K)
end

# for the follower's problem
# check_nlp_sol(x, λ, bop.n2, bop.m2, zeros(bop.m2), bop.g!, bop.∇ₓ₂f!, bop.∇ₓ₂g_size, bop.∇ₓ₂g_rows, bop.∇ₓ₂g_cols, bop.∇ₓ₂g_vals!, bop.∇²ₓ₂L2_size, bop.∇²ₓ₂L2_rows, bop.∇²ₓ₂L2_cols, bop.∇²ₓ₂L2_vals!)
function compute_follow_feas_ind_sets(bop, v; tol=1e-3, verbosity=0, is_forcing_inactive_inds=false)
    Js = Vector{Dict{Int,Set{Int}}}()
    h = zeros(bop.nz)
    bop.h!(h, v)
    z = @view v[bop.inds.v["z"]]
    zl = bop.zl₀
    zu = bop.zu₀
    is_valid, K = check_mcp_sol(bop, bop.nz, h, z, zl, zu; tol, is_forcing_inactive_inds)

    #K = Dict{Int,Vector{Int}}()
    #for j in 1:bop.nz
    #    Kj = Int[]
    #    if j <= bop.n2
    #        push!(Kj, 2) # case 2
    #    else
    #        if h[j] > tol
    #            push!(Kj, 1) # case 1
    #            push!(Kj, 2) # case 2
    #        else
    #            push!(Kj, 1) # case 1
    #        end
    #    end

    #    K[j] = Kj
    #end
    #is_valid = !any(isempty.(Kj for Kj in values(K)))

    #Main.@infiltrate

    if !is_valid
        #if verbosity > 0
        #    print("Failed to compute follower feasible index sets: Not a valid solution!\n")
        #end
        return Js
    end

    # enumerate the ambiguous indexes
    ambigu_inds = [i for i in 1:bop.nz if length(K[i]) > 1]
    single_inds = setdiff(1:bop.nz, ambigu_inds)
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

function convert_J_to_h_z_bounds(J, bop; tol=1e-6)
    hl = zeros(bop.nz)
    hu = zeros(bop.nz)
    zl = zeros(bop.nz)
    zu = zeros(bop.nz)
    zl₀ = bop.zl₀
    zu₀ = bop.zu₀

    for j in J[1] # hⱼ inactive (strictly positive), zⱼ at lb
        hl[j] = 0.0
        hu[j] = Inf
        zl[j] = zl₀[j]
        zu[j] = zl₀[j]
    end
    for j in J[2] # hⱼ active (l ≠ u)
        hl[j] = 0.0
        hu[j] = 0.0
        zl[j] = zl₀[j]
        zu[j] = zu₀[j]
    end
    for j in J[3] # hⱼ inactive (strictly negative), zⱼ at ub
        hl[j] = -Inf
        hu[j] = 0.0
        zl[j] = zu₀[j]
        zu[j] = zu₀[j]
    end
    for j in J[4] # hⱼ free, zⱼ fixed (l = u)
        hl[j] = -Inf
        hu[j] = Inf
        zl[j] = zl₀[j]
        zu[j] = zu₀[j]
    end

    (; hl, hu, zl, zu)
end
