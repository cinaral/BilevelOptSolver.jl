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
function solve_bop(bop; x_init=zeros(bop.nx), param=Float64[], tol=1e-6, fol_feas_set_tol=1e-4, fol_feas_set_tol_max=1e0, x_agree_tol=1e-4, max_iter=50, verbosity=0, init_solver="IPOPT", solver="IPOPT", conv_tol=1e-4, conv_dv_len=3, max_rand_restart_ct=10, do_force_hp_init=false, do_require_all_solved=true, do_require_strict_min=true, do_check_x_agreem=true, do_force_toggle=false, rng=MersenneTwister(123), x_init_min=fill(-1.0, bop.nx), x_init_max=fill(1.0, bop.nx), print_sigdigs=3)

    @assert length(x_init) >= bop.nx "Wrong x_init length!"

    # buffers
    v = zeros(bop.n)
    v[bop.inds.v["x"]] .= x_init[1:bop.nx]
    Λ = zeros(bop.m)
    v_correspond_Js = copy(v)

    is_necc_fol::Vector{Bool} = []
    is_sufc_fol::Vector{Bool} = []
    i_arr::Vector{Int64} = []
    v_arr::Vector{Vector{Float64}} = []
    Λ_arr::Vector{Vector{Float64}} = []
    J2_seen_arr::Vector{Set{Int64}} = []
    v_seen_arr::Vector{Vector{Float64}} = []
    Λ_seen_arr::Vector{Vector{Float64}} = []
    is_solved_seen_arr::Vector{Bool} = []

    # convergence stuff 
    is_sol_valid = false
    is_prev_v_set = false
    is_norm_dv_full = false
    init_counter = 0
    norm_dv_arr = zeros(conv_dv_len)
    norm_dv_cur_idx = 1
    prev_iter_v = zeros(bop.n)

    # restart stuff
    is_done = false
    is_converged = false
    is_initialized = false
    is_random_restart = false
    random_restart_count = 0
    status = "uninitialized"

    iter_count = 0

    while !is_done
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

        if is_random_restart
            if random_restart_count < max_rand_restart_ct
                is_random_restart = false
                random_restart_count += 1
                if verbosity > 0
                    print("Restarting with random init (attempt $random_restart_count)\n")
                end
                v[bop.inds.v["x"]] .= x_init_min .+ rand(rng, bop.nx, 1) .* (x_init_max - x_init_min)
                #@info v
                is_converged = false
                is_initialized = false
                is_prev_v_set = false
                norm_dv_arr .= 0.0
                norm_dv_cur_idx = 1
                empty!(J2_seen_arr)
                empty!(v_seen_arr)
                empty!(Λ_seen_arr)
                empty!(is_solved_seen_arr)
                status = "uninitialized"
                prev_iter_v = zeros(bop.n)
            else
                if verbosity > 0
                    print("Reached maximum number of random restarts! Terminating\n")
                end
                status = "max_random_restarts"
                is_done = true
                break
            end
        end

        #### initialize
        if !is_initialized
            is_initialized = initialize_z!(v, bop; param, verbosity, init_solver, tol, do_force_hp_init)

            norm_dv_arr = zeros(conv_dv_len)
            norm_dv_cur_idx = 1
            prev_iter_v = zeros(bop.n)

            if !is_initialized
                if verbosity > 0
                    print("Iteration $iter_count: Failed to (re)initialize!\n")
                end
                status = "init_fail"
                is_random_restart = true
                continue
            end
        end

        # compute feasible sets
        # by this point, v must at least satisfy the follower's problem
        follow_feas_Js = compute_follow_feas_ind_sets(bop, v; param, tol=fol_feas_set_tol, verbosity, do_force_toggle)
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
                    print("Relaxing follower feasible set tolerance $(round(fol_feas_set_tol,sigdigits=print_sigdigs)) and trying again...\n")
                end
                follow_feas_Js = compute_follow_feas_ind_sets(bop, v; tol=fol_feas_set_tol, do_force_toggle)
                tol_restart_count += 1
            end
        end

        n_J = length(follow_feas_Js)
        if n_J == 0
            status = "max_tol_relax"
            break
        end

        if verbosity > 2
            print("x=$(round.(v[bop.inds.v["x"]],sigdigits=print_sigdigs)), λ=$(round.(v[bop.inds.v["λ"]],sigdigits=print_sigdigs))\n")
        end

        if n_J > 1 && verbosity > 1
            print("Multiple feasible sets ($n_J) detected!\n")
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
                is_solved = solve_sbop_nlp!(v, Λ, bop, hl, hu, zl, zu; tol, solver, param, verbosity, solver_output=0)
                push!(J2_seen_arr, copy(J[2]))
                push!(v_seen_arr, copy(v))
                push!(Λ_seen_arr, copy(Λ))
                push!(is_solved_seen_arr, is_solved)
            end

            if is_solved
                dv = v_correspond_Js - v
                norm_dv = LinearAlgebra.norm(dv)

                if norm_dv > conv_tol
                    has_v_changed = true
                end

                is_fol_nec, is_fol_suf = check_follower_sol(v, bop; param, tol, do_require_strict_min)
                push!(is_necc_fol, is_fol_nec)
                push!(is_sufc_fol, is_fol_suf)
                push!(v_arr, copy(v))
                push!(Λ_arr, copy(Λ))
                push!(i_arr, i)

                if verbosity > 2
                    print("SBOP$i: Solved (follower nc: $is_fol_nec sc: $is_fol_suf) J2=$(J[2]), x=$(round.(v[bop.inds.v["x"]],sigdigits=print_sigdigs)), λ=$(round.(v[bop.inds.v["λ"]],sigdigits=print_sigdigs))\n")
                end

                #is_sbop_nec, is_sbop_suf = check_sbop_sol(v, Λ, bop, hl, hu, zl, zu; param, tol)
                #@info "sbop N $is_sbop_nec / S $is_sbop_suf"
            else
                if verbosity > 1
                    print("SBOP$i: FAILED J2: $(J[2]), x=$(round.(v[bop.inds.v["x"]],sigdigits=print_sigdigs)), λ=$(round.(v[bop.inds.v["λ"]],sigdigits=print_sigdigs))\n")
                end
            end
        end

        ## check if we need to re-initialize
        # if solved none or if none of the follower's are actually solved, we have to re-initialize
        if length(is_necc_fol) == 0 || !any(is_necc_fol)
            is_sol_valid = false
            if verbosity > 1
                print("Failed to solve any SBOPi or no valid follower solutions!\n")
            end
            is_random_restart = true
            continue
        end

        # update v based on cost among the solutions
        Fs = map(v -> bop.F([v[bop.inds.v["x"]]; param]), v_arr)
        fs = map(v -> bop.f([v[bop.inds.v["x"]]; param]), v_arr)
        F_ = Inf
        chosen_i = 0
        for (i, is_fol_nec) in enumerate(is_necc_fol)
            if !is_fol_nec
                continue
            end
            if any(is_sufc_fol) # if checking min, skip non-sufficient solutions unless there's none
                if !is_sufc_fol[i]
                    continue
                end
            end
            if verbosity > 4
                print("Considering SBOP$(i_arr[i]): F(x): $(round.(Fs[i],sigdigits=print_sigdigs)) f(x): $(round.(fs[i],sigdigits=print_sigdigs))\n")
            end
            if Fs[i] < F_  # choose the smallest available F value
                chosen_i = i_arr[i]
                F_ = Fs[i]
                v .= v_arr[i]
                Λ .= Λ_arr[i]
            end
        end

        if verbosity > 1
            print("Chose SBOP$chosen_i x=$(round.(v[bop.inds.v["x"]],sigdigits=print_sigdigs)), λ=$(round.(v[bop.inds.v["λ"]],sigdigits=print_sigdigs))\n")
        end

        #all(is_necc_SBOPi) # this will usually be false
        is_sol_valid = false
        is_all_solved = length(v_arr) == n_J
        if is_all_solved
            if verbosity > 1
                print("All SBOPi ($n_J) solved.\n")
            end
        end

        # if we solved all, there's a good chance we found a solution. now we can check optimality conditions
        if do_check_x_agreem
            is_sol_valid = true
            for (i, vv) in enumerate(v_arr)
                norm_x_err = LinearAlgebra.norm(vv[bop.inds.v["x"]] - v_correspond_Js[bop.inds.v["x"]]) # only checking x error
                if (norm_x_err ≤ x_agree_tol)
                    if !is_sufc_fol[i]
                        is_sol_valid = false
                    end
                end
            end
        else  # out of all solved solutions...
            if all(is_sufc_fol)
                is_sol_valid = true
            end
        end

        if do_require_all_solved && !is_all_solved
            is_sol_valid = false
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

            if !is_norm_dv_full && norm_dv_cur_idx == conv_dv_len
                is_norm_dv_full = true
            end

            # we check convergence by tracking conv_dv_len dv's, if all is less than conv tol we're done 
            if (is_norm_dv_full || norm_dv_cur_idx ≥ conv_dv_len) && all(norm_dv_arr .< conv_tol)
                if verbosity > 3
                    print("Converged in $iter_count iterations! (Last $conv_dv_len norm(dv) is less than conv tol $conv_tol)\n")
                end
                is_converged = true
                status = "converged"
            end

            if verbosity > 1
                print("norm(dv) $norm_dv\n")
            end

            norm_dv_cur_idx += 1
            if norm_dv_cur_idx > conv_dv_len
                norm_dv_cur_idx = 1
            end
        end

        if !has_v_changed
            # are there even any new v's to choose from?
            status = "v_unchanged"
            is_converged = true
        end

        # randomized restart
        if is_converged && is_sol_valid
            is_done = true
        elseif is_converged && !is_sol_valid
            if verbosity > 1
                print("Converged but invalid solution!\n")
            end
            is_random_restart = true
        end
    end

    x = zeros(bop.nx)
    λ = zeros(bop.m2)
    x .= v[bop.inds.v["x"]]
    λ .= v[bop.inds.v["λ"]]

    # final sanity check
    if is_sol_valid
        if !(all(bop.G([x; param]) .≥ 0 - tol) && all(bop.g([x; param]) .≥ 0 - tol))
            if verbosity > 0
                print("Something went VERY wrong!\n")
            end
            status = "very_wrong;"
            is_sol_valid = false
        end
    end
    if verbosity > 1
        print("\n")
    end

    (; is_sol_valid, x, λ, iter_count, status)
end

function initialize_z!(v, bop; param=Float64[], verbosity=0, init_solver="IPOPT", tol=1e-6, max_iter=100, do_force_hp_init=false)
    # if BOPᵢ wasn't solved the low level solution may be invalid, and we have to call the follower nlp
    if verbosity > 1
        print("Initializing...\n")
    end
    x1 = @view v[bop.inds.v["x1"]]
    x2 = @view v[bop.inds.v["x2"]]
    λ = zeros(bop.m2)
    success = false

    if !do_force_hp_init
        (; x2, λ, success) = solve_follower_nlp(bop, x1; x2_init=x2, solver=init_solver, tol, param)
    end

    # if failure it may be that the feasible region of the follower is empty for x₁
    if !success
        if verbosity > 1
            print("Resetting x to a bilevel feasible point using high-point relaxation...\n")
        end
        x, is_x_feasible = solve_high_point_nlp(bop; param, x_init=[x1; x2], tol, max_iter, solver=init_solver)

        x1 = @view x[bop.inds.x["x1"]]
        x2 = @view x[bop.inds.x["x2"]]

        if is_x_feasible # we try again
            v[bop.inds.v["x1"]] .= x1
            (; x2, λ, success) = solve_follower_nlp(bop, x1; x2_init=x2, solver=init_solver, tol, param)
        else
            if verbosity > 1
                print("Failed resetting x to a bilevel feasible point, the problem may be bilevel infeasible! Check your x_init, G(x), and g(x).\n")
            end
        end
    end

    if success
        v[bop.inds.v["x2"]] .= x2
        v[bop.inds.v["λ"]] .= λ
    end

    return success
end

function solve_follower_nlp(bop, x1; x2_init=zeros(bop.fol_nlp.n), param=Float64[], solver="IPOPT", tol=1e-6, max_iter=1000, verbosity=0, is_debug_on=false, solver_output=0)
    success = false

    if solver == "IPOPT"
        nlp = bop.fol_nlp
        x2, λ, solvestat, _ = BilevelOptSolver.solve_NLP(nlp.n, nlp.m, nlp.xl, nlp.xu, nlp.gl, nlp.gu, nlp.f, nlp.g!, nlp.∇ₓf!, nlp.∇ₓg_rows, nlp.∇ₓg_cols, nlp.∇ₓg_vals!, nlp.∇²ₓL_rows, nlp.∇²ₓL_cols, nlp.∇²ₓL_vals!; x_init=x2_init, λ_init=zeros(nlp.m), q=x1, p=param, verbosity, tol, max_iter, is_debug_on, print_level=solver_output, is_using_HSL=false)
        success = solvestat == 0
    elseif solver == "PATH"
        mcp = bop.fol_mcp
        z_init = zeros(mcp.n)
        z_init[bop.inds.z["x2"]] .= x2_init

        z, status, _ = BilevelOptSolver.solve_PATH(mcp.n, mcp.zl, mcp.zu, mcp.h!, mcp.∇h_rows, mcp.∇h_cols, mcp.∇h_vals!; q=x1, p=param, verbosity, is_debug_on, x_init=z_init, tol, max_iter, is_silent=(solver_output == 0))
        x2 = z[bop.inds.z["x2"]]
        λ = z[bop.inds.z["λ"]]

        success = status == PATHSolver.MCP_Solved
    end

    (; x2, λ, success)
end

function solve_sbop_nlp!(v, Λ, bop, hl, hu, zl, zu; tol=1e-6, max_iter=1000, param=Float64[], solver="IPOPT", verbosity=0, is_debug_on=false, solver_output=0)
    if solver == "IPOPT"
        nlp = bop.sbop_nlp

        Γl = fill(0.0, nlp.m)
        Γl[bop.inds.Γ["hl"]] .= hl
        Γl[bop.inds.Γ["zl"]] .= zl
        Γl[bop.inds.Γ["hu"]] .= -hu
        Γl[bop.inds.Γ["zu"]] .= -zu

        v_out, Λ_out, solvestat, _ = BilevelOptSolver.solve_NLP(nlp.n, nlp.m, nlp.xl, nlp.xu, Γl, nlp.gu, nlp.f, nlp.g!, nlp.∇ₓf!, nlp.∇ₓg_rows, nlp.∇ₓg_cols, nlp.∇ₓg_vals!, nlp.∇²ₓL_rows, nlp.∇²ₓL_cols, nlp.∇²ₓL_vals!; x_init=v, λ_init=Λ, p=param, verbosity, tol, max_iter, is_debug_on, print_level=solver_output, is_using_HSL=false)
        v .= v_out
        Λ .= Λ_out
        success = solvestat == 0
    elseif solver == "PATH"
        mcp = bop.sbop_mcp

        θ_init = zeros(mcp.n)
        θl = copy(mcp.zl)
        θ_init[bop.inds.θ["v"]] .= v[bop.inds.v["v"]]
        θ_init[bop.inds.θ["Λ"]] .= Λ
        θl[bop.inds.θ["rzl"]] .= zl
        θl[bop.inds.θ["rhl"]] .= hl
        θl[bop.inds.θ["rzu"]] .= -zu
        θl[bop.inds.θ["rhu"]] .= -hu

        θ, status, _ = BilevelOptSolver.solve_PATH(mcp.n, θl, mcp.zu, mcp.h!, mcp.∇h_rows, mcp.∇h_cols, mcp.∇h_vals!; p=param, verbosity, is_debug_on, x_init=θ_init, tol, max_iter, is_silent=(solver_output == 0))
        v[bop.inds.v["v"]] .= θ[bop.inds.θ["v"]]
        Λ .= θ[bop.inds.θ["Λ"]]

        success = status == PATHSolver.MCP_Solved
    end
end

function solve_high_point_nlp(bop; x_init=zeros(bop.hp_nlp.n), param=Float64[], tol=1e-6, max_iter=1000, solver="PATH", verbosity=0, is_debug_on=false, solver_output=0)
    if solver == "IPOPT"
        nlp = bop.hp_nlp
        x, _, solvestat, _ = BilevelOptSolver.solve_NLP(nlp.n, nlp.m, nlp.xl, nlp.xu, nlp.gl, nlp.gu, nlp.f, nlp.g!, nlp.∇ₓf!, nlp.∇ₓg_rows, nlp.∇ₓg_cols, nlp.∇ₓg_vals!, nlp.∇²ₓL_rows, nlp.∇²ₓL_cols, nlp.∇²ₓL_vals!; x_init, λ_init=zeros(nlp.m), p=param, verbosity, tol, max_iter, is_debug_on, print_level=solver_output, is_using_HSL=false)
        success = solvestat == 0 || solvestat == 1 # this is acceptable too since we only care about feasibility
    elseif solver == "PATH"
        mcp = bop.hp_mcp
        z_init = zeros(bop.hp_mcp.n)
        z_init[1:bop.nx] .= x_init

        z, status, _ = BilevelOptSolver.solve_PATH(mcp.n, mcp.zl, mcp.zu, mcp.h!, mcp.∇h_rows, mcp.∇h_cols, mcp.∇h_vals!; p=param, verbosity, is_debug_on, x_init=z_init, tol, max_iter, is_silent=(solver_output == 0))
        x = z[1:bop.nx]

        success = status == PATHSolver.MCP_Solved
    end

    (; x, success)
end

function check_follower_sol(v, bop; param=Float64[], tol=1e-5, do_require_strict_min=true)
    x = @view v[bop.inds.v["x"]]
    λ = @view v[bop.inds.v["λ"]]
    nlp = bop.fol_nlp

    check_nlp_sol(x, λ, nlp.n, nlp.m, zeros(nlp.m), nlp.g!, nlp.∇ₓf!, nlp.∇ₓg_size, nlp.∇ₓg_rows, nlp.∇ₓg_cols, nlp.∇ₓg_vals!, nlp.∇²ₓL_size, nlp.∇²ₓL_rows, nlp.∇²ₓL_cols, nlp.∇²ₓL_vals!; param, tol, do_require_strict_min)
end

function check_sbop_sol(v, Λ, bop, hl, hu, zl, zu; param=Float64[], tol=1e-5)
    nlp = bop.sbop_nlp
    Γl = fill(0.0, nlp.m)
    Γl[bop.inds.Γ["hl"]] .= hl
    Γl[bop.inds.Γ["zl"]] .= zl
    Γl[bop.inds.Γ["hu"]] .= -hu
    Γl[bop.inds.Γ["zu"]] .= -zu

    check_nlp_sol(v, Λ, nlp.n, nlp.m, Γl, nlp.g!, nlp.∇ₓf!, nlp.∇ₓg_size, nlp.∇ₓg_rows, nlp.∇ₓg_cols, nlp.∇ₓg_vals!, nlp.∇²ₓL_size, nlp.∇²ₓL_rows, nlp.∇²ₓL_cols, nlp.∇²ₓL_vals!; param, tol)
end

"""
n_actual refers to part of v that corresponds to x1 and x2, and without the λ part which always violates SC 
"""
function check_nlp_sol(x, λ, n, m, gl, g!, ∇ₓf!, ∇ₓg_size, ∇ₓg_rows, ∇ₓg_cols, ∇ₓg_vals!, ∇²ₓL_size, ∇²ₓL_rows, ∇²ₓL_cols, ∇²ₓL_vals!; param=Float64[], tol=1e-5, do_require_strict_min=true)
    g = zeros(m)
    ∇ₓf = zeros(n)
    ∇ₓg = sparse(∇ₓg_rows, ∇ₓg_cols, zeros(length(∇ₓg_rows)), ∇ₓg_size[1], ∇ₓg_size[2])
    ∇²ₓL = sparse(∇²ₓL_rows, ∇²ₓL_cols, zeros(length(∇²ₓL_rows)), ∇²ₓL_size[1], ∇²ₓL_size[2])
    g!(g, [x; param])
    ∇ₓf!(∇ₓf, [x; param])
    ∇ₓg_vals!(∇ₓg.nzval, [x; param])
    ∇²ₓL_vals!(∇²ₓL.nzval, [x; param], λ, 1.0)

    # necessary conditions
    is_stationary = all(isapprox.(∇ₓf - ∇ₓg' * λ, 0; atol=tol))
    #is_complement = all(isapprox.(g .* λ, 0; atol=1e1*tol)) # complementarity, the tolerance here is difficult
    is_complement = all(isapprox.(g, 0; atol=tol) .| isapprox.(λ, 0; atol=tol))
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
    is_sufficient = false

    # certification: there exists no direction that is both descent and feasible
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
        if n - r > 0
            Z = V[:, n-r+1:n]
            min_eig = eigmin(Z' * ∇²ₓL * Z)
        elseif n - r == 0 # no feasible directions exist
            is_sufficient = true
        end
    end

    if !is_sufficient
        if min_eig ≥ tol # SOSC
            is_sufficient = true # strict minimum
        elseif isapprox(min_eig, 0; atol=tol)
            if do_require_strict_min # SOSC
                is_sufficient = false
            else # SONC
                is_sufficient = true
            end
        end
    end

    is_necessary = is_stationary && is_complement && is_primal_feas && is_dual_feas
    is_sufficient = is_necessary && is_sufficient

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
function check_mcp_sol(bop, n, h, z, zl, zu; tol=1e-3, do_force_toggle=false)
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
            push!(Kj, 2) # case 2
            if do_force_toggle # Always also consider active constraints inactive
                if j > bop.n2
                    push!(Kj, 1) 
                end
            end
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
function compute_follow_feas_ind_sets(bop, v; param=Float64[], tol=1e-3, verbosity=0, do_force_toggle=false)
    Js = Vector{Dict{Int,Set{Int}}}()
    h = zeros(bop.nz)
    bop.fol_mcp.h!(h, [v; param])
    z = @view v[bop.inds.v["z"]]
    is_valid, K = check_mcp_sol(bop, bop.nz, h, z, bop.fol_mcp.zl, bop.fol_mcp.zu; tol, do_force_toggle)

    if !is_valid
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
    zl₀ = bop.fol_mcp.zl
    zu₀ = bop.fol_mcp.zu

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
