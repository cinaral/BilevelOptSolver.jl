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
function solve_bop(bop; x_init=zeros(bop.nx), tol=1e-6, fol_feas_set_tol_max=1e-0, x_agree_tol=1e-3, max_iter=10, verbosity=0, init_solver="IPOPT", solver="IPOPT", is_checking_x_agree=true, conv_tol=1e-3, norm_dv_len=10, conv_dv_len=2, is_checking_min=true)
    status = "ok"

    # buffer
    #x = copy(x_init)
    v = zeros(bop.n)
    v[bop.v_inds["x"]] .= x_init[1:bop.nx]
    Λ = zeros(bop.m)

    vΛ_J_inds::Vector{Int64} = [] # indexes v and Λ array wrt J
    v_arr::Vector{Vector{Float64}} = []
    Λ_arr::Vector{Vector{Float64}} = []
    #vΛ_status::Vector{Int64} = [] # 0 solved, 1 feasible

    iter_count = 0
    is_converged = false
    is_initialized = false

    ## restart stuff
    fol_feas_set_tol = copy(tol)

    # convergence stuff 
    is_sol_valid = false
    is_prev_v_set = false
    is_norm_dv_full = false
    norm_dv_arr = zeros(norm_dv_len)
    chron_norm_dv_arr = copy(norm_dv_arr)
    norm_dv_cur_idx = 1
    prev_iter_v = zeros(bop.n)


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

        if !is_initialized
            is_initialized = initialize_z!(v, bop; verbosity, init_solver, tol)

            if !is_initialized
                if verbosity > 0
                    print("Iteration $iter_count: Failed to (re)initialize!\n")
                end
                status = "init_fail"
                break
            end
        end

        # by this point, v must at least satisfy the follower's problem
        if verbosity > 5
            print("Computing feasible sets for the follower at v: ")
            display(v)
        end

        follow_feas_Js = compute_follow_feas_ind_sets(bop, v; tol=fol_feas_set_tol, verbosity)

        n_J = length(follow_feas_Js)

        if n_J == 0
            if verbosity > 1
                print("Could not compute follower's feasible sets despite satisfying follower's KKT! Maybe it's a tolerance issue?\n")
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
                follow_feas_Js = compute_follow_feas_ind_sets(bop, v; tol=fol_feas_set_tol)
                tol_restart_count += 1
            end
        end

        n_J = length(follow_feas_Js)
        if n_J == 0
            status = "max_tol_relax"
            break
        end

        fol_feas_set_tol = tol # reset fol_feas_set_tol

        if n_J > 1 && verbosity > 1
            print("Multiple feasible sets ($n_J) detected!\n")
        end

        empty!(vΛ_J_inds)
        empty!(v_arr)
        empty!(Λ_arr)

        for i in 1:n_J
            Ji = follow_feas_Js[i]
            hl, hu, zl, zu = convert_J_to_h_z_bounds(Ji, bop)

            is_solved = solve_SBOPi_nlp!(v, Λ, bop, hl, hu, zu, zl; solver, tol)

            is_valid = false
            if is_solved
                #is_SBOPi_nec, is_SBOPi_suf = check_SBOPi_sol(v, Λ, bop, hl, hu, zu, zl)
                #@info "SBOPi nec $is_SBOPi_nec suf $is_SBOPi_suf"
                is_fol_nec, is_fol_suf = check_follower_sol(v, bop)

                if is_checking_min && is_fol_nec && is_fol_suf
                    is_valid = true
                elseif !is_checking_min && is_fol_nec
                    is_valid = true
                end
            end

            #Main.@infiltrate
            if is_valid
                push!(vΛ_J_inds, i)
                push!(v_arr, copy(v))
                push!(Λ_arr, copy(Λ))

                if verbosity > 2
                    print("i=$i: Valid J1/2/3/4: $(Ji[1])/$(Ji[2])/$(Ji[3])/$(Ji[4])\n")
                end
            else
                if verbosity > 1
                    print("i=$i: NOT VALID J1/2/3/4: $(Ji[1])/$(Ji[2])/$(Ji[3])/$(Ji[4])\n")
                end
            end
        end

        if length(vΛ_J_inds) == n_J && n_J > 0
            is_sol_valid = true
        else
            is_sol_valid = false
        end

        if is_sol_valid
            if verbosity > 4
                print("All $n_J follower solutions have valid v and Λ.\n")
            end
        end

        if length(vΛ_J_inds) == 0
            #Main.@infiltrate
            if verbosity > 1
                print("Failed to find any valid solutions!\n")
            end
            status = "no_valid_sols"
            break
        end

        # optional check
        if is_checking_x_agree && is_sol_valid
            for (i, vv) in enumerate(v_arr)
                if i < n_J
                    norm_x_err = LinearAlgebra.norm(v_arr[i+1][bop.v_inds["x"]] - vv[bop.v_inds["x"]]) # only checking x error

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

        #n_Js = map(v -> length(compute_follow_feas_ind_sets(bop, v; tol=fol_feas_set_tol)), v_arr)
        #Fs = map(v -> bop.F(v[bop.v_inds["x"]]), v_arr)
        #fs = map(v -> bop.f(v[bop.v_inds["x"]]), v_arr)
        #Gs = mapreduce(v -> bop.G(v[bop.v_inds["x"]])', vcat, v_arr)
        #gs = mapreduce(v -> bop.g(v[bop.v_inds["x"]])', vcat, v_arr)
        #max_n_J = maximum(n_Js)
        #fs_sorted_inds = sortperm(fs)
        #@info "$n_Js $Fs $fs"

        F = Inf
        for (i, _) in enumerate(v_arr)
            F_ = bop.F(v_arr[i][bop.v_inds["x"]])   # choose the smallest available F value

            if F_ < F
                F = F_
                v .= v_arr[i]
                Λ .= Λ_arr[i]
            end
        end

        if is_initialized
            is_necessary, is_sufficient = check_follower_sol(v, bop)
            if is_checking_min
                is_initialized = is_necessary && is_sufficient
            else
                is_initialized = is_necessary
            end

            if !is_initialized
                if verbosity > 0
                    print("Iteration $iter_count: Follower is not at minimum, we must initialize again!\n")
                end
                continue
            end
        end

        # by this point is_sol_valid = true: v's are sols and agree, Λs are feasible, all we have to check now is if the solution has converged
        if !is_prev_v_set
            prev_iter_v .= v
            is_prev_v_set = true
        else
            dv = v[bop.v_inds["x"]] - prev_iter_v[bop.v_inds["x"]]
            prev_iter_v .= v
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
                break
            elseif is_norm_dv_full && all(diff(chron_norm_dv_arr) .≤ -tol)
                # also it's worth checking if dv is slowly
                if verbosity > 3
                    print("last $norm_dv_len norm(dv) is monotonously decreasing without meeting conv tol, terminating\n")
                end
                is_dv_mon_decreasing = true
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

        if verbosity > 5
            print("v = ")
            display(v)
        end
    end

    if is_converged
        if verbosity > 1
            print("\n")
        end
    end

    x = @view v[bop.v_inds["x"]]
    λ = @view v[bop.v_inds["λ"]]

    # final sanity check
    is_fol_necessary, is_fol_sufficient = check_follower_sol(v, bop)

    if is_sol_valid
        if !is_fol_necessary || (is_checking_min && !is_fol_sufficient) || !(all(bop.G(x) .≥ 0 - tol) && all(bop.g(x) .≥ 0 - tol))
            if verbosity > 0
                print("Something went VERY wrong!\n")
            end
            status = "very_wrong"
            is_sol_valid = false
        end
    end

    (; x, status, iter_count)
end

function initialize_z!(v, bop; verbosity=0, init_solver="IPOPT", tol=1e-6, max_iter=100)
    # if BOPᵢ wasn't solved the low level solution may be invalid, and we have to call the follower nlp
    if verbosity > 1
        print("Initializing...\n")
    end
    x = zeros(bop.nx)
    x1 = @view v[bop.v_inds["x1"]]
    x2 = @view v[bop.v_inds["x2"]]
    λ = zeros(bop.m2)
    (; x, λ, success) = solve_follower_nlp(bop, x1; x2_init=x2, solver=init_solver, tol)

    # if failure it may be that the feasible region of the follower is empty for x₁
    if !success
        if verbosity > 0
            print("Resetting x to a bilevel feasible point using high-point relaxation...\n")
        end
        x, λ, is_x_feasible = solve_high_point_nlp(bop; x_init=x, tol=1e-6, max_iter=1000)
        x1 = @view x[bop.x_inds["x1"]]
        x2 = @view x[bop.x_inds["x2"]]

        if is_x_feasible # we try again
            (; x, λ, success) = solve_follower_nlp(bop, x1; x2_init=x2, solver=init_solver, tol)
        else
            if verbosity > 0
                print("Failed resetting x to a bilevel feasible point, the problem may be bilevel infeasible! Check your x_init, G(x), and g(x).\n")
            end
        end
    end

    if success
        v[bop.v_inds["x"]] .= x
        v[bop.v_inds["λ"]] .= λ
    end

    return success
end

function solve_high_point_nlp(bop; x_init=zeros(bop.nx), tol=1e-6, max_iter=1000)
    xl = fill(-Inf, bop.n1 + bop.n2)
    xu = fill(Inf, bop.n1 + bop.n2)
    Gg_l = fill(0.0, bop.m1 + bop.m2)
    Gg_u = fill(Inf, bop.m1 + bop.m2)

    solve_high_point_nlp = setup_nlp_solve_IPOPT(bop.nx, bop.m1 + bop.m2, xl, xu, Gg_l, Gg_u, bop.F, bop.Gg!, bop.∇ₓF!, bop.∇ₓGg_rows, bop.∇ₓGg_cols, bop.∇ₓGg_vals!, bop.∇²ₓL₃_rows, bop.∇²ₓL₃_cols, bop.∇²ₓL₃_vals!)

    x, λ, solvestat, _ = solve_high_point_nlp(; x_init, tol, max_iter, print_level=0)

    success = solvestat == 0 || solvestat == 1 # this is acceptable too since we only care about feasibility
    (; x, λ, success)
end

function solve_follower_nlp(bop, x1; x2_init=zeros(bop.n2), solver="IPOPT", tol=1e-6, max_iter=1000)
    x = zeros(bop.nx)
    x[bop.x_inds["x1"]] .= x1
    λ = zeros(bop.m2)
    success = false

    if solver == "IPOPT"
        x2_l = fill(-Inf, bop.n2)
        x2_u = fill(Inf, bop.n2)
        gl = fill(0.0, bop.m2)
        gu = fill(Inf, bop.m2)

        function eval_f(x2::Vector{Float64})
            x[bop.x_inds["x2"]] .= x2
            bop.f(x)
        end
        function eval_g(g::Vector{Float64}, x2::Vector{Float64})
            x[bop.x_inds["x2"]] .= x2
            bop.g!(g, x)
        end
        function eval_∇ₓ₂f(∇ₓ₂f::Vector{Float64}, x2::Vector{Float64})
            x[bop.x_inds["x2"]] .= x2
            bop.∇ₓ₂f!(∇ₓ₂f, x)
        end
        function eval_∇ₓ₂g_vals(∇ₓ₂g_vals::Vector{Float64}, x2::Vector{Float64})
            x[bop.x_inds["x2"]] .= x2
            bop.∇ₓ₂g_vals!(∇ₓ₂g_vals, x)
        end
        function eval_∇²ₓ₂L₂_vals(∇²ₓ₂L₂_vals::Vector{Float64}, x2::Vector{Float64}, λ::Vector{Float64}, obj_factor::Float64)
            x[bop.x_inds["x2"]] .= x2
            bop.∇²ₓ₂L₂_vals!(∇²ₓ₂L₂_vals, x, λ, obj_factor)
        end
        solve = setup_nlp_solve_IPOPT(bop.n2, bop.m2, x2_l, x2_u, gl, gu, eval_f, eval_g, eval_∇ₓ₂f, bop.∇ₓ₂g_rows, bop.∇ₓ₂g_cols, eval_∇ₓ₂g_vals, bop.∇²ₓ₂L₂_rows, bop.∇²ₓ₂L₂_cols, eval_∇²ₓ₂L₂_vals)

        x_out, λ_out, solvestat, _ = solve(; x_init=x2_init, tol, max_iter, print_level=0)

        x[bop.x_inds["x2"]] .= x_out
        λ .= λ_out
        success = solvestat == 0 # || solvestat == 1

    elseif solver == "PATH"
        v = zeros(bop.n)
        v[bop.v_inds["x"]] .= x
        z_init = zeros(bop.nz)
        z_init[bop.z_inds["x2"]] = x2_init

        function eval_F!(h, z::Vector{Float64})
            v[bop.v_inds["z"]] .= z
            bop.h!(h, v, 1.0)
        end
        function eval_J_vals!(∇h, z::Vector{Float64})
            v[bop.v_inds["z"]] .= z
            bop.∇h_vals!(∇h, v, 1.0)
        end
        solve = setup_mcp_solve_PATH(bop.nz, bop.zl₀, bop.zu₀, eval_F!, bop.∇h_rows, bop.∇h_cols, eval_J_vals!)

        z_out, status, _ = solve(; x_init=z_init, tol, max_iter, is_silent=true)

        x[bop.x_inds["x2"]] .= z_out[bop.z_inds["x2"]]
        λ .= z_out[bop.z_inds["λ"]]
        success = status == PATHSolver.MCP_Solved
    end

    (; x, λ, success)
end

# 2025-07-10 in order to reduce allocations we could do...
#function setup_SBOP_nlp_IPOPT()
#end

#function setup_SBOP_nlp_PATH()
#end

function solve_SBOPi_nlp!(v, Λ, bop, hl, hu, zu, zl; solver="IPOPT", tol=1e-6, max_iter=1000)
    #v = zeros(bop.n)
    #v[bop.v_inds["x"]] .= x_init
    #Λ = zeros(bop.m)
    success = false

    if solver == "IPOPT"
        vl = fill(-Inf, bop.n)
        vu = fill(Inf, bop.n)
        Γl = fill(0.0, bop.m)
        Γl[bop.Γ_inds["hl"]] .= hl
        Γl[bop.Γ_inds["hu"]] .= -hu
        Γl[bop.Γ_inds["zl"]] .= zl
        Γl[bop.Γ_inds["zu"]] .= -zu
        Γu₀ = fill(Inf, bop.m) # doesn't change
        solve = setup_nlp_solve_IPOPT(bop.n, bop.m, vl, vu, Γl, Γu₀, bop.Fv, bop.Γ!, bop.∇ᵥF!, bop.∇ᵥΓ_rows, bop.∇ᵥΓ_cols, bop.∇ᵥΓ_vals!, bop.∇²ᵥL₁_rows, bop.∇²ᵥL₁_cols, bop.∇²ᵥL₁_vals!)

        v_out, Λ_out, solvestat, _ = solve(; gl=Γl, x_init=v, λ_init=Λ, tol, max_iter, print_level=0)
        success = solvestat == 0 # || solvestat == 1

    elseif solver == "PATH"
        θ_init = zeros(bop.nθ)
        θ_init[bop.θ_inds["v"]] .= v
        θ_init[bop.θ_inds["Λ"]] .= Λ
        θl = copy(bop.θl₀)
        θl[bop.θ_inds["Λhl"]] .= zl
        θl[bop.θ_inds["Λzl"]] .= hl
        θu = copy(bop.θu₀)
        θu[bop.θ_inds["Λhu"]] .= zu
        θu[bop.θ_inds["Λzu"]] .= hu

        solve = setup_mcp_solve_PATH(bop.nθ, θl, θu, bop.Φ!, bop.∇Φ_rows, bop.∇Φ_cols, bop.∇Φ_vals!)

        θ_out, status, _ = solve(; xl=θl, xu=θu, x_init=θ_init, tol, is_silent=true)
        v_out = θ_out[bop.θ_inds["v"]]
        Λ_out = θ_out[bop.θ_inds["Λ"]]
        success = status == PATHSolver.MCP_Solved
    end

    if success
        v .= v_out
        Λ .= Λ_out
    end
    success
end

function check_nlp_sol(x, λ, n, m, gl, g!, ∇ₓf!, ∇ₓg_size, ∇ₓg_rows, ∇ₓg_cols, ∇ₓg_vals!, ∇²ₓL_size, ∇²ₓL_rows, ∇²ₓL_cols, ∇²ₓL_vals!; tol=1e-3)
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
    is_stationary = all(isapprox.(∇ₓf - ∇ₓg' * λ, 0; atol=2 * tol))
    is_complement = isapprox(g' * λ, 0; atol=2 * tol) # complementarity
    is_primal_feas = all(g .≥ gl .- tol)
    is_dual_feas = all(λ .≥ 0 - tol) # dual feas
    is_necessary = is_stationary && is_complement && is_primal_feas && is_dual_feas

    # sufficient conditions
    # TODO 2025-07-10: tangent space is wrong
    if is_necessary
        is_sufficient = true
        act_inds = isapprox.(g, 0; atol=2 * tol)
        #act_inds = all(λ .≥ 0 + tol) # strictly active indices??
        if any(act_inds)
            for j in findall(act_inds .> 0)
                for k in 1:n
                    if isapprox(∇ₓg[j, k], 0; atol=2 * tol)
                        w = zeros(n)
                        w[k] = 1.0
                    else
                        w = ones(n)
                        w[k] = 0.0
                        offset = w' * ∇ₓg[j, :]
                        w[k] = -offset / ∇ₓg[j, k]
                    end

                    if !isapprox(w' * ∇ₓg[j, :], 0; atol=2 * tol)
                        #Main.@infiltrate
                    end
                    # sanity checks, there might be a smarter way to get the tangent cone 
                    @assert(isapprox(w' * ∇ₓg[j, :], 0; atol=2 * tol)) # w' * ∇ₓgⱼ = 0 for j∈Jᵢ⁺

                    # w ≠ 0 && λ ≠ 0 && w' * ∇²ₓL * w > 0 
                    if all(isapprox.(w, 0; atol=2 * tol)) || isapprox(λ[j], 0; atol=2 * tol) || w' * ∇²ₓL * w ≤ 0.0 
                        is_sufficient = false
                    end
                end
            end
        else
            # if there are no active indices this is the unconstrained problem
            is_sufficient = isposdef(∇²ₓL)
        end
    else
        is_sufficient = false # no need to check the sufficiency condition if necessary cond is not met
    end

    (; is_necessary, is_sufficient)
end

function check_follower_sol(v, bop)
    x = @view v[bop.v_inds["x"]]
    λ = @view v[bop.v_inds["λ"]]

    check_nlp_sol(x, λ, bop.n2, bop.m2, zeros(bop.m2), bop.g!, bop.∇ₓ₂f!, bop.∇ₓ₂g_size, bop.∇ₓ₂g_rows, bop.∇ₓ₂g_cols, bop.∇ₓ₂g_vals!, bop.∇²ₓ₂L₂_size, bop.∇²ₓ₂L₂_rows, bop.∇²ₓ₂L₂_cols, bop.∇²ₓ₂L₂_vals!)
end

function check_SBOPi_sol(v, Λ, bop, hl, hu, zu, zl)
    Γl = fill(0.0, bop.m)
    Γl[bop.Γ_inds["hl"]] .= hl
    Γl[bop.Γ_inds["hu"]] .= -hu
    Γl[bop.Γ_inds["zl"]] .= zl
    Γl[bop.Γ_inds["zu"]] .= -zu

    #Main.@infiltrate
    check_nlp_sol(v, Λ, bop.n, bop.m, Γl, bop.Γ!, bop.∇ᵥF!, bop.∇ᵥΓ_size, bop.∇ᵥΓ_rows, bop.∇ᵥΓ_cols, bop.∇ᵥΓ_vals!, bop.∇²ᵥL₁_size, bop.∇²ᵥL₁_rows, bop.∇²ᵥL₁_cols, bop.∇²ᵥL₁_vals!)
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
function check_mcp_sol(n, h, z, zl, zu; tol=1e-3)
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
# check_nlp_sol(x, λ, bop.n2, bop.m2, zeros(bop.m2), bop.g!, bop.∇ₓ₂f!, bop.∇ₓ₂g_size, bop.∇ₓ₂g_rows, bop.∇ₓ₂g_cols, bop.∇ₓ₂g_vals!, bop.∇²ₓ₂L₂_size, bop.∇²ₓ₂L₂_rows, bop.∇²ₓ₂L₂_cols, bop.∇²ₓ₂L₂_vals!)

function compute_follow_feas_ind_sets(bop, v; tol=1e-3, verbosity=0)
    Js = Vector{Dict{Int,Set{Int}}}()
    h = zeros(bop.nz)
    bop.h!(h, v)
    z = @view v[bop.v_inds["z"]]
    zl = bop.zl₀
    zu = bop.zu₀
    is_valid, K = check_mcp_sol(bop.nz, h, z, zl, zu; tol)

    if !is_valid
        if verbosity > 0
            print("Failed to compute follower feasible index sets: Not a valid solution!\n")
        end
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
        hl[j] = tol
        hu[j] = Inf
        zl[j] = zl₀[j]
        zu[j] = zl₀[j]
    end
    for j in J[2] # hⱼ active (l ≠ u)
        hl[j] = 0.0
        hu[j] = 0.0
        zl[j] = zl₀[j] + tol
        zu[j] = zu₀[j]
    end
    for j in J[3] # hⱼ inactive (strictly negative), zⱼ at ub
        hl[j] = -Inf
        hu[j] = -tol
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


### unused 2025-07-10

"""
```
KKT conditions for BOPᵢ is practically an LP feasibility problem when v is given:
    ∃Λ ∈ R⁽ᵐ¹⁺ᵐʰ⁾, ∃Λ_l ∈ R⁽ⁿᵛ⁾, ∃Λ_u ∈ R⁽ⁿᵛ⁾:
            ∇ᵥF(x) - Λᵀ ∇ᵥGhs(v) = 0              (KKT conditions of BOPᵢ NLP)                               
     Gh_uᵢ ≥ Ghs(v) ≥ Gh_l₀ ⟂   Λ ≥ 0              

    ∃Λ ∈ R⁽ᵐ¹⁺ᵐʰ⁾, ∃Λ_l ∈ R⁽ⁿᵛ⁾, ∃Λ_u ∈ R⁽ⁿᵛ⁾:
        ∇ᵥF(x) - Λᵀ ∇ᵥGhs(v)  = 0                                            
     Gh_uᵢ ≥ Gh(v) ≥ Gh_l₀ ⟂    Λ ≥ 0              (KKT conditions of BOPᵢ NLP)
```

"""
function setup_check_Λ_lp_feas(nv, mΛ, Ghs!, ∇ᵥF!, ∇ᵥGhs_rows, ∇ᵥGhs_cols, ∇ᵥGhs_vals!, ∇ᵥGhs_shape, Ghs_inds, z_inds; tol=1e-6)
    check_feas = setup_lp_feas_check_HiGHS(mΛ)

    function check_Λ_lp_feas(v, z_l, z_u, h_l, h_u)
        Ghs = zeros(mΛ)
        Ghs!(Ghs, v)
        ∇ᵥF = zeros(nv)
        ∇ᵥF!(∇ᵥF, v)
        ∇ᵥGhs = sparse(∇ᵥGhs_rows, ∇ᵥGhs_cols, zeros(Cdouble, length(∇ᵥGhs_rows)), ∇ᵥGhs_shape[1], ∇ᵥGhs_shape[2])
        ∇ᵥGhs_vals!(∇ᵥGhs.nzval, v)

        # this part ensures dual feasibility
        Λ_l = fill(0.0, mΛ)
        Λ_u = fill(Inf, mΛ)
        # Λₕ=z is complement to h, s is complement to λ
        Λ_l[Ghs_inds["h"]] .= z_l
        Λ_u[Ghs_inds["h"]] .= z_u
        Λ_l[Ghs_inds["s"]] .= z_l[z_inds["s"]]
        Λ_u[Ghs_inds["s"]] .= z_u[z_inds["s"]]
        # constraint matrix is column-wise, this part checks :
        # ∇ᵥF - Λ' * ∇ᵥGhs = 0 (stationarity)
        # Λ' * Ghs = 0 (complementarity)
        A_l = [∇ᵥF; 0]
        A_u = [∇ᵥF; 0]
        A = [∇ᵥGhs'; Ghs'] # violates constraint qualifications like this

        check_feas(Λ_l, Λ_u, A_l, A_u, A)
    end
end


"""
is_minimizing = false: BOPᵢ KKT is solved as an MCP
is_minimizing = true: BOPᵢ is solved as an NLP
is_using_HSL = true: Use HSL LP solver back-end when solving the NLP
"""
#function update_vΛ!(v, Λ, Ghs_l, Ghs_u, θ_l, θ_u, θ_init, bop, Ji_bounds; solver="IPOPT", tol=1e-6, max_iter=1000)
#    is_vΛ_valid::Bool = false
#    #update_bounds!(Ghs_l, Ghs_u, θ_l, θ_u, bop, Ji_bounds)

#    v, Λ, success = solve_SBOP_nlp(bop, Ji_bounds.hl, Ji_bounds.hu, Ji_bounds.zu, Ji_bounds.zl; x_init=zeros(bop.nx), solver, tol, max_iter)

#    #if is_minimizing
#    #    v_out, Λ_out, solvestat = bop.solve_BOPᵢ_nlp(g_l=Ghs_l, g_u=Ghs_u, x_init=v, λ_init=Λ; is_using_HSL, tol)
#    #    is_vΛ_valid = solvestat == 0 #|| solvestat == 1
#    #    if is_vΛ_valid
#    #        v .= v_out # even if it's not solved, we update v so we can try to initialize z again
#    #        Λ .= Λ_out
#    #    end
#    #else
#    #    θ_init[bop.θ_inds["v"]] .= v
#    #    θ_init[bop.θ_inds["Λ"]] .= Λ[1:bop.m1+bop.mh]
#    #    θ_out, status, _ = bop.solve_BOPᵢ_KKT_mcp(x_l=θ_l, x_u=θ_u, x_init=θ_init; tol, is_silent=true, verbosity=10)
#    #    v_out = θ_out[bop.θ_inds["v"]]
#    #    Λ_out = θ_out[bop.θ_inds["Λ"]]
#    #    is_vΛ_valid = status == PATHSolver.MCP_Solved
#    #    if is_vΛ_valid
#    #        v .= v_out # even if it's not solved, we update v so we can try to initialize z again
#    #        Λ .= [Λ_out; zeros(bop.m2)]
#    #    end
#    #end
#    is_vΛ_valid
#end

function update_bounds!(Ghs_l, Ghs_u, θ_l, θ_u, bop, Ji_bounds)
    Ghs_l[bop.Ghs_inds["h"]] .= Ji_bounds.h_l
    Ghs_u[bop.Ghs_inds["h"]] .= Ji_bounds.h_u
    Ghs_l[bop.Ghs_inds["s"]] .= Ji_bounds.z_l[bop.z_inds["s"]]
    Ghs_u[bop.Ghs_inds["s"]] .= Ji_bounds.z_u[bop.z_inds["s"]]

    θ_l[bop.θ_inds["z"]] .= Ji_bounds.z_l
    θ_u[bop.θ_inds["z"]] .= Ji_bounds.z_u
    θ_l[bop.θ_inds["rh"]] .= Ji_bounds.h_l
    θ_u[bop.θ_inds["rh"]] .= Ji_bounds.h_u
end



#function check_is_sol_valid(bop, v; tol=1e-3)
#    Ghs = zeros(bop.m1 + bop.mh + bop.m2)
#    bop.Ghs!(Ghs, v)
#    h = @view Ghs[bop.Ghs_inds["h"]]
#    z = @view v[bop.v_inds["z"]]
#    zl = bop.zl₀
#    zu = bop.zu₀
#    is_valid, K = check_mcp_sol(h, z, zl, zu; tol)

#    (; is_sol_valid, K)
#    K = Dict{Int,Vector{Int}}()

#    # note which constraints are active
#    for j in 1:bop.mh
#        Kj = Int[]

#        if isapprox(z_l[j], z_u[j]; atol=2 * tol)
#            push!(Kj, 4) # case 4
#        elseif tol ≤ h[j] && z[j] < z_l[j] + tol
#            push!(Kj, 1) # case 1
#        elseif -tol < h[j] < tol && z[j] < z_l[j] + tol
#            push!(Kj, 1) # case 1 OR  
#            push!(Kj, 2) # case 2
#        elseif -tol < h[j] < tol && z_l[j] + tol ≤ z[j] ≤ z_u[j] - tol
#            push!(Kj, 2) # case 2
#        elseif -tol < h[j] < tol && z_u[j] - tol < z[j]
#            push!(Kj, 2) # case 2 OR 
#            push!(Kj, 3) # case 3
#        elseif h[j] ≤ tol && z_u[j] - tol < z[j]
#            push!(Kj, 2) # case 3
#        end
#        K[j] = Kj
#    end
#    is_sol_valid = !any(isempty.(Kj for Kj in values(K)))
#    is_sol_valid, K
#end

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
function solve_bop2(bop; x_init=zeros(bop.n1 + bop.n2), tol=1e-6, max_iter=100, verbosity=0, is_using_HSL=false, seed=0, max_init_restart_count=10, x_agree_tol=1e-3, is_using_PATH_to_init=false, conv_tol=1e-3, norm_dv_len=10, conv_dv_len=1, fol_feas_set_tol_max=1e-0, is_checking_min=true, is_checking_x_agree=false, max_invalid_sols=10)

    if is_using_HSL && !haskey(ENV, "HSL_PATH")
        is_using_HSL = false
        if verbosity > 0
            print("HSL_PATH not found: Setting is_using_HSL = false! If you would like to use HSL, please obtain a license and download HSL_jll.jl (https://licences.stfc.ac.uk/products/Software/HSL/LibHSL), and set HSL_PATH environment variable to the extracted location.\n")
        end
    end

    if bop.nθ > 300 && !haskey(ENV, "PATH_LICENSE_STRING")
        if verbosity > 0
            print("PATH_LICENSE_STRING not found and problem size is too large: Please obtain a license (https://pages.cs.wisc.edu/~ferris/path/julia/LICENSE), and set PATH_LICENSE_STRING environment variable.\n")
        end
    end

    iter_count = 0
    is_converged = false
    is_sol_valid = false
    is_dv_mon_decreasing = false
    is_norm_dv_full = false
    is_max_iter = false
    is_max_init_restarts = false
    is_max_tol_restarts = false
    is_max_invalid_sols = false
    is_prev_v_set = false
    is_none_J_feas_or_sol = false
    is_very_wrong = false

    # buffer
    x::Vector{Float64} = zeros(bop.nx)
    v = zeros(bop.n)

    Λ = zeros(bop.mΛ)
    Ghs_l = copy(bop.Ghs_l₀)
    Ghs_u = copy(bop.Ghs_u₀)
    θ_l = copy(bop.θ_l₀)
    θ_u = copy(bop.θ_u₀)
    θ_init = zeros(bop.n_θ)

    prev_iter_v = zeros(bop.nv)

    vΛ_J_inds::Vector{Int64} = [] # indexes v and Λ array wrt J
    v_arr::Vector{Vector{Float64}} = []
    Λ_arr::Vector{Vector{Float64}} = []
    vΛ_status::Vector{Int64} = [] # 0 solved, 1 feasible
    is_vΛ_feas = false
    is_vΛ_sol = false
    do_all_x_agree = false

    init_restart_count = 0
    tol_restart_count = 0

    if isnothing(seed)
        seed = 0
    end
    fol_feas_set_tol = deepcopy(tol)

    norm_dv_arr = zeros(norm_dv_len)
    chron_norm_dv_arr = copy(norm_dv_arr)
    norm_dv_cur_idx = 1
    status = 0

    is_all_J_vΛ_feas = false
    is_all_J_vΛ_sol = false
    invalid_sol_count = 0
    #prev_J_i_chosen = 0

    while !is_converged
        if iter_count >= max_iter
            if verbosity > 0
                print("Max iterations reached!\n")
            end
            is_max_iter = true
            break
        end
        iter_count += 1

        if invalid_sol_count >= max_invalid_sols
            if verbosity > 0
                print("We don't seem to be going anywhere!\n")
            end
            is_max_invalid_sols = true
            break
        end

        if verbosity > 1
            print("--Iteration $iter_count\n")
        end

        if !is_vΛ_feas
            if verbosity > 2
                print("v isn't feasible, initializing v\n")
            end
            init_z_success = initialize_z!(v, bop; is_using_HSL, verbosity, is_using_PATH=is_using_PATH_to_init, tol)

            if init_restart_count >= max_init_restart_count
                if verbosity > 0
                    print("Reached maximum init attempt, terminating!\n")
                end
                is_max_init_restarts = true
                break
            end

            if !init_z_success
                if verbosity > 1
                    print("Failed to initialize z! Randomly initializing and restarting...\n")
                end

                v[bop.v_inds["x"]] .= 10^(init_restart_count) * (0.5 .- rand(MersenneTwister(seed + init_restart_count), bop.n1 + bop.n2))
            end
            init_restart_count += 1
        end

        # By this point, v must at least satisfy the follower's problem
        if verbosity > 5
            print("Computing feasible sets for the follower at v: ")
            display(v)
        end



        #is_fol_KKT_ok = is_follower_KKT_satisfied(bop, v)
        #if !is_fol_KKT_ok
        #    if verbosity > 0
        #        print("wtf somehow follower KKT isn't satisfied\n")
        #    end
        #    continue
        #end


        follow_feas_Js = compute_follow_feas_ind_sets(bop, v; tol=fol_feas_set_tol)

        if length(follow_feas_Js) < 1
            if verbosity > 2
                print("Could not compute follower's feasible sets despite satisfying follower's KKT! Maybe it's a tolerance issue?\n")
            end

            tol_restart_count = 1

            while length(follow_feas_Js) < 1
                fol_feas_set_tol = 10^(tol_restart_count) * tol

                if fol_feas_set_tol > fol_feas_set_tol_max
                    if verbosity > 1
                        print("Reached max follow set tol relaxations!\n")
                    end
                    fol_feas_set_tol = fol_feas_set_tol_max
                    is_max_tol_restarts = true
                    break
                end

                if verbosity > 1
                    print("Relaxing follower feasible set tolerance $(round(fol_feas_set_tol,sigdigits=2)) and trying again...\n")
                end
                follow_feas_Js = compute_follow_feas_ind_sets(bop, v; tol=fol_feas_set_tol)
                tol_restart_count += 1
            end
        end

        fol_feas_set_tol = deepcopy(tol)
        n_J = length(follow_feas_Js)

        if verbosity > 1
            if n_J > 1
                print("Multiple feasible sets ($n_J) detected!\n")
            end
        end


        empty!(vΛ_J_inds)
        empty!(v_arr)
        empty!(Λ_arr)
        empty!(vΛ_status)

        # check if there exists Λ for all i: x ∈ Hᵢ
        # this step could be done in parallel
        for i in 1:n_J
            Ji = follow_feas_Js[i]
            Ji_bounds = convert_J_to_h_z_bounds(Ji, bop)

            ## does there exist Λ for v? feasibility problem so is_minimizing=false
            is_vΛ_feas = update_vΛ!(v, Λ, Ghs_l, Ghs_u, θ_l, θ_u, θ_init, bop, Ji_bounds; is_minimizing=false, is_using_HSL, tol)
            # debug
            #is_Λ_feas_2 = bop.check_Λ_lp_feas(v, Ji_bounds.z_l, Ji_bounds.z_u, Ji_bounds.h_l, Ji_bounds.h_u)

            ##check_v_Λ_BOPᵢ_KKT(v, Λ, bop, Ghs_l, Ghs_u)
            #####
            #if (is_vΛ_feas && !is_Λ_feas_2) || (!is_vΛ_feas && is_Λ_feas_2)
            #    #TODO 2025-06-29 jesus 
            #    @info "update_v! returns $is_vΛ_feas but bop.check_Λ_lp_feas returns $is_Λ_feas_2"
            #    
            #end

            if is_vΛ_feas
                if is_follower_minimized(bop, v, Ji)
                    push!(vΛ_J_inds, i)
                    push!(v_arr, copy(v))
                    push!(Λ_arr, copy(Λ))
                    push!(vΛ_status, 1) # feasible
                    if verbosity > 2
                        print("i=$i: Feasible J1/2/3/4: $(Ji[1])/$(Ji[2])/$(Ji[3])/$(Ji[4])\n")
                    end

                end
            else
                if verbosity > 1
                    print("i=$i: INFEASIBLE J1/2/3/4: $(Ji[1])/$(Ji[2])/$(Ji[3])/$(Ji[4])\n")
                end
            end
        end

        if length(vΛ_status) == n_J && n_J > 0
            is_all_J_vΛ_feas = true
        else
            is_all_J_vΛ_feas = false
        end

        if is_all_J_vΛ_feas
            if verbosity > 4
                print("All $n_J follower solutions have feasible v and Λ.\n")
            end
        end

        if is_checking_min
            for i in 1:n_J
                Ji = follow_feas_Js[i]
                Ji_bounds = convert_J_to_h_z_bounds(Ji, bop)
                # if there exists feasible v, Λ, we initialize with it
                vΛ_i = findfirst(vΛ_J_inds .== i)
                if !isnothing(vΛ_i) && vΛ_status[vΛ_i] == 1
                    v .= v_arr[vΛ_i]
                    Λ .= Λ_arr[vΛ_i]
                end

                # is there a Λ for v that solves the problem?
                is_vΛ_sol = update_vΛ!(v, Λ, Ghs_l, Ghs_u, θ_l, θ_u, θ_init, bop, Ji_bounds; is_minimizing=true, is_using_HSL, tol)


                # debug
                #is_Λ_feas_2 = bop.check_Λ_lp_feas(v, Ji_bounds.z_l, Ji_bounds.z_u, Ji_bounds.h_l, Ji_bounds.h_u)

                ####
                #if (is_vΛ_sol && !is_Λ_feas_2) || (!is_vΛ_sol && is_Λ_feas_2)
                #    #TODO 2025-06-29 jesus 
                #    @info "update_v!min returns $is_vΛ_sol but bop.check_Λ_lp_feas returns $is_Λ_feas_2"
                #    
                #end

                if is_vΛ_sol
                    #if !check_v_Λ_BOPᵢ_KKT(v, Λ, bop, Ghs_l, Ghs_u)
                    #    @info "woah there"
                    #end

                    is_vΛ_feas = is_vΛ_sol
                    if !isnothing(vΛ_i)
                        v_arr[vΛ_i] .= v
                        Λ_arr[vΛ_i] .= Λ
                        vΛ_status[vΛ_i] = 0
                    else
                        push!(vΛ_J_inds, i)
                        push!(v_arr, copy(v))
                        push!(Λ_arr, copy(Λ))
                        push!(vΛ_status, 0)
                    end

                    if verbosity > 2
                        print("i=$i: Solved J1/2/3/4: $(Ji[1])/$(Ji[2])/$(Ji[3])/$(Ji[4])\n")
                    end
                else
                    if verbosity > 1
                        print("i=$i: UNSOLVED J1/2/3/4: $(Ji[1])/$(Ji[2])/$(Ji[3])/$(Ji[4])\n")
                    end
                    if !is_checking_x_agree
                        break # no need to solve further
                    end
                end
            end

            is_all_J_vΛ_sol = false
            if length(vΛ_status) == n_J && n_J > 0
                if all(vΛ_status .== 0)
                    is_all_J_vΛ_sol = true
                end
            end

            if is_all_J_vΛ_sol
                if verbosity > 4
                    print("All $n_J follower solutions solved for v and Λ.\n")
                end
            end

            if is_checking_x_agree && is_all_J_vΛ_sol
                do_all_x_agree = true

                for (i, vv) in enumerate(v_arr)
                    if i < n_J
                        norm_x_err = LinearAlgebra.norm(v_arr[i+1][bop.v_inds["x"]] - vv[bop.v_inds["x"]]) # only checking x error
                        if (norm_x_err > x_agree_tol)
                            if verbosity > 1
                                print("BOPᵢ solutions disagree! norm x err: $norm_x_err\n")
                            end
                            do_all_x_agree = false
                            break
                        end
                    end
                end
            end
        end

        if length(vΛ_status) == 0
            is_none_J_feas_or_sol = true
            is_sol_valid = false
            if verbosity > 1
                print("Awkward!!\n")
            end
            break
        end

        # this makes sure next v and Λ are at least feasible, AiyoshiShimizu1984Ex2 happens to converge to optimal with this
        #v .= v_arr[end]
        #Λ .= Λ_arr[end]

        # choose the point the one that leads to lowest cost, but this makes or breaks so idk what to do
        is_there_one_sol = any(vΛ_status .== 0)
        F = Inf
        #n_Js = map(v -> length(compute_follow_feas_ind_sets(bop, v; tol=fol_feas_set_tol)), v_arr)
        #Fs = map(v -> bop.F(v[bop.v_inds["x"]]), v_arr)
        #fs = map(v -> bop.f(v[bop.v_inds["x"]]), v_arr)
        #Gs = mapreduce(v -> bop.G(v[bop.v_inds["x"]])', vcat, v_arr)
        #gs = mapreduce(v -> bop.g(v[bop.v_inds["x"]])', vcat, v_arr)
        #max_n_J = maximum(n_Js)

        #fs_sorted_inds = sortperm(fs)

        #@info "$n_Js $Fs $fs"



        for (i, stat) in enumerate(vΛ_status)
            #if prev_J_i_chosen == i && n_J > 1
            #    continue
            #end
            if is_there_one_sol && is_checking_min && stat != 0# if checking min ignore the feasible points if there's at least one sol, otherwise choose a feasible point
                continue
            end
            # maintain larger sols
            #if n_Js[i] != max_n_J
            #   continue
            #end

            #if !all(Gs[i, :] .≥ 0 - tol) || !all(gs[i, :] .≥ 0 - tol)
            #    @info "this one isn't it chief"
            #    continue
            #end

            # choose the smallest available f value
            F_ = bop.F(v_arr[i][bop.v_inds["x"]])

            #n_J_ = length(compute_follow_feas_ind_sets(bop, v_arr[i]; tol=fol_feas_set_tol))

            if F_ < F
                F = F_
                v .= v_arr[i]
                Λ .= Λ_arr[i]
                #prev_J_i_chosen = i
            end
            #break
        end


        is_sol_valid = false
        if is_checking_min
            if is_all_J_vΛ_sol # there's at least a feas solution
                if is_checking_x_agree
                    if do_all_x_agree
                        is_sol_valid = true
                    end
                else
                    is_sol_valid = true
                end
            end
        elseif is_all_J_vΛ_feas
            is_sol_valid = true
            # this is a good faith assumption if we don't check for solutions
        end


        if !is_sol_valid
            invalid_sol_count += 1
            continue
        else
            invalid_sol_count = 0
        end

        # by this point is_sol_valid = true: v's are sols and agree, Λs are feasible, all we have to check now is if the solution has converged
        if !is_prev_v_set
            prev_iter_v .= v
            is_prev_v_set = true
        else
            dv = v[bop.v_inds["x"]] - prev_iter_v[bop.v_inds["x"]]
            prev_iter_v .= v
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
                break
            elseif is_norm_dv_full && all(diff(chron_norm_dv_arr) .≤ -tol)
                # also it's worth checking if dv is slowly
                if verbosity > 3
                    print("last $norm_dv_len norm(dv) is monotonously decreasing without meeting conv tol, terminating\n")
                end
                is_dv_mon_decreasing = true
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

        if verbosity > 5
            print("v = ")
            display(v)
        end
    end

    if is_converged
        if verbosity > 1
            print("\n")
        end
    end

    x .= v[bop.v_inds["x"]]

    # final sanity check
    if is_sol_valid && !(all(bop.G(x) .≥ 0 - tol) && all(bop.g(x) .≥ 0 - tol))# && is_follower_KKT_satisfied(bop, v))
        if verbosity > 0
            print("Something went VERY wrong, this allegedly valid solution is bilevel infeasible!\n")
        end
        is_very_wrong = true
        is_sol_valid = false
    end

    status = get_status(is_sol_valid, is_converged, is_dv_mon_decreasing, is_max_init_restarts, is_max_tol_restarts, is_max_iter, is_none_J_feas_or_sol, is_max_invalid_sols, is_very_wrong)

    (; x, status, iter_count)
end


"""
```
Status:
    -7: invalid v: other
    -6: invalid v: something went very wrong
    -5: invalid v: max invalid sols
    -4: invalid v: none J is neither feasible or solvable
    -3: invalid v: max iterations
    -2: invalid v: max init restarts
    -1: invalid v: max follower relaxations
    0: success: valid v, converged
    1: valid v: dv is monoton decreasing
    2: valid v: max iterations
    3: valid v: other
```
"""
function get_status(is_sol_valid, is_converged, is_dv_mon_decreasing, is_max_init_restarts, is_max_fol_relaxes, is_max_iter, is_none_J_feas_or_sol, is_max_invalid_sols, is_very_wrong)

    if !is_sol_valid
        if is_max_fol_relaxes
            return -1
        elseif is_max_init_restarts
            return -2
        elseif is_max_iter
            return -3
        elseif is_none_J_feas_or_sol
            return -4
        elseif is_max_invalid_sols
            return -5
        elseif is_very_wrong
            return -6
        else
            return -7
        end
    end

    if is_converged
        return 0
    elseif is_dv_mon_decreasing
        return 1
    elseif is_max_iter
        return 2
    else
        return 3
    end

end