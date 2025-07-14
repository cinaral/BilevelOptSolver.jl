

#function solve_high_point_nlp(bop; x_init=zeros(bop.nx), tol=1e-6, max_iter=1000)
#    xl = fill(-Inf, bop.n1 + bop.n2)
#    xu = fill(Inf, bop.n1 + bop.n2)
#    Gg_l = fill(0.0, bop.m1 + bop.m2)
#    Gg_u = fill(Inf, bop.m1 + bop.m2)

#    solve_high_point_nlp = setup_nlp_solve_IPOPT(bop.nx, bop.m1 + bop.m2, xl, xu, Gg_l, Gg_u, bop.F, bop.Gg!, bop.∇ₓF!, bop.∇ₓGg_rows, bop.∇ₓGg_cols, bop.∇ₓGg_vals!, bop.∇²ₓL3_rows, bop.∇²ₓL3_cols, bop.∇²ₓL3_vals!)

#    x, λ, solvestat, _ = solve_high_point_nlp(; x_init, tol, max_iter, print_level=0)

#    success = solvestat == 0 || solvestat == 1 # this is acceptable too since we only care about feasibility
#    (; x, λ, success)
#end

#function solve_follower_nlp(bop, x1; x2_init=zeros(bop.n2), solver="IPOPT", tol=1e-6, max_iter=1000)
#    x = zeros(bop.nx)
#    x[bop.inds.x["x1"]] .= x1
#    λ = zeros(bop.m2)
#    success = false

#    if solver == "IPOPT"
#        x2_l = fill(-Inf, bop.n2)
#        x2_u = fill(Inf, bop.n2)
#        gl = fill(0.0, bop.m2)
#        gu = fill(Inf, bop.m2)

#        function eval_f(x2::Vector{Float64})
#            x[bop.inds.x["x2"]] .= x2
#            bop.f(x)
#        end
#        function eval_g(g::Vector{Float64}, x2::Vector{Float64})
#            x[bop.inds.x["x2"]] .= x2
#            bop.g!(g, x)
#        end
#        function eval_∇ₓ₂f(∇ₓ₂f::Vector{Float64}, x2::Vector{Float64})
#            x[bop.inds.x["x2"]] .= x2
#            bop.∇ₓ₂f!(∇ₓ₂f, x)
#        end
#        function eval_∇ₓ₂g_vals(∇ₓ₂g_vals::Vector{Float64}, x2::Vector{Float64})
#            x[bop.inds.x["x2"]] .= x2
#            bop.∇ₓ₂g_vals!(∇ₓ₂g_vals, x)
#        end
#        function eval_∇²ₓ₂L2_vals(∇²ₓ₂L2_vals::Vector{Float64}, x2::Vector{Float64}, λ::Vector{Float64}, obj_factor::Float64)
#            x[bop.inds.x["x2"]] .= x2
#            bop.∇²ₓ₂L2_vals!(∇²ₓ₂L2_vals, x, λ, obj_factor)
#        end
#        solve = setup_nlp_solve_IPOPT(bop.n2, bop.m2, x2_l, x2_u, gl, gu, eval_f, eval_g, eval_∇ₓ₂f, bop.∇ₓ₂g_rows, bop.∇ₓ₂g_cols, eval_∇ₓ₂g_vals, bop.∇²ₓ₂L2_rows, bop.∇²ₓ₂L2_cols, eval_∇²ₓ₂L2_vals)

#        x_out, λ_out, solvestat, _ = solve(; x_init=x2_init, tol, max_iter, print_level=0)

#        x[bop.inds.x["x2"]] .= x_out
#        λ .= λ_out
#        success = solvestat == 0 # || solvestat == 1

#    elseif solver == "PATH"
#        v = zeros(bop.n)
#        v[bop.inds.v["x"]] .= x
#        z_init = zeros(bop.nz)
#        z_init[bop.inds.z["x2"]] = x2_init

#        function eval_F!(h, z::Vector{Float64})
#            v[bop.inds.v["z"]] .= z
#            bop.h!(h, v, 1.0)
#        end
#        function eval_J_vals!(∇h, z::Vector{Float64})
#            v[bop.inds.v["z"]] .= z
#            bop.∇h_vals!(∇h, v, 1.0)
#        end
#        solve = setup_mcp_solve_PATH(bop.nz, bop.zl₀, bop.zu₀, eval_F!, bop.∇h_rows, bop.∇h_cols, eval_J_vals!)

#        z_out, status, _ = solve(; x_init=z_init, tol, max_iter, is_silent=true)

#        x[bop.inds.x["x2"]] .= z_out[bop.inds.z["x2"]]
#        λ .= z_out[bop.inds.z["λ"]]
#        success = status == PATHSolver.MCP_Solved
#    end

#    (; x, λ, success)
#end

# 2025-07-10 in order to reduce allocations we could do...
#function setup_SBOP_nlp_IPOPT()
#end

### unused 2025-07-10

#function setup_SBOP_nlp_PATH()
#end

function setup_solve_SBOPi_nlp!(bop; solver="IPOPT")
    if solver == "IPOPT"
        vl = fill(-Inf, bop.n)
        vu = fill(Inf, bop.n)
        Γl = fill(0.0, bop.m)
        Γu₀ = fill(Inf, bop.m) # doesn't change

        solve = setup_nlp_solve_IPOPT(bop.n, bop.m, vl, vu, Γl, Γu₀, bop.Fv, bop.Γ!, bop.∇ᵥF!, bop.∇ᵥΓ_rows, bop.∇ᵥΓ_cols, bop.∇ᵥΓ_vals!, bop.∇²ᵥL1_rows, bop.∇²ᵥL1_cols, bop.∇²ᵥL1_vals!)
    elseif solver == "PATH"
        θ_init = zeros(bop.nθ)
        #θ_init[bop.inds.θ["v"]] .= v
        #θ_init[bop.inds.θ["Λ"]] .= Λ
        θl = copy(bop.θl₀)
        #θl[bop.inds.θ["Λhl"]] .= zl
        #θl[bop.inds.θ["Λzl"]] .= hl
        θu = copy(bop.θu₀)
        #θu[bop.inds.θ["Λhu"]] .= zu
        #θu[bop.inds.θ["Λzu"]] .= hu

        solve = setup_mcp_solve_PATH(bop.nθ, θl, θu, bop.Φ!, bop.∇Φ_rows, bop.∇Φ_cols, bop.∇Φ_vals!)
    end

    (; solve, Γl, θ_init, θl, θu)
end


function update_bounds!(Γl, θ_init, θl, θu, bop, v, Λ, hl, hu, zl, zu)
    Γl[bop.inds.Γ["hl"]] .= hl
    Γl[bop.inds.Γ["hu"]] .= -hu
    Γl[bop.inds.Γ["zl"]] .= zl
    Γl[bop.inds.Γ["zu"]] .= -zu

    θ_init[bop.inds.θ["v"]] .= v
    θ_init[bop.inds.θ["Λ"]] .= Λ
    θl[bop.inds.θ["Λhl"]] .= zl
    θl[bop.inds.θ["Λzl"]] .= hl
    θu[bop.inds.θ["Λhu"]] .= zu
    θu[bop.inds.θ["Λzu"]] .= hu
end



function solve_SBOPi_nlp!(v, Λ, Γl, θ_init, θl, θu, solve; solver="IPOPT", tol=1e-6, max_iter=1000)
    #v = zeros(bop.n)
    #v[bop.inds.v["x"]] .= x_init
    #Λ = zeros(bop.m)
    success = false

    if solver == "IPOPT"
        ##vl = fill(-Inf, bop.n)
        ##vu = fill(Inf, bop.n)
        ##Γl = fill(0.0, bop.m)
        #Γl[bop.inds.Γ["hl"]] .= hl
        #Γl[bop.inds.Γ["hu"]] .= -hu
        #Γl[bop.inds.Γ["zl"]] .= zl
        #Γl[bop.inds.Γ["zu"]] .= -zu
        ##Γu₀ = fill(Inf, bop.m) # doesn't change

        ##solve = setup_nlp_solve_IPOPT(bop.n, bop.m, vl, vu, Γl, Γu₀, bop.Fv, bop.Γ!, bop.∇ᵥF!, bop.∇ᵥΓ_rows, bop.∇ᵥΓ_cols, bop.∇ᵥΓ_vals!, bop.∇²ᵥL1_rows, bop.∇²ᵥL1_cols, bop.∇²ᵥL1_vals!)

        v_out, Λ_out, solvestat, _ = solve(; gl=Γl, x_init=v, λ_init=Λ, tol, max_iter, print_level=0)
        success = solvestat == 0 # || solvestat == 1

    elseif solver == "PATH"
        ##θ_init = zeros(bop.nθ)
        #θ_init[bop.inds.θ["v"]] .= v
        #θ_init[bop.inds.θ["Λ"]] .= Λ
        ##θl = copy(bop.θl₀)
        #θl[bop.inds.θ["Λhl"]] .= zl
        #θl[bop.inds.θ["Λzl"]] .= hl
        ##θu = copy(bop.θu₀)
        #θu[bop.inds.θ["Λhu"]] .= zu
        #θu[bop.inds.θ["Λzu"]] .= hu

        #solve = setup_mcp_solve_PATH(bop.nθ, θl, θu, bop.Φ!, bop.∇Φ_rows, bop.∇Φ_cols, bop.∇Φ_vals!)
        #Main.@infiltrate
        θ_out, status, _ = solve(; xl=θl, xu=θu, x_init=[θ_init], tol, is_silent=true)

        v_out = θ_out[bop.inds.θ["v"]]
        Λ_out = θ_out[bop.inds.θ["Λ"]]
        success = status == PATHSolver.MCP_Solved
    end

    if success
        v .= v_out
        Λ .= Λ_out
    end
    success
end



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
#    #    θ_init[bop.inds.θ["v"]] .= v
#    #    θ_init[bop.inds.θ["Λ"]] .= Λ[1:bop.m1+bop.mh]
#    #    θ_out, status, _ = bop.solve_BOPᵢ_KKT_mcp(x_l=θ_l, x_u=θ_u, x_init=θ_init; tol, is_silent=true, verbosity=10)
#    #    v_out = θ_out[bop.inds.θ["v"]]
#    #    Λ_out = θ_out[bop.inds.θ["Λ"]]
#    #    is_vΛ_valid = status == PATHSolver.MCP_Solved
#    #    if is_vΛ_valid
#    #        v .= v_out # even if it's not solved, we update v so we can try to initialize z again
#    #        Λ .= [Λ_out; zeros(bop.m2)]
#    #    end
#    #end
#    is_vΛ_valid
#end


#function check_is_sol_valid(bop, v; tol=1e-3)
#    Ghs = zeros(bop.m1 + bop.mh + bop.m2)
#    bop.Ghs!(Ghs, v)
#    h = @view Ghs[bop.Ghs_inds["h"]]
#    z = @view v[bop.inds.v["z"]]
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

                v[bop.inds.v["x"]] .= 10^(init_restart_count) * (0.5 .- rand(MersenneTwister(seed + init_restart_count), bop.n1 + bop.n2))
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
                        norm_x_err = LinearAlgebra.norm(v_arr[i+1][bop.inds.v["x"]] - vv[bop.inds.v["x"]]) # only checking x error
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
        #Fs = map(v -> bop.F(v[bop.inds.v["x"]]), v_arr)
        #fs = map(v -> bop.f(v[bop.inds.v["x"]]), v_arr)
        #Gs = mapreduce(v -> bop.G(v[bop.inds.v["x"]])', vcat, v_arr)
        #gs = mapreduce(v -> bop.g(v[bop.inds.v["x"]])', vcat, v_arr)
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
            F_ = bop.F(v_arr[i][bop.inds.v["x"]])

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
            dv = v[bop.inds.v["x"]] - prev_iter_v[bop.inds.v["x"]]
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

    x .= v[bop.inds.v["x"]]

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