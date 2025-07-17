using BilevelOptSolver
include("../src/BOLIB_utils.jl")
using Symbolics
using SparseArrays
using BenchmarkTools
using ProfileView

function solve_bopper()
    #b = BOLIB.DempeDutta2012Ex24()
	b = DempeDutta2012Ex34()
    bop, syms = construct_bop(b.n1, b.n2, b.F, b.G, b.f, b.g, verbosity=0)

    x1 = b.xy_init[bop.inds.x["x1"]]
    #x1 = [0; 0]
    #@info "solve follower at x1=$x1"
    tol = 1e-7

    #(x, λ, success) = BilevelOptSolver.solve_follower_nlp(bop, x1; x2_init=b.xy_init[bop.inds.x["x2"]], solver="IPOPT", tol)
    #v = [x; λ]
    #is_fol_nec, is_fol_suf = BilevelOptSolver.check_follower_sol(v, bop)
    #follow_feas_Js = BilevelOptSolver.compute_follow_feas_ind_sets(bop, v)

    solve_SBOPi_IPOPT! = BilevelOptSolver.setup_solve_SBOPi_IPOPT(bop)
    Λ = zeros(bop.m)

    investigated_J2 = []
    is_v_fol_nec::Vector{Bool} = []
    is_v_fol_suf::Vector{Bool} = []
    v_arr::Vector{Vector{Float64}} = []
    Λ_arr::Vector{Vector{Float64}} = []
    Ff_arr = []

    function update_v_recursive(Js, investigated_J2, is_v_fol_nec, is_v_fol_suf, v_arr, Λ_arr, Ff_arr)
        for J in Js
            if J[2] in investigated_J2
                #@info "we've seen $(J[2])"
                continue
            end
            hl, hu, zl, zu = BilevelOptSolver.convert_J_to_h_z_bounds(J, bop)
            is_solved = solve_SBOPi_IPOPT!(v, Λ, hl, hu, zl, zu; tol)

            is_fol_nec, is_fol_suf = BilevelOptSolver.check_follower_sol(v, bop; tol=1e1 * tol)
            is_SBOPi_nec, is_SBOPi_suf = BilevelOptSolver.check_SBOPi_sol(v, Λ, bop, hl, hu, zu, zl; tol=1e1 * tol) # is_SBOPi_suf 
            print("$(v[bop.inds.v["x"]]) Solved nc: $is_SBOPi_nec sf: $is_SBOPi_suf (follower nc: $is_fol_nec sc: $is_fol_suf) J1: $(J[1]) J2: $(J[2])\n")
            push!(investigated_J2, J[2])
            push!(is_v_fol_nec, is_fol_nec)
            push!(is_v_fol_suf, is_fol_suf)
            push!(v_arr, copy(v))
            push!(Λ_arr, copy(Λ))
            Ff = [bop.F(v[bop.inds.v["x"]]); bop.f(v[bop.inds.v["x"]])]
            push!(Ff_arr, Ff)

            if is_solved
                follow_feas_Js = BilevelOptSolver.compute_follow_feas_ind_sets(bop, v; tol)
                update_v_recursive(follow_feas_Js, investigated_J2, is_v_fol_nec, is_v_fol_suf, v_arr, Λ_arr, Ff_arr)
            end
        end
    end

    x1 = b.xy_init[bop.inds.x["x1"]]
    (x, λ, success) = BilevelOptSolver.solve_follower_nlp(bop, x1; solver="IPOPT")

    if success
        v = [x; λ]
        follow_feas_Js = BilevelOptSolver.compute_follow_feas_ind_sets(bop, v)
        update_v_recursive(follow_feas_Js, investigated_J2, is_v_fol_nec, is_v_fol_suf, v_arr, Λ_arr, Ff_arr)

        @info "checking..."
        F_arr = [F[1] for F in Ff_arr]
        #F_arr[is_v_fol_suf .== false] .= Inf #checking min equals to this
        #Main.@infiltrate
        x_init = v_arr[argmin(F_arr)][bop.inds.v["x"]]
        #v .= v_arr[argmin(F_arr)]
        #Λ .= Λ_arr[argmin(F_arr)]

        ##v = v_arr[argmin(F_arr)]
        ##Λ = Λ_arr[argmin(F_arr)]
        #follow_feas_Js = BilevelOptSolver.compute_follow_feas_ind_sets(bop, v)

        #for J in follow_feas_Js
        #	hl, hu, zl, zu = BilevelOptSolver.convert_J_to_h_z_bounds(J, bop)
        #	is_solved = solve_SBOPi_IPOPT!(v, Λ, hl, hu, zl, zu)

        #	is_fol_nec, is_fol_suf = BilevelOptSolver.check_follower_sol(v, bop; tol=1e1 * tol)
        #	is_SBOPi_nec, is_SBOPi_suf = BilevelOptSolver.check_SBOPi_sol(v, Λ, bop, hl, hu, zu, zl; tol=1e1 * tol) # is_SBOPi_suf 
        #	print("$(v[bop.inds.v["x"]]) Solved nc: $is_SBOPi_nec sf: $is_SBOPi_suf (follower nc: $is_fol_nec sc: $is_fol_suf) J1: $(J[1]) J2: $(J[2])\n")
        #end
        elapsed_time = @elapsed begin
            is_sol_valid, x, λ, iter_count, status = solve_bop(bop; max_iter=50, x_init, verbosity=5, tol, norm_dv_len=10, conv_dv_len=1, is_checking_min=false, is_checking_x_agree=false, init_solver="IPOPT", solver="IPOPT")
        end

        is_optimal, is_best, Ff, Ff_star, rating = rate_BOLIB_result(b, bop, x)
        if is_sol_valid
            print("success ")
        else
            print("FAIL ")
        end
        print("($status), $iter_count iters ($(round(elapsed_time, sigdigits=5)) s),\t" * rating * ",\t x = $(round.(x, sigdigits=5)) -> Ff = $(round.(Ff, sigdigits=5)) (Ff* = $(round.(Ff_star, sigdigits=5)))\n")
    end
end

solve_bopper()



#J1 = follow_feas_Js[1]
#J2 = follow_feas_Js[2]
#J3 = follow_feas_Js[3]
#J4 = follow_feas_Js[4]




#solve_SBOPi_PATH! = BilevelOptSolver.setup_solve_SBOPi_PATH(bop)
#is_solved = solve_SBOPi_PATH!(v, Λ, hl, hu, zl, zu)

#is_necessary, is_sufficient = BilevelOptSolver.check_nlp_sol(x, λ, bop.n2, bop.m2, zeros(bop.m2), bop.g!, bop.∇ₓ₂f!, bop.∇ₓ₂g_size, bop.∇ₓ₂g_rows, bop.∇ₓ₂g_cols, bop.∇ₓ₂g_vals!, bop.∇²ₓ₂L2_size, bop.∇²ₓ₂L2_rows, bop.∇²ₓ₂L2_cols, bop.∇²ₓ₂L2_vals!)

#∇²ₓ₂L2 = sparse(bop.∇²ₓ₂L2_rows, bop.∇²ₓ₂L2_cols, zeros(length(bop.∇²ₓ₂L2_rows)), bop.∇²ₓ₂L2_size[1], bop.∇²ₓ₂L2_size[2])
#bop.∇²ₓ₂L2_vals!(∇²ₓ₂L2.nzval, x, λ, 1.0)

#substitute(syms.f, Dict([syms.x[1] => x[1]]))

#include("../src/forrest_solver.jl")
#using .forrest_solver
#OP1 = forrest_solver.OptimizationProblem(b.n1 + b.n2, 1:b.n1, b.F, b.G, zeros(bop.m1), fill(Inf, bop.m1))
#OP2 = forrest_solver.OptimizationProblem(b.n1 + b.n2, 1:b.n2, b.f, b.g, zeros(bop.m2), fill(Inf, bop.m2))
#bilevel = [OP1; OP2]
#out_bilevel = forrest_solver.solve(bilevel)
#x_forrest = out[1:b.n1+b.n2]
#@info x_forrest 
#@info bop.F(x_forrest)

#nash = [OP1 OP2]
#out_nash = forrest_solver.solve(nash)