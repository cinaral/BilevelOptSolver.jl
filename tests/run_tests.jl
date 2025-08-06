using BilevelOptSolver

# 2025-08-06 TODO: add randomization, add forced randomization, fix x agreement logic

include("./lp-qp/mb_2006_01.jl")
include("./lp-qp/mb_2007_03.jl")
include("./lp-qp/mb_2007_04.jl")
include("./lp-qp/b_1991_02.jl")
include("./lp-qp/as_1984_01.jl")

include("./lp-nlp/mb_2007_05.jl")
include("./lp-nlp/mb_2007_06.jl")
include("./lp-nlp/mb_2007_13.jl")
include("./lp-nlp/mb_2007_13v.jl")


#probs = [mb_2006_01(); mb_2007_04()]
probs = [mb_2007_13()];

function run_tests()

    for p in probs
		print("$(p.name):\n")
        bop, syms = construct_bop(p.n1, p.n2, p.F, p.G, p.f, p.g; verbosity=0, np=0)

        elapsed_time = @elapsed begin
            is_sol_valid, x, λ, iter_count, status = solve_bop(bop; max_iter=50, x_init=p.xy_init, verbosity=5, tol=1e-7, norm_dv_len=10, conv_dv_len=1, is_checking_min=true, is_checking_x_agree=true, is_always_hp=true, init_solver="IPOPT", solver="IPOPT")
        end

        Ff = [bop.F(x); bop.f(x)]
        if is_sol_valid
            print("success ")
        else
            print("FAIL ")
        end
        print("($status), $iter_count iters,\t", "x = $(round.(x, sigdigits=5)) -> Ff = $(round.(Ff, sigdigits=5)) (Ff* = $(round.(b.Ff_optimal[1:2], sigdigits=5)))\n\n")
    end
end

run_tests()