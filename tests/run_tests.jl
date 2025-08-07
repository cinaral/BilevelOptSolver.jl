using BilevelOptSolver
using Test
using Random
# 2025-08-06 TODO: add randomization, add forced randomization, fix x agreement logic

function call_solve_bop(p, x_init; verbosity=0)
    if verbosity > 0
        print("$(p.name):\n")
    end

    bop, _ = construct_bop(p.n1, p.n2, p.F, p.G, p.f, p.g; verbosity=0, np=0)
    is_sol_valid, x, λ, iter_count, status = solve_bop(bop; max_iter=50, x_init, verbosity, tol=1e-7, norm_dv_len=10, conv_dv_len=1, is_checking_min=true, is_checking_x_agree=true, is_always_hp=false, init_solver="IPOPT", solver="IPOPT")

    Ff = [bop.F(x); bop.f(x)]
    if verbosity > 0
        if is_sol_valid
            print("success ")
        else
            print("FAIL ")
        end
        print("($status), $iter_count iters,\t", "x = $(round.(x, sigdigits=5)) -> Ff = $(round.(Ff, sigdigits=5)) (Ff* = $(round.(p.Ff_optimal[1:2], sigdigits=5)))\n\n")
    end

    (; is_sol_valid, x, λ, iter_count, status, Ff)
end

include("./lp-lp/mb_2007_01.jl")
"""
all should solve
"""
function test_mb_2007_01(; tol=1e-6, n=10)
    rng = MersenneTwister(123)
    p = mb_2007_01()
    @testset "mb_2007_01" begin
        for _ = 1:n
            x_init = randn(rng, 1, 1)
            is_sol_valid, x, _, _, _ = call_solve_bop(p, x_init)
            @test is_sol_valid && isapprox(x, p.xy_optimal; atol=tol)
        end
    end
end

include("./lp-lp/mb_2007_02.jl")
"""
no solution, none should solve
"""
function test_mb_2007_02(; tol=1e-6, n=10)
    rng = MersenneTwister(123)
    p = mb_2007_02()
    @testset "mb_2007_02" begin
        for _ = 1:n
            x_init = randn(rng, 1, 1)
            is_sol_valid, _, _, _, _ = call_solve_bop(p, x_init)
            @test !is_sol_valid
        end
    end
end

include("./lp-nlp/mb_2007_05.jl")
"""
"""
function test_mb_2007_05(; tol=1e-6, n=10)
    p = mb_2007_05();
    @testset "mb_2007_05" begin
        # expected to converge to the wrong local min
        x_init = [-0.1]
        is_sol_valid, x, _, _, _ = call_solve_bop(p, x_init)
        @test is_sol_valid && isapprox(x, [-0.5]; atol=tol)

        x_init = [-0.09375]
        is_sol_valid, x, _, _, _ = call_solve_bop(p, x_init)
        @test is_sol_valid && isapprox(x, p.xy_optimal; atol=tol)

        x_init = [0.1]
        is_sol_valid, x, _, _, _ = call_solve_bop(p, x_init)
        @test is_sol_valid && isapprox(x, p.xy_optimal; atol=tol)
    end
end

include("./lp-nlp/mb_2007_06.jl")
"""
"""
function test_mb_2007_06(; tol=1e-6, n=10)
    p = mb_2007_06();
    @testset "mb_2007_06" begin
        x_init = [-0.1]
        is_sol_valid, x, _, _, _ = call_solve_bop(p, x_init)
        @test is_sol_valid && isapprox(x, p.xy_optimal; atol=tol)

        # needs randomization
        x_init = [0.0]
        is_sol_valid, x, _, _, _ = call_solve_bop(p, x_init)
        @test is_sol_valid && isapprox(x, p.xy_optimal; atol=tol)

        x_init = [0.1]
        is_sol_valid, x, _, _, _ = call_solve_bop(p, x_init)
        @test is_sol_valid && isapprox(x, p.xy_optimal; atol=tol)
    end
end


include("./lp-nlp/mb_2007_13.jl")


include("./lp-nlp/mb_2007_13v.jl")

include("./lp-qp/mb_2006_01.jl")
include("./lp-qp/mb_2007_03.jl")
include("./lp-qp/mb_2007_04.jl")
include("./lp-qp/b_1991_02.jl")
include("./lp-qp/as_1984_01.jl")

include("./qp-qp/as_1981_01.jl")


function run_tests()
    test_mb_2007_01()
    test_mb_2007_02()
    test_mb_2007_05()
    test_mb_2007_06()
end

run_tests()