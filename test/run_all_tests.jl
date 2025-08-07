using BilevelOptSolver
using Test
using Random
# 2025-08-06 TODO: add randomization, add forced randomization, fix x agreement logic

function call_solve_bop(p, x_init; verbosity=0, tol=1e-7)
    if verbosity > 0
        print("$(p.name):\n")
    end

    bop, syms = construct_bop(p.n1, p.n2, p.F, p.G, p.f, p.g; verbosity=0, np=0)
    is_sol_valid, x, λ, iter_count, status = solve_bop(bop; max_iter=50, x_init, verbosity, tol, norm_dv_len=10, conv_dv_len=1, is_checking_min=true, is_checking_x_agree=true, is_always_hp=false, init_solver="IPOPT", solver="IPOPT")

    Ff = [bop.F(x); bop.f(x)]
    if verbosity > 0
        if is_sol_valid
            print("success ")
        else
            print("FAIL ")
        end
        print("($status), $iter_count iters,\t", "x = $(round.(x, sigdigits=5)) -> Ff = $(round.(Ff, sigdigits=5)) (Ff* = $(round.(p.Ff_optimal, sigdigits=5)))\n\n")
    end

    (; is_sol_valid, x, λ, iter_count, status, Ff)
end

# LP-LP
include("./lp-lp/mb_2007_01.jl")
"""
all should solve
"""
function test_mb_2007_01(; tol=1e-7, n=10, rng=MersenneTwister(123))
    p = mb_2007_01()
    @testset "mb_2007_01" begin
        for _ = 1:n
            x_init = randn(rng, p.n1 + p.n2)
            is_sol_valid, x = call_solve_bop(p, x_init; tol)
            @test is_sol_valid && isapprox(x, p.xy_optimal; atol=tol)
        end
    end
end

include("./lp-lp/mb_2007_02.jl")
"""
no solution, none should solve
"""
function test_mb_2007_02(; tol=1e-7, n=10, rng=MersenneTwister(123))
    p = mb_2007_02()
    @testset "mb_2007_02" begin
        for _ = 1:n
            x_init = randn(rng, p.n1 + p.n2)
            is_sol_valid, _ = call_solve_bop(p, x_init; tol)
            @test !is_sol_valid
        end
    end
end

include("./lp-lp/bf_1982_02.jl")
"""
doesn't work with smaller tolerances
"""
function test_bf_1982_02(; tol=1e-5, n=10, rng=MersenneTwister(123))
    p = bf_1982_02()
    @testset "bf_1982_02" begin
        for _ = 1:n
            x_init = randn(rng, p.n1 + p.n2)
            is_sol_valid, x = call_solve_bop(p, x_init; tol)
            @test is_sol_valid && isapprox(x, p.xy_optimal; atol=tol)
        end
    end
end

include("./lp-lp/ct_1982_01.jl")
"""
broken
"""
function test_ct_1982_01(; tol=1e-7, n=3, rng=MersenneTwister(123))
    p = ct_1982_01()
    @testset "ct_1982_01" begin
        for _ = 1:n
            x_init = randn(rng, p.n1 + p.n2)
            is_sol_valid, x = call_solve_bop(p, x_init; tol)
            @test is_sol_valid && isapprox(x, p.xy_optimal; atol=tol) broken = true
        end
    end
end

# LP-QP
include("./lp-qp/mb_2006_01.jl")
"""
expected to converge to -1 (correct bilevel solution) and +1 (bilevel feasible)
"""
function test_mb_2006_01(; tol=1e-7, n=10, rng=MersenneTwister(123))
    p = mb_2006_01()
    @testset "mb_2006_01" begin
        x_init = [0.0]
        is_sol_valid, x = call_solve_bop(p, x_init; tol)
        @test is_sol_valid && (isapprox(x, p.xy_optimal; atol=tol) || isapprox(x, [1.0]; atol=tol))

        for _ = 1:n-1
            x_init = randn(rng, p.n1 + p.n2)
            is_sol_valid, x = call_solve_bop(p, x_init; tol)
            @test is_sol_valid && (isapprox(x, p.xy_optimal; atol=tol) || isapprox(x, [1.0]; atol=tol))
        end
    end
end

include("./lp-qp/mb_2007_03.jl")
"""
expected to converge to -1 (correct bilevel solution) and +1 (bilevel feasible)
"""
function test_mb_2007_03(; tol=1e-7, n=10, rng=MersenneTwister(123))
    p = mb_2007_03()
    @testset "mb_2007_03" begin
        x_init = [0.0]
        is_sol_valid, x = call_solve_bop(p, x_init; tol)
        @test is_sol_valid && (isapprox(x, p.xy_optimal; atol=tol) || isapprox(x, [1.0]; atol=tol))

        for _ = 1:n-1
            x_init = randn(rng, p.n1 + p.n2)
            is_sol_valid, x = call_solve_bop(p, x_init; tol)
            @test is_sol_valid && (isapprox(x, p.xy_optimal; atol=tol) || isapprox(x, [1.0]; atol=tol))
        end
    end
end


include("./lp-qp/mb_2007_04.jl")
"""
expected to converge to -.5 (bilevel infeasible, but local solution) and +1 (bilevel solution)
"""
function test_mb_2007_04(; tol=1e-7, n=10, rng=MersenneTwister(123))
    p = mb_2007_04()
    @testset "mb_2007_04" begin
        x_init = [0.0]
        is_sol_valid, x = call_solve_bop(p, x_init; tol)
        @test is_sol_valid && (isapprox(x, p.xy_optimal; atol=tol) || isapprox(x, [-0.5]; atol=tol))

        for _ = 1:n-1
            x_init = randn(rng, p.n1 + p.n2)
            is_sol_valid, x = call_solve_bop(p, x_init; tol)
            @test is_sol_valid && (isapprox(x, p.xy_optimal; atol=tol) || isapprox(x, [-0.5]; atol=tol))
        end
    end
end

include("./lp-qp/b_1991_02.jl")
"""
"""
function test_b_1991_02(; tol=1e-7, n=10, rng=MersenneTwister(123))
    p = b_1991_02()
    @testset "b_1991_02" begin
        for _ = 1:n
            x_init = [2.0; 0; 0] .+ rand(rng, p.n1 + p.n2) .* [2.0; 10.0; 10.0]
            is_sol_valid, x = call_solve_bop(p, x_init; tol)
            @test is_sol_valid && isapprox(x, p.xy_optimal; atol=tol)
        end
    end
end


include("./lp-qp/as_1984_01.jl")
"""
expected to converge to x=[0; 0; -10; -10] or [0; 30; -10; 10] or [20; 0; 0; -10] (local solutions) and x=[25; 30; 5; 10] (global solution
"""
function test_as_1984_01(; tol=1e-7, n=50, rng=MersenneTwister(123))
    p = as_1984_01()
    @testset "as_1984_01" begin
        for _ = 1:n
            x_init = [0.0; 0; -10.0; -10.0] .+ rand(rng, p.n1 + p.n2) .* [50.0; 50.0; 30.0; 30.0]
            is_sol_valid, x = call_solve_bop(p, x_init; tol)
            @test is_sol_valid && (isapprox(x, p.xy_optimal; atol=tol) || isapprox(x, [0.0; 0; -10; -10]; atol=tol) || isapprox(x, [0.0; 30; -10; 10]; atol=tol) || isapprox(x, [20.0; 0; 0; -10]; atol=tol))
        end
    end
end

# LP-NLP
include("./lp-nlp/mb_2007_05.jl")
"""
expected to converge to the wrong local min (-0.5) depending on the initialization
"""
function test_mb_2007_05(; tol=1e-7, n=10, rng=MersenneTwister(123))
    p = mb_2007_05()
    @testset "mb_2007_05" begin
        for _ = 1:n
            x_init = randn(rng, p.n1 + p.n2)
            is_sol_valid, x = call_solve_bop(p, x_init; tol)
            @test is_sol_valid && (isapprox(x, p.xy_optimal; atol=tol) || isapprox(x, [-0.5]; atol=tol))
        end
    end
end

include("./lp-nlp/mb_2007_06.jl")
"""
"""
function test_mb_2007_06(; tol=1e-7, n=10, rng=MersenneTwister(123))
    p = mb_2007_06()
    @testset "mb_2007_06" begin
        x_init = [-0.1]
        is_sol_valid, x = call_solve_bop(p, x_init; tol)
        @test is_sol_valid && isapprox(x, p.xy_optimal; atol=tol)

        # needs randomization
        x_init = [0.0]
        is_sol_valid, x = call_solve_bop(p, x_init; tol)
        @test is_sol_valid && isapprox(x, p.xy_optimal; atol=tol)

        x_init = [0.1]
        is_sol_valid, x = call_solve_bop(p, x_init; tol)
        @test is_sol_valid && isapprox(x, p.xy_optimal; atol=tol)

        for _ = 1:n-3
            x_init = randn(rng, p.n1 + p.n2)
            is_sol_valid, x = call_solve_bop(p, x_init; tol)
            @test is_sol_valid && isapprox(x, p.xy_optimal; atol=tol)
        end
    end
end


include("./lp-nlp/mb_2007_13.jl")
"""
the bilevel feasible set is large, but x=[-1;-1], [0;0] (large flat area so looser tol here), [1;1] are local solutions, while x=[0;1] is the bilevel solution
"""
function test_mb_2007_13(; tol=1e-7, n=50, rng=MersenneTwister(123))
    p = mb_2007_13()
    @testset "mb_2007_13" begin
        for _ = 1:n
            x_init = randn(rng, p.n1 + p.n2)
            is_sol_valid, x = call_solve_bop(p, x_init; tol)
            @test is_sol_valid && (isapprox(x, p.xy_optimal; atol=tol) || isapprox(x, [-1.0; -1]; atol=tol) || isapprox(x, [1.0; 1]; atol=tol) || isapprox(x, [0.0; 0]; atol=1e2 * tol))
        end
    end
end

include("./lp-nlp/mb_2007_13v.jl")
"""
different F from mb_2007_13, so [-1;-1] is the bilevel solution, and we somehow always converge to the optimal solution
"""
function test_mb_2007_13v(; tol=1e-7, n=50, rng=MersenneTwister(123))
    p = mb_2007_13v()
    @testset "mb_2007_13v" begin
        for _ = 1:n
            x_init = randn(rng, p.n1 + p.n2)
            is_sol_valid, x = call_solve_bop(p, x_init; tol)
            @test is_sol_valid && isapprox(x, p.xy_optimal; atol=tol)
        end
    end
end

include("./lp-nlp/ka_2014_01.jl")
"""
different f from mb_2007_13, so x=[-1;-1], [0;0] (large flat area so looser tol here), [1;1] are local solutions, while x=[0;1] is still the bilevel solution. curiously it also converges to [-1; .33], which is still a local solution
"""
function test_ka_2014_01(; tol=1e-7, n=50, rng=MersenneTwister(123))
    p = ka_2014_01()
    @testset "ka_2014_01" begin
        for _ = 1:n
            x_init = randn(rng, p.n1 + p.n2)
            is_sol_valid, x = call_solve_bop(p, x_init; tol)
            @test is_sol_valid && (isapprox(x, p.xy_optimal; atol=tol) || isapprox(x, [-1.0; -1]; atol=tol) || isapprox(x, [-1.0; 1 / 3]; atol=tol))
        end
    end
end

include("./qp-qp/as_1981_01.jl")
"""
broken 

sometimes converges to local sol [7.6271, 3.9406, 11.373, 17.059, 1.5677, 10.0, 28.432, 0], and others. can reach max random restarts
"""
function test_as_1981_01(; tol=1e-7, n=3, rng=MersenneTwister(123))
    p = as_1981_01()
    @testset "as_1981_01" begin
        for _ = 1:n
            x_init = rand(rng, p.n1 + p.n2) .* [10.0; 5; 15; 20; 20; 20; 40; 40]
            is_sol_valid, x = call_solve_bop(p, x_init; tol)
            #print("$x_init->$x\n")
            @test is_sol_valid && (isapprox(x, p.xy_optimal; atol=tol) || isapprox(x, [7.6270726902828185, 3.9406090120498836, 11.372927321048655, 17.05939098652773, 1.567681684379581, 10.00000003777213, 28.432318327395397, 0]; atol=1e2 * tol))
        end
    end
end


include("./qp-qp/d_1978_01.jl")
"""
todo
"""

include("./qp-qp/y_1996_02.jl")
"""
todo
"""

include("./nlp-nlp/ka_2014_02.jl")
"""
"""
function test_ka_2014_02(; tol=1e-7, n=3, rng=MersenneTwister(123))
    p = ka_2014_02()
    @testset "ka_2014_02" begin
        for _ = 1:n
            x_init = randn(rng, p.n1 + p.n2)
            is_sol_valid, x = call_solve_bop(p, x_init; tol)
            print("$x_init->$x\n")
            @test is_sol_valid && isapprox(x, p.xy_optimal; atol=tol) broken = true
        end
    end
end

function run_all_tests()
    # LP-LP:
    print("LP-LP:\n")
    test_mb_2007_01()
    test_mb_2007_02()
    test_bf_1982_02()
    test_ct_1982_01()
    # LP-QP:
    print("LP-QP:\n")
    test_mb_2006_01()
    test_mb_2007_03()
    test_mb_2007_04()
    test_b_1991_02()
    test_as_1984_01()
    # LP-NLP:
    print("LP-NLP:\n")
    test_mb_2007_05()
    test_mb_2007_06()
    test_mb_2007_13()
    test_mb_2007_13v()
    test_ka_2014_01()
    # QP-QP:
    print("QP-QP:\n")
    test_as_1981_01()
    # NLP-NLP:
    print("NLP-NLP:\n")
    test_ka_2014_02()
end

run_all_tests()