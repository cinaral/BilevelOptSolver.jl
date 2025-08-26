### deprecated BASBLib tests
using BilevelOptSolver
using Test
using Random
# 2025-08-06 TODO: add randomization, add forced randomization, fix x agreement logic

function call_solve_bop(p, x_init; verbosity=0, tol=1e-7)
    if verbosity > 0
        print("$(p.name):\n")
    end

    bop, syms = construct_bop(p.n1, p.n2, p.F, p.G, p.f, p.g; verbosity=0, np=0)
    is_sol_valid, x, λ, iter_count, status = solve_bop(bop; max_iter=50, x_init, verbosity, tol, conv_dv_len=3, do_check_x_agreem=true, do_force_hp_init=false, init_solver="PATH", solver="PATH")

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
converges to local sol [0.5, 0.5, 0, 0, 0, 0.48687056090227804 ,0 ,0]
"""
function test_ct_1982_01(; tol=1e-7, n=3, rng=MersenneTwister(123))
    p = ct_1982_01()
    @testset "ct_1982_01" begin
        for _ = 1:n
            x_init = randn(rng, p.n1 + p.n2)
            is_sol_valid, x = call_solve_bop(p, x_init; tol)
            @test is_sol_valid && isapprox(x, p.xy_optimal; atol=tol) broken=true
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
        @test is_sol_valid && isapprox(x, p.xy_optimal; atol=tol)

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
function test_as_1984_01(; tol=1e-6, n=50, rng=MersenneTwister(123))
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
        print("$is_sol_valid $x_init->$x\n")
        @test is_sol_valid && isapprox(x, p.xy_optimal; atol=tol)

        # needs randomization
        x_init = [0.0]
        is_sol_valid, x = call_solve_bop(p, x_init; tol)
        print("$is_sol_valid $x_init->$x\n")
        @test is_sol_valid && isapprox(x, p.xy_optimal; atol=tol)

        x_init = [0.1]
        is_sol_valid, x = call_solve_bop(p, x_init; tol)
        print("$is_sol_valid $x_init->$x\n")
        @test is_sol_valid && isapprox(x, p.xy_optimal; atol=tol)

        for _ = 1:n-3
            x_init = randn(rng, p.n1 + p.n2)
            is_sol_valid, x = call_solve_bop(p, x_init; tol)
            print("$is_sol_valid $x_init->$x\n")
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
    #test_mb_2007_01()
    #test_mb_2007_02()
    #test_bf_1982_02()
    #test_ct_1982_01()
    # LP-QP:
    print("LP-QP:\n")
    #test_mb_2006_01()
    #test_mb_2007_03()
    #test_mb_2007_04()
    #test_b_1991_02()
    #test_as_1984_01()
    # LP-NLP:
    print("LP-NLP:\n")
    #test_mb_2007_05()
    #test_mb_2007_06()
    #test_mb_2007_13()
    #test_mb_2007_13v()
    #test_ka_2014_01()
    # QP-QP:
    print("QP-QP:\n")
    #test_as_1981_01()
    # NLP-NLP:
    print("NLP-NLP:\n")
    #test_ka_2014_02()
end

run_all_tests()


###### 2025-08-20 deprecated reference
struct BilevelOptProb_2
    F::Function # F(x)
    G::Function # G(x)
    f::Function # f(x)
    g::Function # g(x)
    n1::Int # length(x₁)
    n2::Int # length(x₂)
    m1::Int # length(G(x))
    m2::Int # length(g(x))
    nx::Int # length(x) = n₁ + n₂
    nz::Int # length(z) = n₂ + m₂, z := [x₂; λ] s.t. h ⟂ zu ≥ z ≥ zl 
    n::Int  # length(v) = n₁ + nz, v := [x₁; z]
    m::Int  # length(Γ) = m₁ + 4*nz, Γ := [G; h; -h; z; -z] ≥ Γlᵢ (hlᵢ, huᵢ, zlᵢ, zuᵢ are computed later)
    nθ::Int # length(θ) = n + m, θ := [v; Λ] s.t. Φ ⟂ θu ≥ θ ≥ θl  
    np::Int # (optional) # of appended non-decision variables (parameters)
    # follower NLP: min f(x) s.t. g(x) ≥ 0
    g!::Function # g!(out, x) 
    ∇ₓ₂f!::Function # ∇ₓ₂f!(out, x)
    # ∇ₓ₂g
    ∇ₓ₂g_size::Tuple{Int,Int}
    ∇ₓ₂g_rows::Vector{Int}
    ∇ₓ₂g_cols::Vector{Int}
    ∇ₓ₂g_vals!::Function # ∇ₓ₂g_vals!(out, x)
    # ∇²ₓ₂L2, L2 = of * f(x) - λ' g(x)
    ∇²ₓ₂L2_size::Tuple{Int,Int}
    ∇²ₓ₂L2_rows::Vector{Int}
    ∇²ₓ₂L2_cols::Vector{Int}
    ∇²ₓ₂L2_vals!::Function # ∇²ₓ₂L_vals!(out, v, λ, of)
    # follower MCP: z := [x₂; λ] s.t. h ⟂ zu ≥ z ≥ zl
    zl₀::Vector{Float64} # default z bounds 
    zu₀::Vector{Float64}
    h!::Function # h!(out, v)
    # ∇_z h
    ∇h_size::Tuple{Int,Int}
    ∇h_rows::Vector{Int}
    ∇h_cols::Vector{Int}
    ∇h_vals!::Function # ∇h_vals!(out, v)
    # SBOPi NLP: min F(v) s.t. Γ := [G; h; -h; z; -z] ≥ Γlᵢ 
    Fv::Function # F(v)
    Γ!::Function # Γ!(out, v)
    ∇ᵥF!::Function # Γ!(out, v) 
    # ∇ᵥΓ
    ∇ᵥΓ_size::Tuple{Int,Int}
    ∇ᵥΓ_rows::Vector{Int}
    ∇ᵥΓ_cols::Vector{Int}
    ∇ᵥΓ_vals!::Function # ∇ᵥΓ_vals!(out, v)
    # ∇²ᵥL1, L1 = of * F(v) - Λ' Γ(v)
    ∇²ᵥL1_size::Tuple{Int,Int}
    ∇²ᵥL1_rows::Vector{Int}
    ∇²ᵥL1_cols::Vector{Int}
    ∇²ᵥL1_vals!::Function # ∇²ᵥL1_vals!(out, v, Λ, of)
    # SBOPi MCP: θ := [v; Λ; r] s.t. Φ ⟂ θu ≥ θ ≥ θl
    θl₀::Vector{Float64} # default θ bounds
    θu₀::Vector{Float64}
    Φ!::Function # Φ!(out, v, Λ)
    # ∇_θ Φ
    ∇Φ_size::Tuple{Int,Int}
    ∇Φ_rows::Vector{Int}
    ∇Φ_cols::Vector{Int}
    ∇Φ_vals!::Function # ∇Φ_vals!(out, θ)
    # high-point NLP: min F(x) s.t. [G(x); g(x)] ≥ 0
    Gg!::Function # Gg!(out, x) 
    ∇ₓF!::Function # ∇ₓF!(out, x)
    # ∇ₓGg
    ∇ₓGg_size::Tuple{Int,Int}
    ∇ₓGg_rows::Vector{Int}
    ∇ₓGg_cols::Vector{Int}
    ∇ₓGg_vals!::Function # ∇ₓGg_vals!(out, x)
    # ∇²ₓL3, L3 = of * F(x) - Λ₃' Gg(x)
    ∇²ₓL3_size::Tuple{Int,Int}
    ∇²ₓL3_rows::Vector{Int}
    ∇²ₓL3_cols::Vector{Int}
    ∇²ₓL3_vals!::Function # ∇²ₓ₃L_vals!(out, x, Λ₃, of)
    # for convenience
    inds
end

function construct_bop_2(n1, n2, F, G, f, g; np=0, verbosity=0)
    nx = n1 + n2 # length(x)
    xp_dummy = zeros(nx + np) # optionally for p>0, xp would be appended by np parameters
    m1 = length(G(xp_dummy))
    m2 = length(g(xp_dummy))
    nz = n2 + m2
    n = n1 + nz
    m = m1 + 4 * nz
    nθ = n + 2 * m

    x1_sym = Symbolics.@variables(x[1:n1])[1] |> Symbolics.scalarize
    x2_sym = Symbolics.@variables(y[1:n2])[1] |> Symbolics.scalarize
    λ_sym = Symbolics.@variables(λ[1:m2])[1] |> Symbolics.scalarize
    p_sym = Symbolics.@variables(p[1:np])[1] |> Symbolics.scalarize
    x_sym = [x1_sym; x2_sym] # x := [x₁; x₂]
    z_sym = [x2_sym; λ_sym] # z := [x₂; λ]
    v_sym = [x1_sym; z_sym] # z := [x₁; z]
    @assert(nx == length(x_sym))
    @assert(nz == length(z_sym))
    @assert(n == length(v_sym))

    xp_sym = Num[x_sym; p_sym]
    vp_sym = Num[v_sym; p_sym]
    F_sym = F(xp_sym)
    G_sym = G(xp_sym)
    f_sym = f(xp_sym)
    g_sym = g(xp_sym)

    # we define these index dictionaries for convenience later
    inds = define_index_dicts(n1, n2, m1, m2, nx, np, nz, n, m)

    f = Symbolics.build_function(f_sym, xp_sym; expression=Val{false})
    g! = Symbolics.build_function(g_sym, xp_sym; expression=Val{false})[2]
    ∇ₓ₂f_sym = Symbolics.gradient(f_sym, x2_sym)
    ∇ₓ₂f! = Symbolics.build_function(∇ₓ₂f_sym, xp_sym; expression=Val{false})[2]

    # ∇ₓ₂g
    ∇ₓ₂g_sym = Symbolics.sparsejacobian(g_sym, x2_sym)
    ∇ₓ₂g_size = size(∇ₓ₂g_sym)
    (∇ₓ₂g_rows, ∇ₓ₂g_cols, ∇ₓ₂g_vals_sym) = SparseArrays.findnz(∇ₓ₂g_sym)
    ∇ₓ₂g_vals! = Symbolics.build_function(∇ₓ₂g_vals_sym, xp_sym; expression=Val{false})[2]

    # ∇²ₓ₂L2, L2 =  of f(x) - λ' g(x)
    of2_sym = Symbolics.@variables(λ₀)[1] # objective factor
    if isempty(λ_sym)
        L2_sym = of2_sym * f_sym
    else
        L2_sym = of2_sym * f_sym - g_sym' * λ_sym
    end
    ∇ₓ₂L2_sym = Symbolics.gradient(L2_sym, x2_sym)
    ∇²ₓ₂L2_sym = Symbolics.sparsejacobian(∇ₓ₂L2_sym, x2_sym)
    ∇²ₓ₂L2_size = size(∇²ₓ₂L2_sym)
    (∇²ₓ₂L2_rows, ∇²ₓ₂L2_cols, ∇²ₓ₂L2_vals_sym) = SparseArrays.findnz(∇²ₓ₂L2_sym)
    ∇²ₓ₂L2_vals! = Symbolics.build_function(∇²ₓ₂L2_vals_sym, xp_sym, λ_sym, of2_sym; expression=Val{false})[2]

    # follower MCP: z := [x₂; λ] s.t. h ⟂ zu ≥ z ≥ zl
    # by default x₂ is free and λ ≥ 0, but these will be overwritten later
    ∇ₓ₂L2_sym = substitute(∇ₓ₂L2_sym, Dict([of2_sym => 1.0]))
    zl₀ = [fill(-Inf, n2); zeros(m2)] # z lb does not change
    zu₀ = fill(Inf, nz) # default z ub
    if nz > 0
        h_sym = [∇ₓ₂L2_sym; g_sym]
    else
        h_sym = []
    end
    @assert(nz == length(zl₀))
    @assert(nz == length(zu₀))
    @assert(nz == length(h_sym))
    h! = Symbolics.build_function(h_sym, vp_sym; expression=Val(false))[2]

    # ∇_z h
    ∇h_sym = Symbolics.sparsejacobian(h_sym, z_sym) # Jacobian of F for the PATH solver
    ∇h_size = size(∇h_sym)
    (∇h_rows, ∇h_cols, ∇h_vals_sym) = SparseArrays.findnz(∇h_sym)
    ∇h_vals! = Symbolics.build_function(∇h_vals_sym, vp_sym; expression=Val{false})[2]

    # SBOPi NLP: min F(v) s.t. Γ := [G; h; z; -h; -z] ≥ Γlᵢ 
    # we append [-h; -z] to keep the corresponding Λ ≥ 0 without specifying upper bounds
    Γ_sym = [G_sym; h_sym; z_sym; -h_sym; -z_sym]
    @assert(m == length(Γ_sym))
    Fv = Symbolics.build_function(F_sym, vp_sym; expression=Val{false}) # for convenience: F that takes v as argument 
    Γ! = Symbolics.build_function(Γ_sym, vp_sym; expression=Val(false))[2]
    ∇ᵥF_sym = Symbolics.gradient(F_sym, vp_sym)
    ∇ᵥF! = Symbolics.build_function(∇ᵥF_sym, vp_sym; expression=Val{false})[2]

    # ∇ᵥΓ
    ∇ᵥΓ_sym = Symbolics.sparsejacobian(Γ_sym, v_sym)
    ∇ᵥΓ_size = size(∇ᵥΓ_sym)
    (∇ᵥΓ_rows, ∇ᵥΓ_cols, ∇ᵥΓ_vals_sym) = SparseArrays.findnz(∇ᵥΓ_sym)
    ∇ᵥΓ_vals! = Symbolics.build_function(∇ᵥΓ_vals_sym, vp_sym; expression=Val{false})[2]

    # ∇²ᵥL1, L1 = F(v) - Λ' Γ(v)
    Λ_sym = Symbolics.@variables(Λ[1:m])[1] |> Symbolics.scalarize
    of1_sym = Symbolics.@variables(λ₀)[1] # objective factor
    if isempty(λ_sym)
        L1_sym = of1_sym * F_sym
    else
        L1_sym = of1_sym * F_sym - Γ_sym' * Λ_sym
    end
    ∇ᵥL1_sym = Symbolics.gradient(L1_sym, v_sym)
    ∇²ᵥL1_sym = Symbolics.sparsejacobian(∇ᵥL1_sym, v_sym)
    ∇²ᵥL1_size = size(∇²ᵥL1_sym)
    (∇²ᵥL1_rows, ∇²ᵥL1_cols, ∇²ᵥL1_vals_sym) = SparseArrays.findnz(∇²ᵥL1_sym)
    ∇²ᵥL1_vals! = Symbolics.build_function(∇²ᵥL1_vals_sym, vp_sym, Λ_sym, of1_sym; expression=Val{false})[2]

    # SBOPi MCP: θ := [v; Λ; r] s.t. Φ ⟂ θu ≥ θ ≥ θl
    # in order to conveniently specify the bounds the Γ := [G; h; z; -h; -z], we introduce its slack, r ⟂ Λ
    # this way v and Λ is free, and by default r ≥ 0, but r bounds will be overwritten later
    # TODO 2025-07-14: this is quite a few more problem variables than needed, but we wish to keep Λ identical to the NLP interpretation for the purposes of checking optimality conditions
    r_sym = Symbolics.@variables(r[1:m])[1] |> Symbolics.scalarize
    θ_sym = [v_sym; Λ_sym; r_sym]
    @assert(nθ == length(θ_sym))
    θp_sym = Num[v_sym; Λ_sym; r_sym; p_sym]
    θl₀ = [fill(-Inf, n + m); zeros(m)] # default θ lb
    θu₀ = fill(Inf, nθ) # default θ ub
    ∇ᵥL1_sym = substitute(∇ᵥL1_sym, Dict([of1_sym => 1.0]))
    Γ_eq_sym = Γ_sym .- r_sym
    Φ_sym = [∇ᵥL1_sym; Γ_eq_sym; Λ_sym]
    Φ! = Symbolics.build_function(Φ_sym, θp_sym; expression=Val(false))[2]
    ∇Φ_sym = Symbolics.sparsejacobian(Φ_sym, θ_sym)
    ∇Φ_size = size(∇Φ_sym)
    (∇Φ_rows, ∇Φ_cols, ∇Φ_vals_sym) = SparseArrays.findnz(∇Φ_sym)
    ∇Φ_vals! = Symbolics.build_function(∇Φ_vals_sym, θp_sym; expression=Val{false})[2]

    # high-point NLP: min F(x) s.t. [G(x); g(x)] ≥ 0
    ∇ₓF_sym = Symbolics.gradient(F_sym, x_sym)
    ∇ₓF! = Symbolics.build_function(∇ₓF_sym, xp_sym; expression=Val{false})[2]
    Gg_sym = [G_sym; g_sym]
    Gg! = Symbolics.build_function(Gg_sym, xp_sym; expression=Val{false})[2]

    # ∇ₓGg
    ∇ₓGg_sym = Symbolics.sparsejacobian(Gg_sym, x_sym)
    ∇ₓGg_size = size(∇ₓGg_sym)
    (∇ₓGg_rows, ∇ₓGg_cols, ∇ₓGg_vals) = SparseArrays.findnz(∇ₓGg_sym)
    ∇ₓGg_vals! = Symbolics.build_function(∇ₓGg_vals, xp_sym; expression=Val{false})[2]

    # ∇²ₓL3, L3 = F(x) - Λ₃' Gg(x)
    Λ₃_sym = Symbolics.@variables(Λ[1:m1+m2])[1] |> Symbolics.scalarize
    of3_sym = Symbolics.@variables(λ₀)[1] # objective factor
    if isempty(Λ₃_sym)
        L3_sym = of3_sym * F_sym
    else
        L3_sym = of3_sym * F_sym - Gg_sym' * Λ₃_sym
    end
    ∇ₓL3_sym = Symbolics.gradient(L3_sym, x_sym)
    ∇²ₓL3_sym = Symbolics.sparsejacobian(∇ₓL3_sym, x_sym)
    ∇²ₓL3_size = size(∇²ₓL3_sym)
    (∇²ₓL3_rows, ∇²ₓL3_cols, ∇²ₓL3_vals_sym) = SparseArrays.findnz(∇²ₓL3_sym)
    ∇²ₓL3_vals! = Symbolics.build_function(∇²ₓL3_vals_sym, xp_sym, Λ₃_sym, of3_sym; expression=Val{false})[2]

    bop = BilevelOptProb_2(
        F,
        G,
        f,
        g,
        n1,
        n2,
        m1,
        m2,
        nx,
        nz,
        n,
        m,
        nθ,
        np,
        # follower NLP: min f(x) s.t. g(x) ≥ 0
        g!,
        ∇ₓ₂f!,
        # ∇ₓ₂g
        ∇ₓ₂g_size,
        ∇ₓ₂g_rows,
        ∇ₓ₂g_cols,
        ∇ₓ₂g_vals!,
        # ∇²ₓ₂L2, L2 = f(x) - λ' g(x)
        ∇²ₓ₂L2_size,
        ∇²ₓ₂L2_rows,
        ∇²ₓ₂L2_cols,
        ∇²ₓ₂L2_vals!,
        # follower MCP: z := [x₂; λ] s.t. h ⟂ zu ≥ z ≥ zl
        zl₀,
        zu₀,
        h!,
        # ∇_z h
        ∇h_size,
        ∇h_rows,
        ∇h_cols,
        ∇h_vals!,
        # SBOPi NLP: min F(v) s.t. Γ := [G; h; -h; z; -z] ≥ Γlᵢ 
        Fv,
        Γ!,
        ∇ᵥF!,
        # ∇ᵥΓ
        ∇ᵥΓ_size,
        ∇ᵥΓ_rows,
        ∇ᵥΓ_cols,
        ∇ᵥΓ_vals!,
        # ∇²ᵥL1, L1 = F(v) - Λ' Γ(v)
        ∇²ᵥL1_size,
        ∇²ᵥL1_rows,
        ∇²ᵥL1_cols,
        ∇²ᵥL1_vals!,
        # SBOPi MCP: θ := [v; Λ] s.t. Φ ⟂ θu ≥ θ ≥ θl
        θl₀,
        θu₀,
        Φ!,
        # ∇_θ Φ
        ∇Φ_size,
        ∇Φ_rows,
        ∇Φ_cols,
        ∇Φ_vals!,
        # high-point NLP: min F(x) s.t. [G(x); g(x)] ≥ 0
        Gg!,
        ∇ₓF!,
        # ∇ₓGg
        ∇ₓGg_size,
        ∇ₓGg_rows,
        ∇ₓGg_cols,
        ∇ₓGg_vals!,
        # ∇²ₓL3, L3 = F(x) - Λ₃' Gg(x)
        ∇²ₓL3_size,
        ∇²ₓL3_rows,
        ∇²ₓL3_cols,
        ∇²ₓL3_vals!,
        # for convenience
        inds
    )

    syms = (; x=xp_sym, F=F_sym, G=G_sym, f=f_sym, g=g_sym, ∇ₓ₂f=∇ₓ₂f_sym, ∇ₓ₂g=∇ₓ₂g_sym, λ=λ_sym, L2=L2_sym, ∇ₓ₂L2=∇ₓ₂L2_sym, ∇²ₓ₂L2=∇²ₓ₂L2_sym, v=vp_sym, z=z_sym, h=h_sym, ∇h=∇h_sym, Γ=Γ_sym, ∇ᵥF=∇ᵥF_sym, ∇ᵥΓ=∇ᵥΓ_sym, L1=L1_sym, ∇ᵥL1=∇ᵥL1_sym, ∇²ᵥL1=∇²ᵥL1_sym, θ=θp_sym, Φ=Φ_sym, ∇Φ=∇Φ_sym, ∇ₓF=∇ₓF_sym, Gg=Gg_sym, ∇ₓGg=∇ₓGg_sym, L3=L3_sym, ∇ₓL3=∇ₓL3_sym, ∇²ₓL3=∇²ₓL3_sym)

    (; bop, syms)
end


# 2025-08-22

"""
[IPOPT](https://coin-or.github.io/Ipopt/) solves:
```
min     f(x)
 x
s.t.    x_l ≤ x ≤ x_u
        g_l ≤ g(x) ≤ g_u
```
"""
function setup_nlp_solve_IPOPT(n, m, xl, xu, gl, gu, f, g!, ∇ₓf!, ∇ₓg_rows, ∇ₓg_cols, ∇ₓg_vals!, ∇ₓₓL_rows, ∇ₓₓL_cols, ∇ₓₓL_vals!)
    #@assert(n == length(xl), "wrong x_l dimensions")
    #@assert(n == length(xu), "wrong x_u dimensions")
    #@assert(m == length(gl), "wrong g_l dimensions")
    #@assert(m == length(gu), "wrong g_u dimensions")

    n_nz_∇ₓg = length(∇ₓg_rows)
    @assert(n_nz_∇ₓg == length(∇ₓg_cols), "∇ₓg_rows and ∇ₓg_cols dimensions do not match")

    n_nz_∇ₓₓL = length(∇ₓₓL_rows)
    @assert(n_nz_∇ₓₓL == length(∇ₓₓL_cols), "∇ₓₓL_rows and ∇ₓₓL_cols dimensions do not match")

    function eval_f(x::Vector{Float64})
        f(x)
    end

    function eval_g(x::Vector{Float64}, g::Vector{Float64})
        g!(g, x)
    end

    function eval_∇ₓf(x::Vector{Float64}, ∇ₓf::Vector{Float64})
        ∇ₓf!(∇ₓf, x)
    end

    function eval_∇ₓg(x::Vector{Float64}, rows::Vector{Int32}, cols::Vector{Int32}, vals::Union{Nothing,Vector{Float64}})
        if vals === nothing
            rows .= ∇ₓg_rows
            cols .= ∇ₓg_cols
        else
            ∇ₓg_vals!(vals, x)
        end
    end

    function eval_∇ₓₓL(
        x::Vector{Float64},
        rows::Vector{Int32},
        cols::Vector{Int32},
        obj_factor::Float64,
        λ::Vector{Float64},
        vals::Union{Nothing,Vector{Float64}},
    )
        if vals === nothing
            rows .= ∇ₓₓL_rows
            cols .= ∇ₓₓL_cols
        else
            ∇ₓₓL_vals!(vals, x, -λ, obj_factor) # convention change: IPOPT expects L = obj_factor * f(x) + g(x)' λ
        end
    end

    # we return this function
    function solve_nlp(; gl=gl, gu=gu, x_init=zeros(n), λ_init=zeros(m), tol=1e-6, max_iter=1000, print_level=0, is_using_HSL=false)
        ipopt_prob = Ipopt.CreateIpoptProblem(
            n,
            xl,
            xu,
            m,
            gl,
            gu,
            n_nz_∇ₓg,
            n_nz_∇ₓₓL,
            eval_f,
            eval_g,
            eval_∇ₓf,
            eval_∇ₓg,
            eval_∇ₓₓL
        )

        Ipopt.AddIpoptNumOption(ipopt_prob, "tol", tol)
        Ipopt.AddIpoptIntOption(ipopt_prob, "max_iter", max_iter)
        Ipopt.AddIpoptIntOption(ipopt_prob, "print_level", print_level)
        if is_using_HSL && !haskey(ENV, "HSL_PATH")
            is_using_HSL = false
            if verbosity > 0
                print("HSL_PATH not found: Defaulting is_using_HSL = false. If you would like to use HSL, please obtain a license and download HSL_jll.jl (https://licences.stfc.ac.uk/products/Software/HSL/LibHSL), and set HSL_PATH environment variable to the extracted location.\n")
            end
        end    
        if is_using_HSL
            Ipopt.AddIpoptStrOption(ipopt_prob, "hsllib", HSL_jll.libhsl_path)
            Ipopt.AddIpoptStrOption(ipopt_prob, "linear_solver", "ma27")
        end

        ipopt_prob.x .= x_init
        ipopt_prob.mult_g .= -λ_init # convention change!!
        solvestat = Ipopt.IpoptSolve(ipopt_prob)

        x = ipopt_prob.x
        λ = -ipopt_prob.mult_g # convention change!!
        (; x, λ, solvestat, ipopt_prob)
    end
    solve_nlp
end

"""
[PATHSolver](https://coin-or.github.io/Ipopt/) solves:
```
x s.t. at least one of the following holds: 
    1.  F_i(x) = 0, x_l_i <= x_i <= x_u_i
    2.  F_i(x) > 0, x_l_i =  x_i  
    3.  F_i(x) < 0,          x_i = x_u_i

J is the Jacobian: ∇ₓF
```
"""
function setup_mcp_solve_PATH(n, xl, xu, F, J_rows, J_cols, J_vals!; verbosity=0)
    #@assert(n == length(x_l), "wrong x_l dimensions")
    #@assert(n == length(x_u), "wrong x_u dimensions")

    n_nz_J = length(J_rows)
    #@assert(n_nz_J == length(J_cols), "J_rows and J_cols dimensions do not match")

    J_shape = sparse(J_rows, J_cols, ones(Cdouble, n_nz_J), n, n)
    J_col = J_shape.colptr[1:end-1]
    J_len = diff(J_shape.colptr)
    J_row = J_shape.rowval

    # 2025-07-13 Why do I need this??
    #if isempty(J_row)
    #    J_row = 0
    #end

    function eval_F(n, x, vals)
        vals .= 0.0
        F(vals, x)
        Cint(0)
    end

    function eval_J(n, nnz, x, col, len, row, vals)
        vals .= 0.0
        J_vals!(vals, x)
        col .= J_col
        len .= J_len
        row .= J_row
        Cint(0)
    end

    if n > 300 && !haskey(ENV, "PATH_LICENSE_STRING")
        if verbosity > 0
            print("PATH_LICENSE_STRING not found and problem size is too large: Please obtain a license (https://pages.cs.wisc.edu/~ferris/path/julia/LICENSE), and set PATH_LICENSE_STRING environment variable.\n")
        end
    end

    if haskey(ENV, "PATH_LICENSE_STRING")
        PATHSolver.c_api_License_SetString(ENV["PATH_LICENSE_STRING"])
    end

    function solve_mcp(; xl=xl, xu=xu, x_init=zeros(n), tol=1e-6, max_iter=500, is_silent=true, verbosity=5)
        status, x, info = PATHSolver.solve_mcp(
            eval_F,
            eval_J,
            xl,
            xu,
            x_init,
            silent=is_silent,
            nnz=n_nz_J,
            jacobian_structure_constant=true,
            jacobian_data_contiguous=true,
            cumulative_iteration_limit=20_000,
            major_iteration_limit=max_iter,
            time_limit=5,
            convergence_tolerance=tol,
            lemke_rank_deficiency_iterations=30, # fixes silent crashes " ** SOLVER ERROR ** Lemke: invertible basis could not be computed."
            preprocess=1,
            presolve=1,
            output=verbosity,
            output_options=0,
            output_errors=0,
            output_warnings=0,
            output_final_summary=1
        )

        (; x, status, info)
    end
    solve_mcp
end

"""
[HiGHS](https://ergo-code.github.io/HiGHS/dev/) solves:
```
    min     cᵀx + d   
    x
    s.t.    x_l ≤ x ≤ x_u
            A_l ≤ A x ≤ A_u
```
Because we're only checking for LP feasiblity: c = d = 0.

Here's a HiGHS encoding example because the API is kinda obscure:
```
    Min     x₀ +  x₁ + 3
    x
    s.t.    0 ≤ x₀ ≤ 4
            1 ≤ x₁   
                      x₁ ≤ 7
            5 ≤ x₀ + 2x₁ ≤ 15
            6 ≤ 3x₀ + 2x₁
```
This would be encoded like this:
```
c = [1.0, 1.0] (col_cost)
offset = 3.;
x_l = [0.0, 1.0] (col_lower)
x_u = [4.0, Inf] (col_upper)
A_l = [-Inf, 5.0, 6.0] (row_lower)
A_u = [7.0, 15.0, Inf] (row_upper)
a_start = [0, 2] (column-wise, column start index, the first index is always zero)
a_index = [1, 2, 0, 1, 2]
a_value = [1.0, 3.0, 1.0, 2.0, 2.0]
```
"""
function setup_lp_feas_check_HiGHS(n; primal_feas_tol=1e-6, zero_tol=1e-3, output_flag=false, verbosity=0)
    c = zeros(n)
    offset = 0.0

    model = HiGHS.Highs_create()
    HiGHS.Highs_setDoubleOptionValue(model, "primal_feasibility_tolerance",
        primal_feas_tol)
    HiGHS.Highs_setBoolOptionValue(model, "output_flag", output_flag)
    #HiGHS.Highs_setBoolOptionValue(model, "log_to_console", true)

    function check_feas(x_l, x_u, A_l, A_u, A)
        # constraint matrix is column-wise:
        a_start = A.colptr[1:end-1] .- 1
        a_index = A.rowval .- 1
        a_value = A.nzval

        n_x = convert(Cint, n)
        n_A_l = convert(Cint, size(A_l, 1))
        n_nz = convert(Cint, size(a_index, 1))
        c = convert(Array{Cdouble}, c)
        x_l = convert(Array{Cdouble}, x_l)
        x_u = convert(Array{Cdouble}, x_u)
        offset = convert(Cdouble, offset)
        A_l = convert(Array{Cdouble}, A_l)
        A_u = convert(Array{Cdouble}, A_u)
        a_start = convert(Array{Cint}, a_start)
        a_index = convert(Array{Cint}, a_index)
        a_value = convert(Array{Cdouble}, a_value)

        status = HiGHS.Highs_passLp(
            model,
            n_x,
            n_A_l,
            n_nz,
            HiGHS.kHighsMatrixFormatColwise,
            HiGHS.kHighsObjSenseMinimize,
            offset,
            c,
            x_l,
            x_u,
            A_l,
            A_u,
            a_start,
            a_index,
            a_value
        )

        if status != HiGHS.kHighsStatusOk
            if verbosity > 1
                print("Failed passing model to HiGHS!\n")
            end
            return false
        end

        if HiGHS.Highs_run(model) != HiGHS.kHighsStatusOk
            if verbosity > 1
                print("Failed running the model!\n")
            end
            return false
        end

        primal_sol_status = Ref{Cint}(0)
        HiGHS.Highs_getIntInfoValue(model, "primal_solution_status", primal_sol_status)

        is_primal_feasible = primal_sol_status[] == HiGHS.kHighsSolutionStatusFeasible

        if is_primal_feasible
            return true
        else
            return false
        end
    end
    check_feas
end




###### 


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