using Symbolics
using SparseArrays
using Ipopt
using BenchmarkTools
using ProfileView

import Pkg
Pkg.develop(path = "./HSL_jll.jl.v2024.11.28")

import HSL_jll

# example from
# https://nlopt.readthedocs.io/en/latest/NLopt_Tutorial/#example-nonlinearly-constrained-problem
n::Int64 = 2

function f(x)
    sqrt(x[2])
end

function g(x)
    [
        x[2] - (2 * x[1] + 0)^3
        x[2] - (-1 * x[1] + 1)^3
    ]
end

function create_prob(n, f, g)
    x_dummy = zeros(n)
    m = length(g(x_dummy))

    x_sym = Symbolics.@variables(x[1:n])[1] |> Symbolics.scalarize
    f_sym = f(x_sym)
    g_sym = g(x_sym)

    λ_sym = Symbolics.@variables(λ[1:m])[1] |> Symbolics.scalarize

    eval_f = Symbolics.build_function(f_sym, x_sym; expression=Val{false})
    eval_g! = Symbolics.build_function(g_sym, x_sym; expression=Val{false})[2]
    ∇ₓf = Symbolics.gradient(f_sym, x)
    eval_∇ₓf! = Symbolics.build_function(∇ₓf, x_sym; expression=Val{false})[2]

    ∇ₓg = Symbolics.sparsejacobian(g_sym, x)
    (∇ₓg_rows, ∇ₓg_cols, ∇ₓg_vals) = SparseArrays.findnz(∇ₓg)
    eval_∇ₓg_vals! = Symbolics.build_function(∇ₓg_vals, x_sym; expression=Val{false})[2]
    eval_∇ₓg = (; shape=size(∇ₓg), rows=∇ₓg_rows, cols=∇ₓg_cols, vals=eval_∇ₓg_vals!)

    obj_factor = Symbolics.@variables(σf)[1]
    if isempty(λ_sym)
        L_follow = obj_factor * f_sym # WARN: IPOPT convention: ∇²ₓ₂f(x) + λᵀ ∇²ₓ₂ g(x)
    else
        L_follow = obj_factor * f_sym + g_sym' * λ_sym # WARN: IPOPT convention: ∇²ₓ₂f(x) + λᵀ ∇²ₓ₂ g(x)
    end

    ∇ₓL = Symbolics.gradient(L_follow, x)
    ∇²ₓL = Symbolics.sparsejacobian(∇ₓL, x)
    (∇²ₓL_rows, ∇²ₓL_cols, ∇²ₓL_vals) = SparseArrays.findnz(∇²ₓL)
    eval_∇²ₓL_vals! = Symbolics.build_function(∇²ₓL_vals, x_sym, obj_factor, λ_sym; expression=Val{false})[2]
    eval_∇²ₓL = (; shape=size(∇²ₓL), rows=∇²ₓL_rows, cols=∇²ₓL_cols, vals=eval_∇²ₓL_vals!) # hessian of L

    (; n, m, eval_f, eval_g!, eval_∇ₓf!, eval_∇ₓg, eval_∇²ₓL)
end

function create_ipopt(prob; tol=1e-6, max_iter=1000, verbosity=5)
    x = zeros(prob.n)
    x_l = [-Inf; 0]
    x_u = fill(Inf, prob.n)
    g_l = fill(0.0, prob.m)
    g_u = fill(Inf, prob.m)
    nele_jac_g = length(prob.eval_∇ₓg.rows)
    nele_hess = length(prob.eval_∇²ₓL.rows)

    function eval_f(x::Vector{Float64})
        prob.eval_f(x)
    end

    function eval_g(x::Vector{Float64}, g::Vector{Float64})
        prob.eval_g!(g, x)
    end

    function grad_f(x::Vector{Float64}, grad_f::Vector{Float64})
        prob.eval_∇ₓf!(grad_f, x)
    end

    function jac_g(x::Vector{Float64}, rows::Vector{Int32}, cols::Vector{Int32}, vals::Union{Nothing,Vector{Float64}})
        if vals === nothing
            rows .= prob.eval_∇ₓg.rows
            cols .= prob.eval_∇ₓg.cols
        else
            prob.eval_∇ₓg.vals(vals, x)
        end
    end

    function hess(
        x::Vector{Float64},
        rows::Vector{Int32},
        cols::Vector{Int32},
        obj_factor::Float64,
        λ::Vector{Float64},
        values::Union{Nothing,Vector{Float64}},
    )
        if values === nothing
            rows .= prob.eval_∇²ₓL.rows
            cols .= prob.eval_∇²ₓL.cols
        else
            prob.eval_∇²ₓL.vals(values, x, obj_factor, λ)
        end
    end

    ipopt_prob = Ipopt.CreateIpoptProblem(
        prob.n,
        x_l,
        x_u,
        prob.m,
        g_l,
        g_u,
        nele_jac_g,
        nele_hess,
        eval_f,
        eval_g,
        grad_f,
        jac_g,
        hess
    )
    Ipopt.AddIpoptNumOption(ipopt_prob, "tol", tol)
    Ipopt.AddIpoptIntOption(ipopt_prob, "max_iter", max_iter)
    Ipopt.AddIpoptIntOption(ipopt_prob, "print_level", verbosity)
    Ipopt.AddIpoptStrOption(ipopt_prob, "hsllib", HSL_jll.libhsl_path)
    Ipopt.AddIpoptStrOption(ipopt_prob, "linear_solver", "ma97")

	function solve(x_init=zeros(n))
		ipopt_prob.x = x_init
		#Main.@infiltrate
		solvestat = Ipopt.IpoptSolve(ipopt_prob)
		(; ipopt_prob.x, status=solvestat)
	end

    #solvestat = Ipopt.IpoptSolve(ipopt_prob)
	solve
end