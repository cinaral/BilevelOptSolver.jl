using BilevelOptSolver
using Test
using BenchmarkTools

n::Int64 = 2
nq::Int64 = 0
np::Int64 = 1

function f(x)
    x[nq+1]^2 + x[nq+2]^2
end

function g(x) # ≥ 0 
    [
        x[nq+1] + x[nq+2] - 1.0;
    ]
end

nlp, nlp_sym = BilevelOptSolver.convert_to_NLP(n, f, g; nq, np, xl=zeros(n), xu=fill(Inf, n));
mcp, mcp_sym = BilevelOptSolver.convert_NLP_to_MCP(nlp.n, nlp.m, nlp_sym.x, nlp_sym.g, nlp_sym.λ, nlp_sym.of, nlp_sym.∇ₓL, nlp_sym.q, nlp_sym.p; xl=fill(-Inf, nlp.n), xu=fill(Inf, nlp.n))

x, λ, solvestat, ipopt_prob = BilevelOptSolver.solve_NLP(nlp.n, nlp.m, nlp.xl, nlp.xu, nlp.gl, nlp.gu, nlp.f, nlp.g!, nlp.∇ₓf!, nlp.∇ₓg_rows, nlp.∇ₓg_cols, nlp.∇ₓg_vals!, nlp.∇²ₓL_rows, nlp.∇²ₓL_cols, nlp.∇²ₓL_vals!; x_init=zeros(nlp.n), λ_init=zeros(nlp.m), q=Float64[], p=Float64[], verbosity=0, tol=1e-8, max_iter=1000, print_level=0, is_debug_on=false, is_using_HSL=false)
@test isapprox(x, [0.5; 0.5])

z, status, info = BilevelOptSolver.solve_PATH(mcp.n, mcp.zl, mcp.zu, mcp.h!, mcp.∇h_rows, mcp.∇h_cols, mcp.∇h_vals!; q=zeros(nq), p=zeros(np), verbosity=0, is_debug_on=false, x_init=zeros(mcp.n), tol=1e-8, max_iter=500, is_silent=true)
@test isapprox(z, [0.5; 0.5; 1.0])


