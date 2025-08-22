"""
[IPOPT](https://coin-or.github.io/Ipopt/) solves:
```
min     f(x)
 x
s.t.    xl ≤ x ≤ xu
        gl ≤ g(x) ≤ gu
```
- L(x, λ, α) := α f(x) - λ' g(x)
- g! and ∇ₓf! are in-place
- ∇ₓg and ∇ₓₓL are in-place sparse matrices (rows, cols, vals!)
"""
function solve_NLP(n, m, xl, xu, gl, gu, f, g!, ∇ₓf!, ∇ₓg_rows, ∇ₓg_cols, ∇ₓg_vals!, ∇ₓₓL_rows, ∇ₓₓL_cols, ∇ₓₓL_vals!; x_init=zeros(n), λ_init=zeros(m), q=Float64[], p=Float64[], verbosity=0, tol=1e-6, max_iter=1000, print_level=0, is_debug_on=false, is_using_HSL=false, solver="IPOPT")

    if is_using_HSL && !haskey(ENV, "HSL_PATH")
        is_using_HSL = false
        if verbosity > 0
            print("HSL_PATH not found: Defaulting is_using_HSL = false. If you would like to use HSL, please obtain a license and download HSL_jll.jl (https://licences.stfc.ac.uk/products/Software/HSL/LibHSL), and set HSL_PATH environment variable to the extracted location.\n")
        end
    end
    n_∇ₓₓL_nz = length(∇ₓₓL_rows)
    n_∇ₓg_nz = length(∇ₓg_rows)

    if is_debug_on
        @assert(n == length(xl), "n and xl dim mismatch")
        @assert(n == length(xu), "n and xu dim mismatch")
        @assert(m == length(gl), "n and gl dim mismatch")
        @assert(m == length(gu), "m and gu dim mismatch")
        @assert(n_∇ₓg_nz == length(∇ₓg_cols), "∇ₓg row and col dim mismatch")
        @assert(n_∇ₓₓL_nz == length(∇ₓₓL_cols), "∇ₓₓL row and col dim mismatch")
    end

    qxp = [q; zeros(n); p]
    nq = length(q)
    x_inds = nq+1:nq+n

    function eval_f(x)
        qxp[x_inds] .= x
        f(qxp)
    end

    function eval_g(x, g)
        qxp[x_inds] .= x
        g!(g, qxp)
    end

    function eval_∇ₓf(x, ∇ₓf)
        qxp[x_inds] .= x
        ∇ₓf!(∇ₓf, qxp)
    end

    function eval_∇ₓg(x, rows, cols, vals)
        if vals === nothing
            rows .= ∇ₓg_rows
            cols .= ∇ₓg_cols
        else
            qxp[x_inds] .= x
            ∇ₓg_vals!(vals, qxp)
        end
    end

    # 2025-08-20 TODO: non-zero allocations from convention change
    function eval_∇ₓₓL(x, rows, cols, obj_factor, λ, vals)
        if vals === nothing
            rows .= ∇ₓₓL_rows
            cols .= ∇ₓₓL_cols
        else
            qxp[x_inds] .= x
            λ .= -λ # convention change: IPOPT expects L = obj_factor * f(x) + g(x)' λ
            ∇ₓₓL_vals!(vals, qxp, λ, obj_factor)
        end
    end

    ipopt_prob = Ipopt.CreateIpoptProblem(
        n,
        xl,
        xu,
        m,
        gl,
        gu,
        n_∇ₓg_nz,
        n_∇ₓₓL_nz,
        eval_f,
        eval_g,
        eval_∇ₓf,
        eval_∇ₓg,
        eval_∇ₓₓL
    )
    if is_using_HSL
        Ipopt.AddIpoptStrOption(ipopt_prob, "hsllib", HSL_jll.libhsl_path)
        Ipopt.AddIpoptStrOption(ipopt_prob, "linear_solver", "ma27")
    end

    Ipopt.AddIpoptNumOption(ipopt_prob, "tol", tol)
    Ipopt.AddIpoptIntOption(ipopt_prob, "max_iter", max_iter)
    Ipopt.AddIpoptIntOption(ipopt_prob, "print_level", print_level)

    ipopt_prob.x .= x_init
    ipopt_prob.mult_g .= -λ_init # convention change!!
    solvestat = Ipopt.IpoptSolve(ipopt_prob)

    x = ipopt_prob.x
    λ = -ipopt_prob.mult_g # convention change!!
    (; x, λ, solvestat, ipopt_prob)
end


"""
[PATHSolver](https://coin-or.github.io/Ipopt/) solves:
```
x s.t. at least one of the following holds: 
    1.  F_i(x) = 0, xl_i ≤ x_i ≤ xu_i
    2.  F_i(x) > 0, xl_i = x_i  
    3.  F_i(x) < 0,        x_i = xu_i
```
- J := ∇ₓF
"""
function solve_PATH(n, xl, xu, F, J_rows, J_cols, J_vals!; q=Float64[], p=Float64[], verbosity=0, is_debug_on=false, x_init=zeros(n), tol=1e-6, max_iter=500, is_silent=true)
    if n > 300 && !haskey(ENV, "PATH_LICENSE_STRING")
        if verbosity > 0
            print("PATH_LICENSE_STRING not found and problem size is too large: Please obtain a license (https://pages.cs.wisc.edu/~ferris/path/julia/LICENSE), and set PATH_LICENSE_STRING environment variable.\n")
        end
    end
    if haskey(ENV, "PATH_LICENSE_STRING")
        PATHSolver.c_api_License_SetString(ENV["PATH_LICENSE_STRING"])
    end

    n_nz_J = length(J_rows)

    if is_debug_on
        @assert(n == length(xl), "wrong x_l dimensions")
        @assert(n == length(xu), "wrong x_u dimensions")
        @assert(n_nz_J == length(J_cols), "J_rows and J_cols dimensions do not match")
    end

    J_shape = sparse(J_rows, J_cols, ones(Cdouble, n_nz_J), n, n)
    J_col = J_shape.colptr[1:end-1]
    J_len = diff(J_shape.colptr)
    J_row = J_shape.rowval

    # 2025-07-13 Is this needed?
    #if isempty(J_row)
    #    J_row = 0
    #end
    qxp = [q; zeros(n); p]
    nq = length(q)
    x_inds = nq+1:nq+n

    function eval_F(n, x, vals)
        vals .= 0.0
        qxp[x_inds] .= x
        F(vals, qxp)
        Cint(0)
    end

    function eval_J(n, nnz, x, col, len, row, vals)
        vals .= 0.0
        qxp[x_inds] .= x
        J_vals!(vals, qxp)
        col .= J_col
        len .= J_len
        row .= J_row
        Cint(0)
    end

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
        lemke_rank_deficiency_iterations=30, # can fix silent crashes " ** SOLVER ERROR ** Lemke: invertible basis could not be computed."
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

