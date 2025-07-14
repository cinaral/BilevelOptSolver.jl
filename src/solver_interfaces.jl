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
function setup_mcp_solve_PATH(n, xl, xu, F, J_rows, J_cols, J_vals!)
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
