"""
2025-05-06 Andrew

Assume, 
F: X->R, 
f: X->R, 
G: X->Rᵐ¹, 
g: X->Rᵐ², 
are all continuous and twice differentiable.

x = [x₁, x₂] ∈ X ⊆ Rⁿ
n = n₁ + n₂

____Bilevel Optimization Problem (BOP), Standard (Optimistic):____
min     F(x)
 x
s.t.    x₁ ∈ X₁
        G(x) ≥ 0                                            (BOP)
        x₂ ∈ S := { x₂ : x₂ ∈ arg min   f(x₁, y) 
                                 y   
                                s.t.    y ∈ Rⁿ²
                                        g(x₁, y) ≥ 0 }           
__________________________________________________________________

Consider H ⊂ S under Slater's constraint qualifications (NOTE: ∇ₓ₂ is Λ_x₂, not ∇²ₓ):
    {    x₂ ∈ X₂: ∃λ ∈ Rᵐ²           } 
H = {    ∇ₓ₂f(x) - λᵀ ∇ₓ₂g(x)) = 0   }    
    {    g(x) ≥ 0 ⟂ λ ≥ 0            }


We can replace the complementarity constraints for all i ∈ {1, …, 2ᵐ}:
                {  x₂ ∈ X₂: ∃λ ∈ Rᵐ²         }
H = ∪ Hᵢ :=     { ∇ₓ₂f(x) - λᵀ ∇ₓ₂g(x)) = 0  }   
    i           { gⱼ(x) ≥ 0, λⱼ = 0, j ∈ J⁻ᵢ  }           
                { gⱼ(x) = 0, λⱼ ≥ 0, j ∈ J⁺ᵢ  }, 
where for each i:
J⁻ᵢ ∪ J⁺ᵢ = {1, …, 2ᵐ}
J⁻ᵢ ∩ J⁺ᵢ = ∅
J⁺ᵢ ≠ J⁺ₖ, i ≠ k        

**Theorem**: x* is a local optimimum of BOP if and only if x* is a local optimium of BOPᵢ for all i: x* ∈ Hᵢ         

____Sub Bilevel Optimization Problem (BOPᵢ):____                                      
min     F(x)
x,λ,s
s.t.    x₁ ∈ X₁         (BOPᵢ)
        G(x) ≥ 0
        x₂ ∈ Hᵢ
________________________________________________

Let's define                                                      Λₕ := [x; λ; s] 
                [ ∇ₓ₂f(x) - λᵀ ∇ₓ₂g(x))  ] = 0                        -∞ ≤ x₂ ≤ ∞    
h(x, λ, s) :=   [      g(x, y) - s       ] = 0              ⟂         -∞ ≤ λ ≤ ∞  
                [        λ               ] ≥/=? 0                      0 ≤ s ≤ s_ubᵢ  (upper bounds depend on i)  

So now we can replace x₂ ∈ Hᵢ with the following constraints for BOPᵢ:
h(x, λ, s) ≥ 0

KKT conditions for BOPᵢ is:
    ∇ₓF(x) - Λ₁ᵀ ∇ₓG(x) - Λₕᵀ ∇ₓhᵢ(x, λ, s) = 0  
                   G(x) ≥ 0                      ⟂ Λ₁ ≥ 0
  hᵢ_ub(i) ≥ h(x, λ, s) ≥ 0                      ⟂ Λₕ ≥ 0

The tricky part is dealing with hᵢ_ub:
hᵢ_ub(i) :=  [  if i in J⁺ᵢ => 0, else inf ]

Note that this is practically an LP feasibility problem when x and λ are given, let Λ = [Λ₁; Λₕ]: A Λ = b, Λ ≥ 0

____Solution procedure:____
0. Create BilevelOptProb (symbolics in setup_bop)
1. Initialize by finding some x, λ (solve_follower_nlp())
2. while True:
    2.1. Find all i such that x ∈ Hᵢ (find_feas_index_sets())
    2.2. Construct Hᵢ sets to get BOPᵢ
    2.3. Check if x, λ can solve BOPᵢ for all i: x ∈ Hᵢ (setup_Λ_feas_solver()):
        If True:
            2.3.1. (Optional*) Check if x, λ is a minimum
            2.3.2. Success
        Otherwise: Solve for x and λ for the violating i (setup_leader_nlp()), go to 2.1.

*After we find a feasible Λ, we could check if (x, λ, Λ) is indeed a minimum of (BOPᵢ) by writing the sufficient conditions or using an NLP solver again.
___________________________
"""
struct BilevelOptProb # do I need to parematrize by type of functions?
    nₓ::Int # n₁ + n₂, length(x)
    n₁::Int # length(x₁)
    n₂::Int # length(x₂)
    m₁::Int # length(G(x))
    m₂::Int # length(g(x))
    F::Function # F(x)
    G::Function # G(x)
    f::Function # f(x)
    g::Function # g(x)
    # used for solving BOPᵢ LP feasibility and NLP, let v = [x, λ, s]
    nᵥ::Int  # nₓ + m₂ + m₂, length(v)
    m::Int   # m₁ + mₕ, length(Λ) = length([G(v); h(v)])
    mₕ::Int   # n₂ + m₂ + m₂, length(h(v))
    eval_h!::Function    # eval_h!(out, v)
    eval_∇ᵥh             # (; shape, rows, cols, vals!(out, v))
    eval_Gh!::Function   # [G(v); h(v)], (; shape, rows, cols, vals!(out, v))
    eval_∇ᵥF!::Function  # eval_∇ᵥF!(out, v)
    eval_∇ᵥGh            # (; shape, rows, cols, vals!(out, v))
    eval_∇²ᵥL # (; shape, rows, cols, vals!(out, v, σf, Λ = [Λ₁; Λₕ])), WARN: IPOPT convention: ∇²ᵥF(v) + Λᵀ ∇²ᵥGh(v)  
    # used for solving follower's NLP
    eval_g!::Function    # eval_g!(out, x)
    eval_∇ₓ₂f!::Function # eval_∇ₓ₂f!(out, x)
    eval_∇ₓ₂g            # (; shape, rows, cols, vals!(out, x))
    eval_∇²ₓ₂L_follow     # (; shape, rows, cols, vals!(out, v, σf, λ)), WARN: IPOPT convention
    # used for checking bilevel feasibility
    eval_Gg!::Function   #  eval_Gg!(out, x)
    eval_∇ₓGg            # (; shape, rows, cols, vals!(out, x))
    eval_∇²ₓL_feas       # WARN: IPOPT convention
end

"""
Construct BilevelOptProb
"""
function construct_bop(n₁, n₂, F, G, f, g)
    nₓ = n₁ + n₂ # length(x)

    x_dummy = zeros(nₓ)
    m₁ = length(G(x_dummy))
    m₂ = length(g(x_dummy))

    # used for solving BOPᵢ LP feasibility and NLP, let v = [x, λ, s] ∈ Rⁿᵛ
    nᵥ = nₓ + m₂ + m₂ # length(v)
    mₕ = n₂ + m₂ + m₂ # length(h(v))
    m = m₁ + mₕ       # length([G(v); h(v)])

    x_sym = Symbolics.@variables(x[1:nₓ])[1] |> Symbolics.scalarize
    F_sym = F(x_sym)
    G_sym = G(x_sym)
    f_sym = f(x_sym)
    g_sym = g(x_sym)
    x₂ = @view x_sym[n₁+1:end] # follower's variables
    λ_sym = Symbolics.@variables(λ[1:m₂])[1] |> Symbolics.scalarize
    s_sym = Symbolics.@variables(s[1:m₂])[1] |> Symbolics.scalarize
    v = [x_sym; λ_sym; s_sym]
    @assert(nᵥ == length(v))

    #Main.@infiltrate
    ∇ₓ₂f = Symbolics.gradient(f_sym, x₂)
    ∇ₓ₂g = Symbolics.sparsejacobian(g_sym, x₂)
    h = [∇ₓ₂f - ∇ₓ₂g' * λ_sym; g_sym - s_sym; λ_sym]
    @assert(mₕ == length(h))
    eval_h! = Symbolics.build_function(h, v; expression=Val{false})[2]

    ∇ᵥh = Symbolics.sparsejacobian(h, v)
    (∇ᵥh_rows, ∇ᵥh_cols, ∇ᵥh_vals) = SparseArrays.findnz(∇ᵥh)
    ∇ᵥh_vals! = Symbolics.build_function(∇ᵥh_vals, v; expression=Val{false})[2]
    eval_∇ᵥh = (; shape=size(∇ᵥh), rows=∇ᵥh_rows, cols=∇ᵥh_cols, vals=∇ᵥh_vals!)

    Gh = [G_sym; h] # BOPᵢ constraints
    @assert(m == length(Gh))
    eval_Gh! = Symbolics.build_function(Gh, v; expression=Val{false})[2]

    ∇ᵥF = Symbolics.gradient(F_sym, v)
    eval_∇ᵥF! = Symbolics.build_function(∇ᵥF, v; expression=Val{false})[2]

    ∇ᵥGh = Symbolics.sparsejacobian(Gh, v)
    (∇ᵥGh_rows, ∇ᵥGh_cols, ∇ᵥGh_vals) = SparseArrays.findnz(∇ᵥGh)
    eval_∇ᵥGh_vals! = Symbolics.build_function(∇ᵥGh_vals, v; expression=Val{false})[2]
    eval_∇ᵥGh = (; shape=size(∇ᵥGh), rows=∇ᵥGh_rows, cols=∇ᵥGh_cols, vals=eval_∇ᵥGh_vals!)

    Λ_sym = Symbolics.@variables(Λ[1:m])[1] |> Symbolics.scalarize
    obj_factor = Symbolics.@variables(σf)[1]
    L = obj_factor * f_sym + Gh' * Λ_sym # WARN: IPOPT convention: ∇²ᵥF(v) + Λᵀ ∇²ᵥ[G(v); h(v)]
    ∇ᵥL = Symbolics.gradient(L, v)
    ∇²ᵥL = Symbolics.sparsejacobian(∇ᵥL, v)
    (∇²ᵥL_rows, ∇²ᵥL_cols, ∇²ᵥL_vals) = SparseArrays.findnz(∇²ᵥL)
    eval_∇²ᵥL_vals! = Symbolics.build_function(∇²ᵥL_vals, v, obj_factor, Λ_sym; expression=Val{false})[2]
    eval_∇²ᵥL = (; shape=size(∇²ᵥL), rows=∇²ᵥL_rows, cols=∇²ᵥL_cols, vals=eval_∇²ᵥL_vals!) # hessian of L

    # used for solving follower's NLP
    eval_g! = Symbolics.build_function(g_sym, x_sym; expression=Val{false})[2]
    ∇ₓ₂f = Symbolics.gradient(f_sym, x₂)
    eval_∇ₓ₂f! = Symbolics.build_function(∇ₓ₂f, x_sym; expression=Val{false})[2]

    ∇ₓ₂g = Symbolics.sparsejacobian(g_sym, x₂)
    (∇ₓ₂g_rows, ∇ₓ₂g_cols, ∇ₓ₂g_vals) = SparseArrays.findnz(∇ₓ₂g)
    eval_∇ₓ₂g_vals! = Symbolics.build_function(∇ₓ₂g_vals, x_sym; expression=Val{false})[2]
    eval_∇ₓ₂g = (; shape=size(∇ₓ₂g), rows=∇ₓ₂g_rows, cols=∇ₓ₂g_cols, vals=eval_∇ₓ₂g_vals!)

    L_follow = obj_factor * f_sym + g_sym' * λ_sym # WARN: IPOPT convention: ∇²ₓ₂f(x) + λᵀ ∇²ₓ₂ g(x)
    ∇ₓ₂L_follow = Symbolics.gradient(L_follow, x₂)
    ∇²ₓ₂L_follow = Symbolics.sparsejacobian(∇ₓ₂L_follow, x₂)
    (∇²ₓ₂L_follow_rows, ∇²ₓ₂L_follow_cols, ∇²ₓ₂L_follow_vals) = SparseArrays.findnz(∇²ₓ₂L_follow)
    eval_∇²ₓ₂L_follow_vals! = Symbolics.build_function(∇²ₓ₂L_follow_vals, x_sym, obj_factor, λ_sym; expression=Val{false})[2]
    eval_∇²ₓ₂L_follow = (; shape=size(∇²ₓ₂L_follow), rows=∇²ₓ₂L_follow_rows, cols=∇²ₓ₂L_follow_cols, vals=eval_∇²ₓ₂L_follow_vals!) # hessian of L_follow

    # used for checking bilevel feasibility
    Gg = [G_sym; g_sym]
    eval_Gg! = Symbolics.build_function(Gg, x_sym; expression=Val{false})[2]

    ∇ₓGg = Symbolics.sparsejacobian(Gg, x_sym)
    (∇ₓGg_rows, ∇ₓGg_cols, ∇ₓGg_vals) = SparseArrays.findnz(∇ₓGg)
    eval_∇ₓGg_vals! = Symbolics.build_function(∇ₓGg_vals, x_sym; expression=Val{false})[2]
    eval_∇ₓGg = (; shape=size(∇ₓGg), rows=∇ₓGg_rows, cols=∇ₓGg_cols, vals=eval_∇ₓGg_vals!)

    Λf_sym = Symbolics.@variables(Λf[1:m₁+m₂])[1] |> Symbolics.scalarize
    L_feas = Λf_sym' * Gg # WARN: IPOPT convention: λᵀ ∇²ₓ₂ g(x) because no cost
    ∇ₓL_feas = Symbolics.gradient(L_feas, x_sym)
    ∇²ₓL_feas = Symbolics.sparsejacobian(∇ₓL_feas, x_sym)
    (∇²ₓL_feas_rows, ∇²ₓL_feas_cols, ∇²ₓL_feas_vals) = SparseArrays.findnz(∇²ₓL_feas)
    eval_∇²ₓL_feas_vals! = Symbolics.build_function(∇²ₓL_feas_vals, x_sym, Λf_sym; expression=Val{false})[2]
    eval_∇²ₓL_feas = (; shape=size(∇²ₓL_feas), rows=∇²ₓL_feas_rows, cols=∇²ₓL_feas_cols, vals=eval_∇²ₓL_feas_vals!) # hessian of L_feas

    BilevelOptProb(
        nₓ,
        n₁,
        n₂,
        m₁,
        m₂,
        F,
        G,
        f,
        g,
        nᵥ,
        m,
        mₕ,
        eval_h!,
        eval_∇ᵥh,
        eval_Gh!,
        eval_∇ᵥF!,
        eval_∇ᵥGh,
        eval_∇²ᵥL,
        eval_g!,
        eval_∇ₓ₂f!,
        eval_∇ₓ₂g,
        eval_∇²ₓ₂L_follow,
        eval_Gg!,
        eval_∇ₓGg,
        eval_∇²ₓL_feas
    )
end

function solve_bop(bop; x₁_init=zeros(bop.n₁), x₂_init=zeros(bop.n₂), tol=1e-6, max_iter=100)
    iters = 0
    converged = false
    success = false

    x::Vector{Float64} = zeros(bop.nₓ)
    v = zeros(bop.nᵥ)
    x_init = [x₁_init; x₂_init]

    # try to solving the follower's problem for the given x_init first, but this may fail if the feasible set of the follower is empty for x₁_init
    try
        x₂, λ, s = solve_follower_nlp(bop, x₁_init; y_init=x₂_init)
        x = [x₁_init; x₂]
        v = [x; λ; s]
    catch
        # first solve the bilevel feasibility problem, then solve the follower's problem
        #@info "x₁_init is not feasible, solving the bilevel feasibility problem"
        x₁, x₂ = find_bilevel_feas_pt(bop; x_init)
        #@info "x_init set to $([x₁; x₂])"
        x₂, λ, s = solve_follower_nlp(bop, x₁; y_init=x₂)
        x = [x₁; x₂]
        v = [x; λ; s]
    end

    # doesn't matter much it but could happen
    if any(bop.G(x) .≤ -tol) || any(bop.g(x) .≤ -tol)
        #@warn("starting x is not feasible")
    end

    solve_Λ_feas = setup_Λ_feas_LP(bop)
    solve_BOPᵢ_nlp = setup_BOPᵢ_NLP(bop)
    s = zeros(bop.m₂)
    prev_v = v

    while !converged
        iters += 1
        converged = true
        # TODO: this can make it go faster... or slower
        #x₁ = v[1:bop.n₁]
        #x₂ = v[bop.n₁+1:bop.n₁+bop.n₂]
        #x₂, λ, s = solve_follower_nlp(bop, x₁; y_init=x₂)
        #x = [x₁; x₂]
        #v = [x; λ; s]
        feasible_J = find_feas_index_sets(bop, v)

        for Ji in feasible_J
            # TODO 2025-05-11 FIX solve_Λ_feas
            sol = solve_Λ_feas(v, Ji)
            if sol.is_solved @warn "I actually work sometimes" end 

            if !sol.is_solved
                v, Λ = solve_BOPᵢ_nlp(Ji; v_init=v)
            end

            dv = v - prev_v

            if (dv' * dv > tol)
                converged = false
                prev_v = v
                break
            end
        end

        if converged
            success = true
            break
        end

        if iters > max_iter
            @warn "reached max iterations!"
            break
        end
    end

    if success
        x = v[1:bop.nₓ]
        λ = v[bop.nₓ+1:bop.nₓ+bop.m₂]
        s = v[bop.nₓ+bop.m₂+1:bop.nₓ+bop.m₂+bop.m₂]
    end
    x
end


"""
Check viable ways to satisfy h ⟂ Λₕ  

This function constructs:
K[1] = { j : hⱼ(v) > 0, Λₕⱼ = 0 } (case 1: inactive, J-ᵢ)
K[2] = { j : hⱼ(v) = 0, Λₕⱼ = 0 } (case 1 or 2: ambiguous)
K[3] = { j : hⱼ(v) = 0, Λₕⱼ > 0 } (case 2: active, J+ᵢ) 

Then expands ambiguous cases and creates J from K
e.g.
K = Dict{Int64, Vector{Int64}} with 3 entries:
    1 => [1]
    2 => [1, 2] <- ambiguous (multiples)
    3 => [2]

J = resolve_ambiguous_index_sets(K)
2-element Vector{Dict{Int64, Set{Int64}}}:
    Dict(2 => Set([3]), 1 => Set([2, 1]))
    Dict(2 => Set([2, 3]), 1 => Set([1]))
"""
function find_feas_index_sets(bop, v; tol=1e-3)
    h = zeros(bop.mₕ)
    bop.eval_h!(h, v)

    Λₕ = v[bop.n₁+1:end]

    # which indexes could go into J[1] and/or J[2]
    K = Dict{Int,Vector{Int}}()

    # note which constraints are active
    for j in 1:bop.mₕ
        Kj = Int[]

        if h[j] > tol && tol ≥ Λₕ[j] # case 1 inactive
            push!(Kj, 1) # J[1] candidate
        elseif tol ≥ h[j] && tol ≥ Λₕ[j]  # case 1, 2
            push!(Kj, 1) # J[1] candidate
            push!(Kj, 2) # J[2] candidate
        elseif tol ≥ h[j] && Λₕ[j] > tol # case 2 active
            push!(Kj, 2) # J[2] candidate
        end
        K[j] = Kj
    end
    valid_solution = !any(isempty.(Kj for Kj in values(K)))

    if !valid_solution
        throw(error("not a valid solution"))
    end

    J = Vector{Dict{Int,Set{Int}}}()
    K_len = length(K)
    multiples = [i for i in 1:K_len if length(K[i]) > 1]
    singles = setdiff(1:K_len, multiples)
    It = Iterators.product([K[i] for i in multiples]...)

    for assignment in It
        JJ = Dict(j => Set{Int}() for j in 1:2)
        for (e, ej) in enumerate(assignment)
            push!(JJ[ej], multiples[e])
        end
        for e in singles
            push!(JJ[K[e][1]], e)
        end
        push!(J, JJ)
    end
    J
end

"""
This function interfaces the follower NLP with IPOPT:
min   f(x₁, y) 
 y    
s.t.   y ∈ Rⁿ²
       g(x₁, y) ≥ 0

IPOPT (https://coin-or.github.io/Ipopt/) solves:
min     f(y)
 y
s.t.    yₗ ≤ y ≤ yᵤ
        gₗ ≤ g(y) ≤ gᵤ

This may return bilevel infeasible x₂, but we only use this to find a guess for the λ_init for the feasible index sets.
"""
function solve_follower_nlp(bop, x₁; y_init=zeros(bop.n₂), tol=1e-6, max_iter=1000, verbosity=0)
    y_lb = fill(-Inf, bop.n₂)
    y_ub = fill(Inf, bop.n₂)
    g_lb = fill(0.0, bop.m₂)
    g_ub = fill(Inf, bop.m₂)
    nele_jac_g = length(bop.eval_∇ₓ₂g.rows)
    nele_hess = length(bop.eval_∇²ₓ₂L_follow.rows)

    function eval_f(y::Vector{Float64})
        x = [x₁; y]
        bop.f(x)
    end

    function eval_g(y::Vector{Float64}, g::Vector{Float64})
        x = [x₁; y]
        bop.eval_g!(g, x)
    end

    function grad_f(y::Vector{Float64}, grad_f::Vector{Float64})
        x = [x₁; y]
        bop.eval_∇ₓ₂f!(grad_f, x)
    end

    function jac_g(y::Vector{Float64}, rows::Vector{Int32}, cols::Vector{Int32}, vals::Union{Nothing,Vector{Float64}})
        if vals === nothing
            rows .= bop.eval_∇ₓ₂g.rows
            cols .= bop.eval_∇ₓ₂g.cols
        else
            x = [x₁; y]
            bop.eval_∇ₓ₂g.vals(vals, x)
        end
    end

    function hess(
        y::Vector{Float64},
        rows::Vector{Int32},
        cols::Vector{Int32},
        obj_factor::Float64,
        λ::Vector{Float64},
        values::Union{Nothing,Vector{Float64}},
    )
        if values === nothing
            rows .= bop.eval_∇²ₓ₂L_follow.rows
            cols .= bop.eval_∇²ₓ₂L_follow.cols
        else
            x = [x₁; y]
            bop.eval_∇²ₓ₂L_follow.vals(values, x, obj_factor, λ)

        end
    end

    ipopt_prob = Ipopt.CreateIpoptProblem(
        bop.n₂,
        y_lb,
        y_ub,
        bop.m₂,
        g_lb,
        g_ub,
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

    ipopt_prob.x = y_init
    solvestat = Ipopt.IpoptSolve(ipopt_prob)

    if solvestat != 0
        throw(error("Failed to solve follower NLP, problem may be infeasible for given x₁"))
    end

    x₂ = ipopt_prob.x
    λ = -ipopt_prob.mult_g # convention change!!
    s = -bop.g([x₁; x₂])
    (; x₂, λ, s, solvestat)
end

"""
This function interfaces the BOPᵢ NLP with IPOPT, let v = [x; λ; s]:

min     F(x)
 v
s.t.    G(x) ≥ 0
        ∇ₓ₂f(x) - λᵀ ∇ₓ₂g(x)) = 0
        g(x) - s = 0
        λ_ubᵢ ≥ λ ≥ 0
        s_ubᵢ ≥ s ≥ 0
"""
function setup_BOPᵢ_NLP(bop; tol=1e-6, max_iter=1000, verbosity=0)
    v_lb = [fill(-Inf, bop.nₓ); zeros(bop.m₂ + bop.m₂)]
    v_ub = fill(Inf, bop.nᵥ) # will be overwritten wrt Js
    Gh_lb = fill(0.0, bop.m)
    Gh_ub = fill(Inf, bop.m) # will be overwritten wrt Js

    nele_jac_g = length(bop.eval_∇ᵥGh.rows)
    nele_hess = length(bop.eval_∇²ᵥL.rows)

    function eval_f(v::Vector{Float64})
        x = @view v[1:bop.nₓ]
        bop.F(x)
    end

    function eval_Gh(v::Vector{Float64}, Gh::Vector{Float64})
        bop.eval_Gh!(Gh, v)
    end

    function grad_F(v::Vector{Float64}, grad_F::Vector{Float64})
        x = @view v[1:bop.nₓ]
        bop.eval_∇ᵥF!(grad_F, x)
    end

    function jac_Gh(v::Vector{Float64}, rows::Vector{Int32}, cols::Vector{Int32}, vals::Union{Nothing,Vector{Float64}})
        if vals === nothing
            rows .= bop.eval_∇ᵥGh.rows
            cols .= bop.eval_∇ᵥGh.cols
        else
            bop.eval_∇ᵥGh.vals(vals, v)
        end
    end

    function hess(
        v::Vector{Float64},
        rows::Vector{Int32},
        cols::Vector{Int32},
        obj_factor::Float64,
        Λ::Vector{Float64},
        values::Union{Nothing,Vector{Float64}},
    )
        if values === nothing
            rows .= bop.eval_∇²ᵥL.rows
            cols .= bop.eval_∇²ᵥL.cols
        else
            bop.eval_∇²ᵥL.vals(values, v, obj_factor, Λ)
        end
    end

    function solve(Ji; v_init=zeros(n))
        for j in Ji[1] # inactive Λₕ = [x₂; λ; s]ⱼ = 0, j ∈ j∈J⁻ᵢ  
            v_ub[bop.n₁+j] = 0
            #  Gh_ul[bop.m₁+j] = tol?
        end

        for j in Ji[2] # active hⱼ = 0, j ∈ j∈J⁺ᵢ  
            #  v_lb[bop.n₁+j] = tol?
            Gh_ub[bop.m₁+j] = 0
        end

        ipopt_prob = Ipopt.CreateIpoptProblem(
            bop.nᵥ,
            v_lb,
            v_ub,
            bop.m,
            Gh_lb,
            Gh_ub,
            nele_jac_g,
            nele_hess,
            eval_f,
            eval_Gh,
            grad_F,
            jac_Gh,
            hess
        )

        Ipopt.AddIpoptNumOption(ipopt_prob, "tol", tol)
        Ipopt.AddIpoptIntOption(ipopt_prob, "max_iter", max_iter)
        Ipopt.AddIpoptIntOption(ipopt_prob, "print_level", verbosity)

        ipopt_prob.x = v_init
        solvestat = Ipopt.IpoptSolve(ipopt_prob)

        if solvestat != 0
            throw(error("Failed to solve BOPᵢ NLP!"))
        end

        v = ipopt_prob.x
        Λ = -ipopt_prob.mult_g # convention change!!
        (; v, Λ, solvestat)
    end
    solve
end


"""
KKT conditions for BOPᵢ is:
    ∇ᵥF(x) - Λ₁ᵀ ∇ᵥG(x) - Λₕᵀ ∇ᵥh(v) = 0  
                   G(x) ≥ 0  (no Λ dependency)              
                     Λ₁ ≥ 0
               Λ₁ᵀ G(x) ≥ 0
                   h(v) ≥ 0  (no Λ dependency)
              0 ≥ hⱼ(v)    , j∈J⁻ᵢ  (no Λ dependency)             
                     Λₕⱼ ≥ 0
                 0 ≥ Λₕⱼ    , j∈J⁺ᵢ  
                Λₕᵀ h(v) ≥ 0  (already taken care of)

Note that this is practically an LP feasibility problem when x and λ are given, let Λ = [Λ₁; Λₕ]: A Λ = b, Λ ≥ 0 

This function interfaces the above LP feasibility problem with HiGHS. HiGHS solves (https://ergo-code.github.io/HiGHS/dev/):
    min     cᵀΛ + d
     Λ
    s.t.    Λₗ ≤ Λ ≤ Λᵤ
            Aₗ ≤ A Λ ≤ Aᵤ
 
So basically: A = [∇ᵥGhᵀ; G(x)ᵀ 0ᵀ], b = [∇ᵥF(x); 0], Aₗ = Aᵤ = b, Λₗ = 0, Λ₁ᵤ = Inf, Λₕᵤⱼ = 0 j∈J⁻ᵢ, Λₕᵤⱼ = Inf j∈J⁺ᵢ      

Here's a HiGHS encoding example because the API is kinda obscure:
______________________________
Min     x₀ +  x₁ + 3
 x
s.t.    0 ≤ x₀ ≤ 4
        1 ≤ x₁   
                  x₁ ≤ 7
        5 ≤ x₀ + 2x₁ ≤ 15
        6 ≤ 3x₀ + 2x₁
______________________________
This would be encoded like this:
___________________________________
col_cost = [1.0, 1.0]
offset = 3.;
col_lower = [0.0, 1.0] (xₗ)
col_upper = [4.0, Inf] (xᵤ)
row_lower = [-Inf, 5.0, 6.0] (Aₗ)
row_upper = [7.0, 15.0, Inf] (Aᵤ)
a_start = [0, 2] (column-wise, column start index, the first index is always zero)
a_index = [1, 2, 0, 1, 2]
a_value = [1.0, 3.0, 1.0, 2.0, 2.0]
___________________________________
"""
function setup_Λ_feas_LP(bop; primal_feas_tol=1e-6, zero_tol=1e-3, verbosity=0)
    model = HiGHS.Highs_create()
    HiGHS.Highs_setDoubleOptionValue(model, "primal_feasibility_tolerance",
        primal_feas_tol)
    HiGHS.Highs_setBoolOptionValue(model, "output_flag", verbosity)

    Gh = zeros(bop.m)
    ∇ᵥF = zeros(bop.nᵥ)
    ∇ᵥGh = sparse(bop.eval_∇ᵥGh.rows, bop.eval_∇ᵥGh.cols, zeros(Cdouble, length(bop.eval_∇ᵥGh.rows)), bop.eval_∇ᵥGh.shape[1], bop.eval_∇ᵥGh.shape[2])
    h = zeros(bop.mₕ)
    ∇ᵥh = sparse(bop.eval_∇ᵥh.rows, bop.eval_∇ᵥh.cols, zeros(Cdouble, length(bop.eval_∇ᵥh.rows)), bop.eval_∇ᵥh.shape[1], bop.eval_∇ᵥh.shape[2])

    function solve(v, Ji)
        bop.eval_Gh!(Gh, v)
        bop.eval_∇ᵥF!(∇ᵥF, v)
        bop.eval_∇ᵥGh.vals(∇ᵥGh.nzval, v)
        bop.eval_h!(h, v)
        bop.eval_∇ᵥh.vals(∇ᵥh.nzval, v)

        # Reminder: A = [∇ᵥGhᵀ; G(x)ᵀ 0ᵀ], b = [∇ᵥF(x); 0], Aₗ = Aᵤ = b, Λₗ = 0, Λ₁ᵤ = Inf, Λₕᵤⱼ = 0 j∈J⁻ᵢ, Λₕᵤⱼ = Inf j∈J⁺ᵢ       
        # Define the column costs, lower bounds and upper bounds
        col_cost = zeros(bop.m) # feasibility problem
        col_lower = zeros(bop.m) # Λ ≥ Λₗ = 0
        col_upper = fill(Inf, bop.m) #  Λᵤ ≥ Λ some need to be set zero based on Ji:
        for j in Ji[1] # active hⱼ > 0, in j∈J⁺ᵢ, Λₕᵤⱼ = 0 j∈J⁻ᵢ 
            col_upper[bop.m₁+j] = 0
        end

        offset = 0.0 # unused
        # Define the row lower bounds and upper bounds Aₗ = Aᵤ = b
        row_lower = [∇ᵥF; 0]
        row_upper = [∇ᵥF; 0]
        # constraint matrix is column-wise:
        G = @view Gh[1:bop.m₁]
        A = vcat(∇ᵥGh', sparse([G' zeros(bop.mₕ)']))
        a_start = A.colptr[1:end-1] .- 1
        a_index = A.rowval .- 1
        a_value = A.nzval
        
        n_col = convert(Cint, size(col_cost, 1))
        n_row = convert(Cint, size(row_lower, 1))
        n_nz = convert(Cint, size(a_index, 1))
        col_cost = convert(Array{Cdouble}, col_cost)
        col_lower = convert(Array{Cdouble}, col_lower)
        col_upper = convert(Array{Cdouble}, col_upper)
        offset = convert(Cdouble, offset)
        row_lower = convert(Array{Cdouble}, row_lower)
        row_upper = convert(Array{Cdouble}, row_upper)
        a_start = convert(Array{Cint}, a_start)
        a_index = convert(Array{Cint}, a_index)
        a_value = convert(Array{Cdouble}, a_value)

        status = HiGHS.Highs_passLp(
            model,
            n_col,
            n_row,
            n_nz,
            HiGHS.kHighsMatrixFormatColwise,
            HiGHS.kHighsObjSenseMinimize,
            offset,
            col_cost,
            col_lower,
            col_upper,
            row_lower,
            row_upper,
            a_start,
            a_index,
            a_value
        )

        if status != HiGHS.kHighsStatusOk
            throw(error("Failed passing model to HiGHS!"))
        end

        if HiGHS.Highs_run(model) != HiGHS.kHighsStatusOk
            throw(error("Failed LP solve!"))
        end

        model_status = HiGHS.Highs_getModelStatus(model)
        is_solved = model_status[] == HiGHS.kHighsModelStatusOptimal

        # The vector x is col_value, the vector of duals for the variables x is col_dual
        col_value, col_dual =
            (Array{Cdouble,1}(undef, n_col), Array{Cdouble,1}(undef, n_col))
        # The vector g is row_value, the vector of duals for the variables g is row_dual
        row_value, row_dual =
            (Array{Cdouble,1}(undef, n_row), Array{Cdouble,1}(undef, n_row))

        HiGHS.Highs_getSolution(model, col_value, col_dual, row_value, row_dual)

        primal_sol_status = Ref{Cint}(0)
        dual_sol_status = Ref{Cint}(0)
        HiGHS.Highs_getIntInfoValue(model, "primal_solution_status", primal_sol_status)
        HiGHS.Highs_getIntInfoValue(model, "dual_solution_status", dual_sol_status)

        is_primal_feasible = primal_sol_status[] == HiGHS.kHighsSolutionStatusFeasible
        is_dual_feasible = dual_sol_status[] == HiGHS.kHighsSolutionStatusFeasible

        (; is_solved, is_primal_feasible, is_dual_feasible, col_value, col_dual, row_value, row_dual)
    end

    solve
end


"""
Find x such that
    g(x) ≥ 0 
    G(x) ≥ 0 
"""
function find_bilevel_feas_pt(bop; x_init=zeros(bop.nₓ), tol=1e-6, max_iter=1000, verbosity=0)
    # no need to do anything if it's already feasible
    if all(bop.G(x_init) .≥ -tol) && all(bop.g(x_init) .≥ -tol)
        return x_init
    end

    x_lower_bound = fill(-Inf, bop.nₓ)
    x_upper_bound = fill(Inf, bop.nₓ)
    m_Gg = bop.m₁ + bop.m₂
    Gg_lower_bound = fill(0.0, m_Gg)
    Gg_upper_bound = fill(Inf, m_Gg)
    nele_jac_g = length(bop.eval_∇ₓGg.rows)
    nele_hess = length(bop.eval_∇²ₓL_feas.rows)

    # feasibility problem: there's no cost
    function eval_f(x::Vector{Float64})
        0.0
    end

    function eval_Gg(x::Vector{Float64}, g::Vector{Float64})
        bop.eval_Gg!(g, x)
    end

    function grad_f(x::Vector{Float64}, grad_f::Vector{Float64})
        grad_f .= 0.0
    end

    function jac_g(
        x::Vector{Float64},
        rows::Vector{Int32},
        cols::Vector{Int32},
        values::Union{Nothing,Vector{Float64}},
    )
        if values === nothing
            rows .= bop.eval_∇ₓGg.rows
            cols .= bop.eval_∇ₓGg.cols
        else
            bop.eval_∇ₓGg.vals(values, x)
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
            rows .= bop.eval_∇²ₓL_feas.rows
            cols .= bop.eval_∇²ₓL_feas.cols
        else
            bop.eval_∇²ₓL_feas.vals(values, x, λ)
        end
    end

    ipopt_prob = Ipopt.CreateIpoptProblem(
        bop.nₓ,
        x_lower_bound,
        x_upper_bound,
        m_Gg,
        Gg_lower_bound,
        Gg_upper_bound,
        nele_jac_g,
        nele_hess,
        eval_f,
        eval_Gg,
        grad_f,
        jac_g,
        hess
    )

    Ipopt.AddIpoptNumOption(ipopt_prob, "tol", tol)
    Ipopt.AddIpoptIntOption(ipopt_prob, "max_iter", max_iter)
    Ipopt.AddIpoptIntOption(ipopt_prob, "print_level", verbosity)

    ipopt_prob.x = x_init
    solvestat = Ipopt.IpoptSolve(ipopt_prob)

    if solvestat != 0
        throw(error("Failed initialization, follower problem may be infeasible for the given x₁!"))
    end

    x₁ = ipopt_prob.x[1:bop.n₁]
    x₂ = ipopt_prob.x[bop.n₁+1:bop.n₁+bop.n₂]
    (; x₁, x₂, solvestat)
end