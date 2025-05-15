"""
2025-05-06 Andrew

Assume, 
F: Rⁿ¹->R, 
f: Rⁿ²->R, 
G: X->Rᵐ¹, 
g: X->Rᵐ², 
are all continuous and twice differentiable.

x = [x₁, x₂] ∈ X ⊆ Rⁿ
n = n₁ + n₂

____Bilevel Optimization Problem (BOP), Standard (Optimistic):____
min     F(x)
 x
s.t.    x₁ ∈ Rⁿ¹
        G(x) ≥ 0                                            (BOP)
        x₂ ∈ S := { x₂ : x₂ ∈ arg min   f(x₁, y) 
                                 y   
                                s.t.    y ∈ Rⁿ²
                                        g(x₁, y) ≥ 0 }           
__________________________________________________________________

Consider H ⊂ S under Slater's constraint qualifications (NOTE: ∇ₓ₂ is Λ_x₂, not ∇²ₓ):
    {    x₂ ∈ Rⁿ²: ∃λ ∈ Rᵐ²           } 
H = {    ∇ₓ₂f(x) - λᵀ ∇ₓ₂g(x)) = 0   }    
    {    g(x) ≥ 0 ⟂ λ ≥ 0            }


We can replace the complementarity constraints for all i ∈ {1, …, 2ᵐ}:
                {  x₂ ∈ Rⁿ²: ∃λ ∈ Rᵐ²         }
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
s.t.    x₁ ∈ Rⁿ¹         (BOPᵢ)
        G(x) ≥ 0
        x₂ ∈ Hᵢ
________________________________________________

x₂ ∈ Hᵢ is tricky due to the complamentarity constraints Let v := [x; λ; s], we define: h(v) ⟂ z := [x₂; λ; s]
          [ ∇ₓ₂f(x) - λᵀ ∇ₓ₂g(x))  ] = 0                        -∞ ≤ x₂ ≤ ∞    
h(v) :=   [      g(x, y) - s       ] = 0              ⟂         -∞ ≤ λ ≤ ∞  
          [        λ               ] ≥/=? 0                     0 ≤ s ≤ s_ubᵢ  (upper bounds depend on i)  

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
struct BilevelOptProb
    nₓ::Int # n₁ + n₂, length(x)
    n₁::Int # length(x₁)
    n₂::Int # length(x₂)
    m₁::Int # length(G(x))
    m₂::Int # length(g(x))
    F::Function # F(x)
    G::Function # G(x)
    f::Function # f(x)
    g::Function # g(x)
    # used for solving BOPᵢ LP feasibility and NLP, let v := [x₁; z], z := [x₂; λ; s]
    nᵥ::Int  # nₓ + m₂ + m₂, length(v)
    m::Int   # m₁ + mₕ, length(Λ) = length([G(v); h(v)])
    mₕ::Int   # n₂ + m₂ + m₂, length(h(v))
    v_l₀::Vector{Float64} # default v bounds 
    v_u₀::Vector{Float64}
    Gh_l₀::Vector{Float64} # default Gh bounds
    Gh_u₀::Vector{Float64}
    eval_F::Function     # F(v)
    eval_Gh!::Function   # [G(v); h(v)], (; shape, rows, cols, vals!(out, v))
    eval_∇ᵥF!::Function  # eval_∇ᵥF!(out, v)
    eval_∇ᵥGh            # (; shape, rows, cols, vals!(out, v))
    eval_∇²ᵥL # (; shape, rows, cols, vals!(out, v, σf, Λ = [Λ₁; Λₕ])), WARN: IPOPT convention: ∇²ᵥF(v) + Λᵀ ∇²ᵥGh(v)  
    # used for solving follower's NLP
    eval_f::Function    # eval_f(x)
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
    x₁ = @view x_sym[1:n₁] # leader variables
    x₂ = @view x_sym[n₁+1:end] # follower variables
    λ_sym = Symbolics.@variables(λ[1:m₂])[1] |> Symbolics.scalarize
    s_sym = Symbolics.@variables(s[1:m₂])[1] |> Symbolics.scalarize

    z = [x₂; λ_sym; s_sym]
    @assert(mₕ == length(z))

    v = [x₁; z]
    @assert(nᵥ == length(v))

    ∇ₓ₂f = Symbolics.gradient(f_sym, x₂)
    ∇ₓ₂g = Symbolics.sparsejacobian(g_sym, x₂)
    h = [∇ₓ₂f - ∇ₓ₂g' * λ_sym; g_sym - s_sym; λ_sym]
    @assert(mₕ == length(h))

    # h(v) ⟂ z_lb ≤ z ≤ z_ub
    # x₂, λ free,  0 ≤ s ≤ s_ubᵢ
    z_l₀ = [fill(-Inf, n₂ + m₂); zeros(m₂)]
    z_u₀ = fill(Inf, n₂ + m₂ + m₂) # default z ub
    @assert(mₕ == length(z_l₀))
    @assert(mₕ == length(z_u₀))

    # v = [x₁; z], x₁ free, z bounds computed later from follower's feasible set 
    v_l₀ = [fill(-Inf, n₁); z_l₀] # if there were any leader x₁ bounds could be added here
    v_u₀ = [fill(Inf, n₁); z_u₀]
    Gh_l₀ = [fill(0, m₁); zeros(mₕ)]
    Gh_u₀ = fill(Inf, m)

    eval_F = Symbolics.build_function(F_sym, v; expression=Val{false})

    Gh = [G_sym; h] # BOPᵢ constraints
    @assert(m == length(Gh))
    eval_Gh! = Symbolics.build_function(Gh, v; expression=Val{false})[2]

    ∇ᵥF = Symbolics.gradient(F_sym, v)
    eval_∇ᵥF! = Symbolics.build_function(∇ᵥF, v; expression=Val{false})[2]

    ∇ᵥGh = Symbolics.sparsejacobian(Gh, v)

    (∇ᵥGh_rows, ∇ᵥGh_cols, ∇ᵥGh_vals) = SparseArrays.findnz(∇ᵥGh)
    eval_∇ᵥGh_vals! = Symbolics.build_function(∇ᵥGh_vals, v; expression=Val{false})[2]
    eval_∇ᵥGh = (; shape=size(∇ᵥGh), rows=∇ᵥGh_rows, cols=∇ᵥGh_cols, vals=eval_∇ᵥGh_vals!)

    # Λ(v) would be wrong!! Λ ≠ [Λ₁_sym; z] 
    Λ = Symbolics.@variables(Λ[1:m])[1] |> Symbolics.scalarize
    obj_factor = Symbolics.@variables(σf)[1]
    L = obj_factor * F_sym + Gh' * Λ # WARN: IPOPT convention: ∇²ᵥF(v) + Λᵀ ∇²ᵥ[G(v); h(v)]
    ∇ᵥL = Symbolics.gradient(L, v)
    ∇²ᵥL = Symbolics.sparsejacobian(∇ᵥL, v)
    (∇²ᵥL_rows, ∇²ᵥL_cols, ∇²ᵥL_vals) = SparseArrays.findnz(∇²ᵥL)
    eval_∇²ᵥL_vals! = Symbolics.build_function(∇²ᵥL_vals, v, obj_factor, Λ; expression=Val{false})[2]
    eval_∇²ᵥL = (; shape=size(∇²ᵥL), rows=∇²ᵥL_rows, cols=∇²ᵥL_cols, vals=eval_∇²ᵥL_vals!) # hessian of L

    # used for solving follower's NLP
    eval_f = Symbolics.build_function(f_sym, x_sym; expression=Val{false})
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
    L_feas = Λf_sym' * Gg # WARN: IPOPT convention: λᵀ ∇²ₓGg(x) because zero cost
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
        v_l₀,
        v_u₀,
        Gh_l₀,
        Gh_u₀,
        eval_F,
        eval_Gh!,
        eval_∇ᵥF!,
        eval_∇ᵥGh,
        eval_∇²ᵥL,
        eval_f,
        eval_g!,
        eval_∇ₓ₂f!,
        eval_∇ₓ₂g,
        eval_∇²ₓ₂L_follow,
        eval_Gg!,
        eval_∇ₓGg,
        eval_∇²ₓL_feas
    )
end

function solve_bop(bop; x_init=zeros(bop.nₓ), tol=1e-6, max_iter=10)
    iter_count = 0
    is_all_J_feas = false
    is_converged = false

    x₁_init = x_init[1:bop.n₁]
    x₂_init = x_init[bop.n₁+1:bop.n₁+bop.n₂]

    x::Vector{Float64} = zeros(bop.nₓ)
    λ = zeros(bop.m₂)
    s = zeros(bop.m₂)
    v = zeros(bop.nᵥ)

    # try to solving the follower's problem for the given x_init first, but this may fail if the feasible set of the follower is empty for x₁_init
    try
        x₂, λ, s = solve_follower_nlp(bop, x₁_init; y_init=x₂_init)
        v = [x₁_init; x₂; λ; s]
    catch
        # first solve the bilevel feasibility problem, then solve the follower's problem
        x₁, x₂ = find_bilevel_feas_pt(bop; x_init)
        @warn "Failed to find follower solution for x₁_init=$(x₁_init)! Changed x_init=$([x₁; x₂])"
        x₂, λ, s = solve_follower_nlp(bop, x₁; y_init=x₂)
        v = [x₁; x₂; λ; s]
    end

    # doesn't matter too much it but could happen if try above succeeds
    #if any(bop.G(x) .≤ -tol) || any(bop.g(x) .≤ -tol)
    #    @warn("Initial x is not bilevel feasible!")
    #end
    #Main.@infiltrate
    solve_Λ_feas = setup_Λ_feas_LP(bop)
    solve_BOPᵢ_nlp = setup_BOPᵢ_NLP(bop)

    prev_v = v

    while !is_converged
        iter_count += 1
        @info "iteration $iter_count"

        # WARN: this makes it go super slow
        #x₁ = v[1:bop.n₁]
        #x₂ = v[bop.n₁+1:bop.n₁+bop.n₂]
        #x₂, λ, s = solve_follower_nlp(bop, x₁; y_init=x₂)
        #x = [x₁; x₂]
        #v = [x; λ; s]

        follow_feas_Js = find_follow_feas_ind_sets(bop, v)
        is_all_J_feas = true # is there a feasible Λ for all follower feasible Js?

        if length(follow_feas_Js) > 1
            @info "multiple feasible sets detected"
        end

        for Ji in follow_feas_Js
            Ji_bounds = convert_J_to_bounds(Ji, bop)

            is_Λ_feasible = false
            try
                _, is_Λ_feasible = solve_Λ_feas(v, Ji_bounds) # does there exist Λ_all = [Λ; Λ_v_l; Λ_v_u] for BOPᵢ given v
                #@info "solve_Λ_feas successful: $is_Λ_feasible"
            catch
                is_Λ_feasible = false
                #@info "solve_Λ_feas failed: $is_Λ_feasible"
            end
            # if for some J there's no feasible Λ, we need to update v:
            if !is_Λ_feasible
                # @info "infeasible Ji detected"
                is_all_J_feas = false
                v, _ = solve_BOPᵢ_nlp(Ji_bounds; v_init=v) # update v
                break
            end
        end

        if is_all_J_feas
            is_converged = true
            #@info "all Js are feasible"

            for Ji in follow_feas_Js
                Ji_bounds = convert_J_to_bounds(Ji, bop)
                v, _, is_BOPᵢ_solved = solve_BOPᵢ_nlp(Ji_bounds; v_init=v) # check if it's actually a minimum for all feasible regions

                if is_BOPᵢ_solved  # check if all solutions agree
                    dv = v - prev_v

                    if (LinearAlgebra.norm(dv) > tol)
                        is_converged = false
                        prev_v = v
                        break
                    end
                else
                    is_converged = false
                    break
                end
            end
        end

        if iter_count > max_iter
            @error "Reached max iterations!"
            break
        end
    end

    if is_converged
        @info "Success in $iter_count iterations"
    end
    x = v[1:bop.nₓ]
    λ = v[bop.nₓ+1:bop.nₓ+bop.m₂]
    s = v[bop.nₓ+bop.m₂+1:bop.nₓ+bop.m₂+bop.m₂]

    (; x, is_converged)
end


"""
Check viable ways to satisfy h ⟂ z (aka Λₕ)

This function constructs, j∈{1,…,mₕ}:
        { j: hⱼ > 0, l = z    } case 1  : hⱼ inactive (positive), zⱼ at lb
        { j: hⱼ = 0, l = z    } case 1/2: ambiguous
K[j] =  { j: hⱼ = 0, l < z < u} case 2  : hⱼ active (l ≠ u)
        { j: hⱼ = 0,     z = u} case 2/3: ambiguous 
        { j: hⱼ < 0,     z = u} case 3  : hⱼ inactive (negative), zⱼ at ub
        { j: hⱼ,     l = z = u} case 4  : hⱼ free, zⱼ fixed

Then we sort out ambiguities by enumerating and collect them in J[1], J[2], J[3] and J[4]
"""
function find_follow_feas_ind_sets(bop, v; tol=1e-3)
    Gh = zeros(bop.m)
    bop.eval_Gh!(Gh, v)
    h = @view Gh[bop.m₁+1:end]
    z = @view v[bop.n₁+1:end] # bop.v_inds["z"]
    z_l = @view bop.v_l₀[bop.n₁+1:end]
    z_u = @view bop.v_u₀[bop.n₁+1:end]
    K = Dict{Int,Vector{Int}}()

    # note which constraints are active
    for j in 1:bop.mₕ
        Kj = Int[]

        if isapprox(z_l[j], z_u[j]; atol=2 * tol)
            push!(Kj, 4) # case 4
        elseif tol ≤ h[j] && z[j] < z_l[j] + tol
            push!(Kj, 1) # case 1
        elseif -tol < h[j] < tol && z[j] < z_l[j] + tol
            push!(Kj, 1) # case 1 OR  
            push!(Kj, 2) # case 2
        elseif -tol < h[j] < tol && z_l[j] + tol ≤ z[j] ≤ z_u[j] - tol
            push!(Kj, 2) # case 2
        elseif -tol < h[j] < tol && z_u[j] - tol < z[j]
            push!(Kj, 2) # case 2 OR 
            push!(Kj, 3) # case 3
        elseif h[j] ≤ tol && z_u[j] - tol < z[j]
            push!(Kj, 2) # case 3
        end
        K[j] = Kj
    end
    is_sol_valid = !any(isempty.(Kj for Kj in values(K)))

    if !is_sol_valid
        #Main.@infiltrate
        throw(error("Not a valid solution!"))
    end

    #if any([length(Ki) for Ki in values(K)] .> 1)
    #    Main.@infiltrate
    #end
    # enumerate the ambiguous indexes
    Js = Vector{Dict{Int,Set{Int}}}()
    ambigu_inds = [i for i in 1:bop.mₕ if length(K[i]) > 1]
    single_inds = setdiff(1:bop.mₕ, ambigu_inds)
    It = Iterators.product([K[i] for i in ambigu_inds]...)

    for assignment in It
        J = Dict(k => Set{Int}() for k in 1:4) # number of cases = 4
        for (j, ej) in enumerate(assignment)
            push!(J[ej], ambigu_inds[j])
        end
        for j in single_inds
            push!(J[K[j][1]], j)
        end
        push!(Js, J)
    end
    Js
end

function convert_J_to_bounds(J, bop)
    h_l = zeros(bop.mₕ)
    h_u = zeros(bop.mₕ)
    z_l = zeros(bop.mₕ)
    z_u = zeros(bop.mₕ)
    z_l₀ = bop.v_l₀[bop.n₁+1:end]
    z_u₀ = bop.v_u₀[bop.n₁+1:end]

    for j in J[1] # hⱼ inactive (positive), Λₕⱼ at lb
        h_l[j] = 0.0
        h_u[j] = Inf
        z_l[j] = z_l₀[j]
        z_u[j] = z_l₀[j]
    end
    for j in J[2] # hⱼ active (l ≠ u)
        h_l[j] = 0.0
        h_u[j] = 0.0
        z_l[j] = z_l₀[j]
        z_u[j] = z_u₀[j]
    end
    for j in J[3] # hⱼ inactive (negative), Λₕⱼ at ub
        h_l[j] = -Inf
        h_u[j] = 0.0
        z_l[j] = z_u₀[j]
        z_u[j] = z_u₀[j]
    end
    for j in J[4] # hⱼ free, Λₕⱼ fixed
        h_l[j] = -Inf
        h_u[j] = Inf
        z_l[j] = z_l₀[j]
        z_u[j] = z_u₀[j]
    end

    (; h_l, h_u, z_l, z_u)
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
s.t.    y_l ≤ y ≤ y_u
        g_l ≤ g(y) ≤ g_u

This may return bilevel infeasible x₂, but we only use this to find a guess for the λ_init for the feasible index sets.
"""
function solve_follower_nlp(bop, x₁; y_init=zeros(bop.n₂), tol=1e-6, max_iter=1000, verbosity=0)
    y_l = bop.v_l₀[bop.n₁+1:bop.n₁+bop.n₂]
    y_u = bop.v_u₀[bop.n₁+1:bop.n₁+bop.n₂]
    g_l = fill(0.0, bop.m₂)
    g_u = fill(Inf, bop.m₂)
    nele_jac_g = length(bop.eval_∇ₓ₂g.rows)
    nele_hess = length(bop.eval_∇²ₓ₂L_follow.rows)

    function eval_f(y::Vector{Float64})
        x = [x₁; y]
        bop.eval_f(x)
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
        y_l,
        y_u,
        bop.m₂,
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

    ipopt_prob.x = y_init
    solvestat = Ipopt.IpoptSolve(ipopt_prob)

    if solvestat != 0
        throw(error("Failed to solve follower NLP, problem may be infeasible for given x₁"))
    end

    x₂ = ipopt_prob.x
    λ = -ipopt_prob.mult_g # convention change!!
    s = bop.g([x₁; x₂])

    (; x₂, λ, s, solvestat)
end

"""
This function interfaces the BOPᵢ NLP with IPOPT, let v = [x; λ; s], Gh(v) = [G(x); h(v)]:

min     F(x)
 v
s.t.    Gh_u ≥ Gh(x) ≥ Gh_l
        v_u ≥ v ≥ v_l
"""
function setup_BOPᵢ_NLP(bop; tol=1e-6, max_iter=1000, verbosity=0)
    # these will be overwritten wrt Js, but only the follower's Λₕ=[x₂;λ;s] and h
    v_l = bop.v_l₀
    v_u = bop.v_u₀
    Gh_l = bop.Gh_l₀
    Gh_u = bop.Gh_u₀

    nele_jac_g = length(bop.eval_∇ᵥGh.rows)
    nele_hess = length(bop.eval_∇²ᵥL.rows)

    function eval_F(v::Vector{Float64})
        bop.eval_F(v)
    end

    function eval_Gh(v::Vector{Float64}, Gh::Vector{Float64})
        bop.eval_Gh!(Gh, v)
    end

    function grad_F(v::Vector{Float64}, grad_F::Vector{Float64})
        bop.eval_∇ᵥF!(grad_F, v)
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

    function solve(Ji_bounds; v_init=zeros(n))
        v_l[bop.n₁+1:bop.n₁+bop.mₕ] = Ji_bounds.z_l
        v_u[bop.n₁+1:bop.n₁+bop.mₕ] = Ji_bounds.z_u
        Gh_l[bop.m₁+1:bop.m₁+bop.mₕ] = Ji_bounds.h_l
        Gh_u[bop.m₁+1:bop.m₁+bop.mₕ] = Ji_bounds.h_u

        ipopt_prob = Ipopt.CreateIpoptProblem(
            bop.nᵥ,
            v_l,
            v_u,
            bop.m,
            Gh_l,
            Gh_u,
            nele_jac_g,
            nele_hess,
            eval_F,
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

        success = solvestat == 0

        v = ipopt_prob.x
        Λ = -ipopt_prob.mult_g # convention change!!
        (; v, Λ, success, ipopt_prob)
    end
    solve
end

"""
min     F(x)
 v
s.t.    G(x) ≥ 0
        ∇ₓ₂f(x) - λᵀ ∇ₓ₂g(x))          x₂ free
        g(x) - s                 ⟂     λ free
        λ                           0 ≤ s ≤ s_ubᵢ

from IPOPT: ∇ᵥF - ∇ᵥGh' * Λ - ipopt_prob.mult_x_L + ipopt_prob.mult_x_U = 0        
Λ_all = [Λ; Λ_l; Λ_u]
Ghb = [G; h; v; v]
KKT conditions for BOPᵢ is:
    ∇ᵥF(x) - Λ_allᵀ ∇ᵥGhb(v) - Λ_v_lᵀ + Λ_v_uᵀ = 0  
                    Λ₁  ≥ 0
               Λ₁ᵀ G(x) = 0
                    Λₕ bounds depend on convert_J_to_bounds() 

Note G(x) ≥ 0 and h(v) ≥ 0 was redacted because it has no Λ dependency, and also Λₕᵀ h(v) because it's manually taken care of using convert_J_to_bounds(). 

This is practically an LP feasibility problem when x, λ, and s are given

This function interfaces the above LP feasibility problem with HiGHS. HiGHS solves (https://ergo-code.github.io/HiGHS/dev/):
    min     cᵀΛ + d
     Λ
    s.t.    Λ_l ≤ Λ ≤ Λ_u
            A_l ≤ A Λ ≤ A_u

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
col_lower = [0.0, 1.0] (x_l)
col_upper = [4.0, Inf] (x_u)
row_lower = [-Inf, 5.0, 6.0] (A_l)
row_upper = [7.0, 15.0, Inf] (A_u)
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

    function solve(v, Ji_bounds)
        bop.eval_Gh!(Gh, v)
        bop.eval_∇ᵥF!(∇ᵥF, v)
        bop.eval_∇ᵥGh.vals(∇ᵥGh.nzval, v)
        # v bound constraints added, Λ_all = [Λ; Λ_v_l; Λ_v_u]
        ∇ᵥGhb = [∇ᵥGh; LinearAlgebra.I(bop.nᵥ); -LinearAlgebra.I(bop.nᵥ)]

        # Define the column costs, lower bounds and upper bounds
        col_cost = zeros(bop.m + 2 * bop.nᵥ) # feasibility problem
        col_lower = [fill(-Inf, bop.m); zeros(2 * bop.nᵥ)] # Λ ≥ Λ_l = 0
        col_upper = [fill(Inf, bop.m); fill(Inf, 2 * bop.nᵥ)] #  Λ_u ≥ Λ some need to be set zero based on Ji:

        #Ji_bounds = convert_J_to_bounds(Ji, bop)
        col_lower[bop.m₁+1:bop.m₁+bop.mₕ] = Ji_bounds.z_l
        col_upper[bop.m₁+1:bop.m₁+bop.mₕ] = Ji_bounds.z_u

        offset = 0.0 # unused
        # Define the row lower bounds and upper bounds A_l = A_u = b
        row_lower = [∇ᵥF; 0]
        row_upper = [∇ᵥF; 0]
        # TODO: verify this
        row_lower[bop.n₁+1:bop.n₁+bop.mₕ] = Ji_bounds.h_l
        row_upper[bop.n₁+1:bop.n₁+bop.mₕ] = Ji_bounds.h_u

        # constraint matrix is column-wise:
        G = @view Gh[1:bop.m₁]
        A = vcat(∇ᵥGhb', sparse([G' zeros(bop.mₕ + 2 * bop.nᵥ)'])) # Λ₁ᵀ G = 0 (leader complementarity) added
        a_start = A.colptr[1:end-1] .- 1
        a_index = A.rowval .- 1
        a_value = A.nzval

        """
        debug
        """
        #solve_BOPᵢ_nlp = setup_BOPᵢ_NLP(bop)
        #v, Λ, is_BOPᵢ_solved, ipopt_prob = solve_BOPᵢ_nlp(Ji_bounds; v_init=v)
        #Gh = zeros(bop.m)
        #∇ᵥF = zeros(bop.nᵥ)
        #∇ᵥGh = sparse(bop.eval_∇ᵥGh.rows, bop.eval_∇ᵥGh.cols, zeros(Cdouble, length(bop.eval_∇ᵥGh.rows)), bop.eval_∇ᵥGh.shape[1], bop.eval_∇ᵥGh.shape[2])
        #bop.eval_Gh!(Gh, v)
        #bop.eval_∇ᵥF!(∇ᵥF, v)
        #bop.eval_∇ᵥGh.vals(∇ᵥGh.nzval, v)

        #if !all(Gh[bop.m₁+1:end] .≥ Ji_bounds.h_l .- 1e-6) || !all(Gh[bop.m₁+1:end] .≤ Ji_bounds.h_u .+ 1e-6) || !all(Λ[bop.m₁+1:end] .≥ Ji_bounds.z_l .- 1e-6) || !all(Λ[bop.m₁+1:end] .≤ Ji_bounds.z_u .+ 1e-6)
        #    @error "invalid solution"
        #end

        #∇ᵥGhb = [∇ᵥGh; LinearAlgebra.I(bop.nᵥ); -LinearAlgebra.I(bop.nᵥ)]
        #Λ_all = [Λ; ipopt_prob.mult_x_L; ipopt_prob.mult_x_U]
        #@info ∇ᵥF - ∇ᵥGhb' * Λ_all

     

        #if !all(A * Λ_all .≥ row_lower .- 1e-6) || !all(A * Λ_all .≤ row_upper .+ 1e-6) || !all(Λ_all .≥ col_lower .- 1e-6) || !all(Λ_all .≤ col_upper .+ 1e-6)
        #    Main.@infiltrate
        #    @error "you have a bug here mate"
        #end
        """
        debug end
        """

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

        is_Λ_feasible = primal_sol_status[] == HiGHS.kHighsSolutionStatusFeasible
        is_Λ_dual_feasible = dual_sol_status[] == HiGHS.kHighsSolutionStatusFeasible

        (is_solved, is_Λ_feasible, is_Λ_dual_feasible, col_value, col_dual, row_value, row_dual)
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
        throw(error("Failed finding a bilevel feasible point, the problem may be bilevel infeasible! Check your x_init, G(x), and g(x)."))
    end

    x₁ = ipopt_prob.x[1:bop.n₁]
    x₂ = ipopt_prob.x[bop.n₁+1:bop.n₁+bop.n₂]
    (; x₁, x₂, solvestat)
end