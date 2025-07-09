struct BilevelOptProb
    n1::Int # length(x₁)
    n2::Int # length(x₂)
    m1::Int # length(G(x))
    m2::Int # length(g(x))
    F::Function # F(x)
    G::Function # G(x)
    f::Function # f(x)
    g::Function # g(x)
    solve_follower_nlp::Function
    solve_follower_KKT_mcp::Function
    find_bilevel_feas_pt::Function
    nv::Int # length(v) = n₁ + n₂ + m₂ + m₂ 
    mΛ::Int # length(Λ) = m₁ + mₕ + m₂
    mh::Int # length(h) = n₂ + m₂ + m₂ 
    v_inds::Dict{String,UnitRange{Int64}} # indexing v := [x₁; z]
    z_inds::Dict{String,UnitRange{Int64}} # indexing z := [x₂; λ; s]
    z_l₀::Vector{Float64} # default z bounds 
    z_u₀::Vector{Float64} # default z bounds 
    Ghs_inds::Dict{String,UnitRange{Int64}}
    Ghs_l₀::Vector{Float64} # default [G; h] bounds 
    Ghs_u₀::Vector{Float64}
    Ghs!::Function # Ghs(out, v), out := [G(v); h(v); s]
    solve_BOPᵢ_nlp::Function
    info_BOPᵢ::Function
    eval_BOPᵢ!::Function
    check_Λ_lp_feas::Function
    n_θ::Int
    θ_inds::Dict{String,UnitRange{Int64}} # indexing θ := [v; Λ; r], v := [x₁; z], z := [x₂; λ; s]
    θ_l₀::Vector{Float64} # default θ bounds
    θ_u₀::Vector{Float64}
    solve_BOPᵢ_KKT_mcp::Function
    deriv_funs
    sym_derivs
    ∇ₓ₂ₓ₂L_rows
    ∇ₓ₂ₓ₂L_cols
    ∇ₓ₂ₓ₂L
    ∇ₓ₂ₓ₂L_vals!
end

"""
Construct BilevelOptProb
```
Assume, 
    F: Rⁿ¹->R, 
    f: Rⁿ²->R, 
    G: X->Rᵐ¹, 
    g: X->Rᵐ², 
    are all continuous and twice differentiable.

Decision variables:
    x = [x₁, x₂] ∈ X ⊆ Rⁿ
        ^    ^
    Lead / Follower variables

Bilevel Optimization Problem (BOP), Standard (Optimistic):
    min     F(x)
    x
    s.t.    x₁ ∈ Rⁿ¹
            G(x) ≥ 0                                            (BOP)
            x₂ ∈ S := { x₂ : x₂ ∈ arg min   f(x₁, y) 
                                    y   
                                    s.t.    y ∈ Rⁿ²
                                            g(x₁, y) ≥ 0 }           

Consider H ⊆ S under Slater's constraint qualifications:
        { x₂ ∈ Rⁿ²: ∃λ ∈ Rᵐ²        } 
    H = { ∇ₓ₂f(x) - λᵀ ∇ₓ₂g(x)) = 0 }     (KKT conditions of follower)  
        { g(x) ≥ 0 ⟂ λ ≥ 0          }
```
**Note**: ```∇ₓ₂ is ∇_(x₂), and not ∇ₓ²```
```
We can replace the complementarity constraints by enumerating (local) manifolds i ∈ I:
                {  x₂ ∈ Rⁿ²: ∃λ ∈ Rᵐ²        }
    H = ∪ Hᵢ := { ∇ₓ₂f(x) - λᵀ ∇ₓ₂g(x)) = 0  }   
        i       { gⱼ(x) ≥ 0, λⱼ = 0, j ∈ J⁻ᵢ }   <- Inactive follower constraints           
                { gⱼ(x) = 0, λⱼ ≥ 0, j ∈ J⁺ᵢ },  <- Active follower constraints
    where for all i-th manifold, the index sets (J⁻ᵢ, J⁺ᵢ) satisfy:
        J⁻ᵢ ∪ J⁺ᵢ = {1, …, 2⁽ᵐ²⁾}
        J⁻ᵢ ∩ J⁺ᵢ = ∅
        J⁺ᵢ ≠ J⁺ₖ, i ≠ k        
```
**Theorem**: ```x* is a local optimimum of BOP if and only if x* is a local optimium of BOPᵢ for all i: x* ∈ Hᵢ```
```
Reduced Bilevel Optimization Problem:                                     
     min     F(x)
    x,λ,s
     s.t.    x₁ ∈ Rⁿ¹           (BOPᵢ)
             G(x) ≥ 0
             x₂ ∈ Hᵢ(x, λ, s)

Let v := [x₁; x₂; λ; s], and h(v) be the complement of z := [x₂; λ; s]:
                                         h(v)     ⟂       z
                                         ^:=              ^:=
                  0 = [ ∇ₓ₂f(x) - λᵀ ∇ₓ₂g(x)) ]   ⟂ -∞ ≤ [ x₂ ] ≤ ∞ (free)
                  0 = [      g(x, y) - s      ]   ⟂ -∞ ≤ [ λ  ] ≤ ∞ (free) 
                  0 ≤ [        λ              ]   ⟂  0 ≤ [ s  ] ≤ s_uᵢ  (upper bounds depend on manifold)  

If we replace x₂ ∈ Hᵢ with more constraints, we can solve BOPᵢ as an NLP:              
    min     F(x)
     v
    s.t.     x₁ ∈ Rⁿ¹            (BOPᵢ NLP)
              0 ≤ G(x) 
              0 ≤ hᵢ(v) ≤ h_uᵢ
           z_lᵢ ≤   z   ≤ z_uᵢ
    where the bounds are selected by first constructing J⁻ᵢ and J⁺ᵢ index sets.

Let Ghᵢ := [G; hᵢ], and Λ := [Λ₁; Λₕ], the KKT conditions of this BOPᵢ NLP is:
    ∃Λ ∈ R⁽ᵐ¹⁺ᵐʰ⁾, ∃Λ_l ∈ R⁽ⁿᵛ⁾, ∃Λ_u ∈ R⁽ⁿᵛ⁾:
        ∇ᵥF(x) - Λᵀ ∇ᵥGhᵢ(v) - Λ_v_lᵀ + Λ_v_uᵀ = 0                                            
                Ghᵢ(v) ≥ 0 ⟂   Λ ≥ 0              (KKT conditions of BOPᵢ NLP)
              v - v_lᵢ ≥ 0 ⟂ Λ_v_lᵀ ≥ 0                
             -v + v_uᵢ ≥ 0 ⟂ Λ_v_u ≥ 0 
```
**Note**: ```When v is given, this is an LP feasibility problem: A Λ_all = b, Λ_all ≥ 0 where Λ_all = [Λ; Λ_l; Λ_u]```    
"""
function construct_bop(n1, n2, F, G, f, g; verbosity=0)
    nx = n1 + n2 # length(x)
    x_dummy = zeros(nx)
    m1 = length(G(x_dummy))
    m2 = length(g(x_dummy))
    nv = n1 + n2 + m2 + m2 # length(v)
    mh = n2 + m2 + m2 # length(h)
    x1_sym = Symbolics.@variables(x1[1:n1])[1] |> Symbolics.scalarize
    x2_sym = Symbolics.@variables(x2[1:n2])[1] |> Symbolics.scalarize
    #λ0_sym = Symbolics.@variables(lamb0[1])[1] |> Symbolics.scalarize
    λ_sym = Symbolics.@variables(lamb[1:m2])[1] |> Symbolics.scalarize
    s_sym = Symbolics.@variables(s[1:m2])[1] |> Symbolics.scalarize
    z_sym = [x2_sym; λ_sym; s_sym] # z = [x₂; λ; s] ∈ R⁽ᵐʰ⁾
    v_sym = [x1_sym; z_sym] # v = [x₁; x₂; λ; s] ∈ R⁽ⁿᵛ⁾
    @assert(mh == length(z_sym))
    @assert(nv == length(v_sym))

    # defined for convenience
    z_inds = Dict([
        ("x2", 1:n2),
        ("λ", n2+1:n2+m2),
        ("s", n2+m2+1:n2+m2+m2),
    ])
    v_inds = Dict([
        ("x1", 1:n1),
        ("z", n1+1:n1+mh),
        ("x", 1:nx),
        ("x2", n1+1:nx),
        ("λ", nx+1:nx+m2),
        ("s", nx+m2+1:nx+2*m2)
    ])

    x_sym = [x1_sym; x2_sym]
    F_sym = F(x_sym)
    G_sym = G(x_sym)
    f_sym = f(x_sym)
    g_sym = g(x_sym)

    ∇ₓ₂f = Symbolics.gradient(f_sym, x2_sym)
    ∇ₓ₂g = Symbolics.sparsejacobian(g_sym, x2_sym)

    g_s_sym = g_sym - s_sym

    if isempty(λ_sym)
        h_sym = [∇ₓ₂f; g_s_sym; λ_sym]
    else
        h_sym = [∇ₓ₂f - ∇ₓ₂g' * λ_sym; g_s_sym; λ_sym]
    end
    @assert(mh == length(h_sym))

    # hᵢ ⟂ z_lᵢ ≤ z ≤ z_uᵢ
    # x₂, λ are free, but 0 ≤ s ≤ s_ubᵢ
    z_l₀ = fill(-Inf, mh) # default z lb
    z_u₀ = fill(Inf, mh) # default z ub
    z_l₀[z_inds["s"]] .= zeros(m2)

    g_l = fill(0.0, m2)
    g_u = fill(Inf, m2)
    solve_follower_nlp, ∇ₓ₂ₓ₂L_rows, ∇ₓ₂ₓ₂L_cols, ∇ₓ₂ₓ₂L, ∇ₓ₂ₓ₂L_vals! = setup_follower_nlp(n1, n2, m2, g_l, g_u, x_sym, λ_sym, x2_sym, f_sym, g_sym)
    solve_follower_KKT_mcp = setup_follower_KKT_mcp(n1, n2, m2, x1_sym, x2_sym, λ_sym, s_sym, f_sym, g_sym, z_inds)
    find_bilevel_feas_pt = setup_find_bile_feas_pt(n1, n2, m1, m2, x_sym, F_sym, G_sym, g_sym)

    mΛ = m1 + mh + m2 # length([G; h; s])
    Ghs_inds = Dict([ # defined for convenience
        ("G", 1:m1),
        ("h", m1+1:m1+mh),
        ("s", m1+mh+1:mΛ),
        ("∇ₓ₂L", m1+1:m1+n2), # part of h that corresponds to follower's stationary condtion
        ("g_s", m1+n2+1:m1+n2+m2),  # part of h that corresponds to follower's constraints
        ("λ", m1+n2+m2+1:m1+mh)  # part of h that corresponds to follower's dual
    ])
    Ghs_l₀ = fill(-Inf, mΛ) # default Gh lower bound (shouldn't change)
    Ghs_u₀ = fill(Inf, mΛ) # default Gh upper bound (this depends on follower solution)
    # stationary condition and constraint - slack is zero, their complements x₂ and λ are free 
    Ghs_l₀[Ghs_inds["G"]] .= 0
    Ghs_l₀[Ghs_inds["∇ₓ₂L"]] .= 0
    Ghs_u₀[Ghs_inds["∇ₓ₂L"]] .= 0
    Ghs_l₀[Ghs_inds["g_s"]] .= 0
    Ghs_u₀[Ghs_inds["g_s"]] .= 0
    Ghs_l₀[Ghs_inds["λ"]] .= 0
    Ghs_l₀[Ghs_inds["s"]] .= 0

    Ghs_sym = [G_sym; h_sym; s_sym] # BOPᵢ constraints
    @assert(mΛ == length(Ghs_sym))

    solve_BOPᵢ_nlp, info_BOPᵢ, eval_BOPᵢ!, Ghs!, ∇ᵥF!, ∇ᵥGhs_rows, ∇ᵥGhs_cols, ∇ᵥGhs_vals!, ∇ᵥGhs_shape = setup_BOPᵢ_nlp(nv, mΛ, Ghs_l₀, Ghs_u₀, v_sym, F_sym, Ghs_sym)

    check_Λ_lp_feas = setup_check_Λ_lp_feas(nv, mΛ, Ghs!, ∇ᵥF!, ∇ᵥGhs_rows, ∇ᵥGhs_cols, ∇ᵥGhs_vals!, ∇ᵥGhs_shape, Ghs_inds, z_inds)

    # F_path ⟂ θ_l ≤ θ ≤ θ_u
    mΛ_mcp = m1 + mh
    Gh_sym = [G_sym; h_sym]
    @assert(mΛ_mcp == length(Gh_sym))

    nθ = nv + 2 * mΛ_mcp

    θ_inds = Dict([ # defined for convenience
        ("v" => 1:nv),
        ("Λ" => nv+1:nv+mΛ_mcp),
        ("r" => nv+mΛ_mcp+1:nv+2*mΛ_mcp),
        ("ΛG" => nv+1:nv+m1), # part of Λ that corresponds to G
        ("Λh" => nv+m1+1:nv+m1+mh), # part of Λ that corresponds to h (follower's KKT)
        ("z" => n1+1:n1+mh), # part of v that corresponds to z
        ("rG" => nv+mΛ_mcp+1:nv+mΛ_mcp+m1), # part of r that corresponds to G
        ("rh" => nv+mΛ_mcp+m1+1:nv+2*mΛ_mcp), # part of r that corresponds to h (follower's KKT)
    ])
    # θ := [v_sym; Λ_sym; r_sym] 
    θ_l₀ = fill(-Inf, nθ)
    θ_u₀ = fill(Inf, nθ)
    θ_l₀[θ_inds["rh"]] .= z_l₀
    θ_u₀[θ_inds["rh"]] .= z_u₀
    θ_l₀[θ_inds["rG"]] .= 0

    Λ_mcp_sym = Symbolics.@variables(Λ[1:mΛ_mcp])[1] |> Symbolics.scalarize # this Λ ⟂ [G; h]
    r_mcp_sym = Symbolics.@variables(r[1:mΛ_mcp])[1] |> Symbolics.scalarize # slacks for Gh
    θ_sym = [v_sym; Λ_mcp_sym; r_mcp_sym]

    solve_BOPᵢ_KKT_mcp = setup_BOPᵢ_KKT_mcp(nθ, m1, mh, θ_l₀, θ_u₀, F_sym, Gh_sym, v_sym, θ_sym, Λ_mcp_sym, r_mcp_sym)

    #Gh! = Symbolics.build_function(Gh_sym, v_sym; expression=Val{false})[2]
    # TODO: used only for verification, this should be optional to save construction time in the future
    (deriv_funs, sym_derivs) = generate_derivatives(n1, n2, m1, m2, x_sym, F_sym, G_sym, f_sym, g_sym)

    BilevelOptProb(
        n1,
        n2,
        m1,
        m2,
        F,
        G,
        f,
        g,
        solve_follower_nlp,
        solve_follower_KKT_mcp,
        find_bilevel_feas_pt,
        nv,
        mΛ,
        mh,
        v_inds,
        z_inds,
        z_l₀,
        z_u₀,
        Ghs_inds,
        Ghs_l₀,
        Ghs_u₀,
        Ghs!,
        solve_BOPᵢ_nlp,
        info_BOPᵢ,
        eval_BOPᵢ!,
        check_Λ_lp_feas,
        nθ,
        θ_inds,
        θ_l₀,
        θ_u₀,
        solve_BOPᵢ_KKT_mcp,
        deriv_funs,
        sym_derivs,
        ∇ₓ₂ₓ₂L_rows, 
        ∇ₓ₂ₓ₂L_cols, 
        ∇ₓ₂ₓ₂L, 
        ∇ₓ₂ₓ₂L_vals!
    )
end

function setup_follower_nlp(n1, n2, m2, g_l, g_u, x_sym, λ_sym, x2_sym, f_sym, g_sym)
    x2_l = fill(-Inf, n2)
    x2_u = fill(Inf, n2)
    f = Symbolics.build_function(f_sym, x_sym; expression=Val{false})
    g! = Symbolics.build_function(g_sym, x_sym; expression=Val{false})[2]

    ∇ₓ₂f_sym = Symbolics.gradient(f_sym, x2_sym)
    ∇ₓ₂f! = Symbolics.build_function(∇ₓ₂f_sym, x_sym; expression=Val{false})[2]

    ∇ₓ₂g_sym = Symbolics.sparsejacobian(g_sym, x2_sym)
    (∇ₓ₂g_rows, ∇ₓ₂g_cols, ∇ₓ₂g_vals) = SparseArrays.findnz(∇ₓ₂g_sym)
    ∇ₓ₂g_vals! = Symbolics.build_function(∇ₓ₂g_vals, x_sym; expression=Val{false})[2]

    obj_factor = Symbolics.@variables(of)[1]
    if isempty(λ_sym)
        L = obj_factor * f_sym
    else
        L = obj_factor * f_sym + g_sym' * λ_sym # WARN: IPOPT convention: ∇²ₓ₂f(x) + λᵀ ∇²ₓ₂ g(x)
    end

    ∇ₓ₂L = Symbolics.gradient(L, x2_sym)
    ∇ₓ₂ₓ₂L = Symbolics.sparsejacobian(∇ₓ₂L, x2_sym)
    (∇ₓ₂ₓ₂L_rows, ∇ₓ₂ₓ₂L_cols, ∇ₓ₂ₓ₂L_vals_sym) = SparseArrays.findnz(∇ₓ₂ₓ₂L)
    ∇ₓ₂ₓ₂L_vals! = Symbolics.build_function(∇ₓ₂ₓ₂L_vals_sym, x_sym, obj_factor, λ_sym; expression=Val{false})[2]

    x = zeros(n1 + n2)

    function solve_follower_nlp(x1; x2_init=zeros(n2), tol=1e-6, max_iter=1000, print_level=0, is_using_HSL=false)
        x[1:n1] .= x1
        x2_inds = n1+1:n1+n2

        function eval_f(x2::Vector{Float64})
            x[x2_inds] .= x2
            f(x)
        end

        function eval_g(g::Vector{Float64}, x2::Vector{Float64})
            x[x2_inds] .= x2
            g!(g, x)
        end

        function eval_∇ₓ₂f(∇ₓ₂f::Vector{Float64}, x2::Vector{Float64})
            x[x2_inds] .= x2
            ∇ₓ₂f!(∇ₓ₂f, x)
        end

        function eval_∇ₓ₂g_vals(∇ₓ₂g_vals::Vector{Float64}, x2::Vector{Float64})
            x[x2_inds] .= x2
            ∇ₓ₂g_vals!(∇ₓ₂g_vals, x)
        end

        function eval_∇ₓ₂ₓ₂L_vals(∇ₓ₂ₓ₂L_vals::Vector{Float64}, x2::Vector{Float64}, obj_factor::Float64, λ::Vector{Float64})
            x[x2_inds] .= x2
            ∇ₓ₂ₓ₂L_vals!(∇ₓ₂ₓ₂L_vals, x, obj_factor, λ)
        end

        solve = setup_nlp_solve_IPOPT(n2, m2, x2_l, x2_u, g_l, g_u, eval_f, eval_g, eval_∇ₓ₂f, ∇ₓ₂g_rows, ∇ₓ₂g_cols, eval_∇ₓ₂g_vals, ∇ₓ₂ₓ₂L_rows, ∇ₓ₂ₓ₂L_cols, eval_∇ₓ₂ₓ₂L_vals)
        solve(; x_init=x2_init, tol, max_iter, print_level, is_using_HSL)
    end

    (; solve_follower_nlp, ∇ₓ₂ₓ₂L_rows, ∇ₓ₂ₓ₂L_cols, ∇ₓ₂ₓ₂L, ∇ₓ₂ₓ₂L_vals!)
end


function setup_follower_KKT_mcp(n1, n2, m2, x1_sym, x2_sym, λ_sym, s_sym, f_sym, g_sym, z_inds)
    z_sym = [x2_sym; λ_sym; s_sym]
    n_z = length(z_sym)
    z_l = fill(-Inf, n_z)
    z_u = fill(Inf, n_z)
    #θ_l[θ_inds["λ"]] .= 0.0
    z_l[z_inds["s"]] .= 0.0

    if isempty(λ_sym)
        L = f_sym
    else
        L = f_sym - g_sym' * λ_sym
    end

    g_s_sym = g_sym .- s_sym

    ∇ₓ₂L_sym = Symbolics.gradient(L, x2_sym)
    PF_sym = [∇ₓ₂L_sym; g_s_sym; λ_sym] # P(ATH)F function for the PATH solver, not to be confused with leader's cost

    x1_z_sym = [x1_sym; x2_sym; λ_sym; s_sym]
    PF! = Symbolics.build_function(PF_sym, x1_z_sym; expression=Val(false))[2]
    PJ_sym = Symbolics.sparsejacobian(PF_sym, z_sym) # Jacobian of F for the PATH solver

    (PJ_rows, PJ_cols, PJ_vals_sym) = SparseArrays.findnz(PJ_sym)
    PJ_vals! = Symbolics.build_function(PJ_vals_sym, x1_z_sym; expression=Val{false})[2]

    function solve_follower_KKT_mcp(x1; z_init=zeros(n_z), tol=1e-6, max_iter=1000, is_silent=true)
        x1_z = zeros(length(x1_z_sym))
        x1_z[1:n1] .= x1
        z_inds = n1+1:n1+n2+2*m2

        function eval_PF!(F, z::Vector{Float64})
            x1_z[z_inds] .= z
            PF!(F, x1_z)
        end

        function eval_PJ_vals!(J_vals, z::Vector{Float64})
            x1_z[z_inds] .= z
            PJ_vals!(J_vals, x1_z)
        end

        solve = setup_mcp_solve_PATH(n_z, z_l, z_u, eval_PF!, PJ_rows, PJ_cols, eval_PJ_vals!)
        solve(; x_init=z_init, tol, max_iter, is_silent)
    end
end

function setup_find_bile_feas_pt(n1, n2, m1, m2, x_sym, F_sym, G_sym, g_sym)
    x_l = fill(-Inf, n1 + n2)
    x_u = fill(Inf, n1 + n2)
    Gg_l = fill(0.0, m1 + m2)
    Gg_u = fill(Inf, m1 + m2)

    F = Symbolics.build_function(F_sym, x_sym; expression=Val{false})

    ∇ₓF_sym = Symbolics.gradient(F_sym, x_sym)
    ∇ₓF = Symbolics.build_function(∇ₓF_sym, x_sym; expression=Val{false})[2]
    #function f_zero(x::Vector{Float64})
    #    F(x)
    #end

    Gg_sym = [G_sym; g_sym]
    Gg! = Symbolics.build_function(Gg_sym, x_sym; expression=Val{false})[2]

    #function ∇ₓf_zero(∇ₓf::Vector{Float64}, x::Vector{Float64})
    #    ∇ₓf .= 0.0
    #end

    ∇ₓGg_sym = Symbolics.sparsejacobian(Gg_sym, x_sym)
    (∇ₓGg_rows, ∇ₓGg_cols, ∇ₓGg_vals) = SparseArrays.findnz(∇ₓGg_sym)
    ∇ₓGg_vals! = Symbolics.build_function(∇ₓGg_vals, x_sym; expression=Val{false})[2]

    Λ_sym = Symbolics.@variables(lamb[1:m1+m2])[1] |> Symbolics.scalarize
    if isempty(Λ_sym)
        L_feas = 0.0
    else
        L_feas = Gg_sym' * Λ_sym # WARN: IPOPT convention: λᵀ ∇²ₓGg(x) because zero cost
    end

    ∇ₓL = Symbolics.gradient(L_feas, x_sym)
    ∇ₓₓL = Symbolics.sparsejacobian(∇ₓL, x_sym)
    (∇ₓₓL_rows, ∇ₓₓL_cols, ∇ₓₓL_vals_sym) = SparseArrays.findnz(∇ₓₓL)
    ∇ₓₓL_feas_vals! = Symbolics.build_function(∇ₓₓL_vals_sym, x_sym, Λ_sym; expression=Val{false})[2]

    setup_nlp_solve_IPOPT(n1 + n2, m1 + m2, x_l, x_u, Gg_l, Gg_u, F, Gg!, ∇ₓF, ∇ₓGg_rows, ∇ₓGg_cols, ∇ₓGg_vals!, ∇ₓₓL_rows, ∇ₓₓL_cols, ∇ₓₓL_feas_vals!)
end

function setup_BOPᵢ_nlp(nv, mΛ, Ghs_l₀, Ghs_u₀, v_sym, F_sym, Ghs_sym)
    v_l₀ = fill(-Inf, nv)
    v_u₀ = fill(Inf, nv)
    F = Symbolics.build_function(F_sym, v_sym; expression=Val{false})

    Ghs! = Symbolics.build_function(Ghs_sym, v_sym; expression=Val{false})[2]

    ∇ᵥF_sym = Symbolics.gradient(F_sym, v_sym)
    ∇ᵥF! = Symbolics.build_function(∇ᵥF_sym, v_sym; expression=Val{false})[2]

    ∇ᵥGhs_sym = Symbolics.sparsejacobian(Ghs_sym, v_sym)
    (∇ᵥGhs_rows, ∇ᵥGhs_cols, ∇ᵥGhs_vals_sym) = SparseArrays.findnz(∇ᵥGhs_sym)
    ∇ᵥGhs_vals! = Symbolics.build_function(∇ᵥGhs_vals_sym, v_sym; expression=Val{false})[2]

    Λ_sym = Symbolics.@variables(lamb[1:mΛ])[1] |> Symbolics.scalarize
    obj_factor = Symbolics.@variables(of)[1]
    L_sym = obj_factor * F_sym + Ghs_sym' * Λ_sym # WARN: IPOPT convention: ∇ᵥᵥF(v) + Λᵀ ∇ᵥᵥ[G(v); h(v)]
    ∇ᵥL_sym = Symbolics.gradient(L_sym, v_sym)
    ∇ᵥᵥL_sym = Symbolics.sparsejacobian(∇ᵥL_sym, v_sym)
    (∇ᵥᵥL_rows, ∇ᵥᵥL_cols, ∇ᵥᵥL_vals_sym) = SparseArrays.findnz(∇ᵥᵥL_sym)
    ∇ᵥᵥL_vals! = Symbolics.build_function(∇ᵥᵥL_vals_sym, v_sym, obj_factor, Λ_sym; expression=Val{false})[2]

    solve = setup_nlp_solve_IPOPT(nv, mΛ, v_l₀, v_u₀, Ghs_l₀, Ghs_u₀, F, Ghs!, ∇ᵥF!, ∇ᵥGhs_rows, ∇ᵥGhs_cols, ∇ᵥGhs_vals!, ∇ᵥᵥL_rows, ∇ᵥᵥL_cols, ∇ᵥᵥL_vals!)

    # used to verify solutions
    function info()
        (; ∇ᵥGhs_rows, ∇ᵥGhs_cols, ∇ᵥGhs_shape=size(∇ᵥGhs_sym), ∇ᵥᵥL_rows, ∇ᵥᵥL_cols, ∇ᵥᵥL_shape=size(∇ᵥᵥL_sym))
    end

    function eval!(Ghs, ∇ᵥF, ∇ᵥGh_vals, ∇ᵥᵥL_vals, v, Λ)
        Ghs!(Ghs, v)
        ∇ᵥF!(∇ᵥF, v)
        ∇ᵥGhs_vals!(∇ᵥGh_vals, v)
        ∇ᵥᵥL_vals!(∇ᵥᵥL_vals, v, 1.0, Λ)
    end

    (; solve, info, eval!, Ghs!, ∇ᵥF!, ∇ᵥGhs_rows, ∇ᵥGhs_cols, ∇ᵥGhs_vals!, ∇ᵥGhs_shape=size(∇ᵥGhs_sym))
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
        Λ_l = fill(0., mΛ)
        Λ_u = fill(Inf, mΛ)
        # Λₕ=z is complement to h, s is complement to λ
        Λ_l[Ghs_inds["h"]] .= z_l 
        Λ_u[Ghs_inds["h"]] .= z_u
        Λ_l[Ghs_inds["s"]] .= z_l[z_inds["s"]]
        Λ_u[Ghs_inds["s"]] .= z_u[z_inds["s"]]
# Main.@infiltrate
        # constraint matrix is column-wise, this part checks :
        # ∇ᵥF - Λ' * ∇ᵥGhs = 0 (stationarity)
        # Λ' * Ghs = 0 (complementarity)
        A_l = [∇ᵥF; 0;]
        A_u = [∇ᵥF; 0;]
        A = [∇ᵥGhs'; Ghs';] # violates constraint qualifications like this
        #Main.@infiltrate

        check_feas(Λ_l, Λ_u, A_l, A_u, A)
    end
end


function setup_BOPᵢ_KKT_mcp(nθ, m1, mh, θ_l₀, θ_u₀, F_sym, Gh_sym, v_sym, θ_sym, Λ_sym, r_sym)
    L = F_sym - Gh_sym' * Λ_sym
    ∇ᵥL_sym = Symbolics.gradient(L, v_sym)
    Gh_r_sym = Gh_sym .- r_sym
    # F ⟂ θ = [v_sym; Λ_sym; r_sym] 
    PF_sym = [∇ᵥL_sym; Gh_r_sym; Λ_sym]

    PF! = Symbolics.build_function(PF_sym, θ_sym; expression=Val(false))[2]
    PJ = Symbolics.sparsejacobian(PF_sym, θ_sym)

    (PJ_rows, PJ_cols, PJ_vals) = SparseArrays.findnz(PJ)
    PJ_vals! = Symbolics.build_function(PJ_vals, θ_sym; expression=Val{false})[2]

    setup_mcp_solve_PATH(nθ, θ_l₀, θ_u₀, PF!, PJ_rows, PJ_cols, PJ_vals!)
end

"""
Generate derivative functions (currently unused 2025-06-20)
"""
function generate_derivatives(n1, n2, m1, m2, x_sym, F_sym, G_sym, f_sym, g_sym)
    nx = n1 + n2 # length(x)
    # First-order derivatives --- ∇ₓ
    ∇ₓF_sym = Symbolics.gradient(F_sym, x_sym)
    ∇ₓf_sym = Symbolics.gradient(f_sym, x_sym)
    ∇ₓG_sym = Symbolics.sparsejacobian(G_sym, x_sym)
    ∇ₓg_sym = Symbolics.sparsejacobian(g_sym, x_sym)
    @assert(nx == length(∇ₓF_sym))
    @assert(nx == length(∇ₓf_sym))
    @assert((m1, nx) == size(∇ₓG_sym))
    @assert((m2, nx) == size(∇ₓg_sym))
    # ∇ₓF!(val, x)
    ∇ₓF! = Symbolics.build_function(∇ₓF_sym, x_sym; expression=Val{false})[2]
    # ∇ₓf!(val, x)
    ∇ₓf! = Symbolics.build_function(∇ₓf_sym, x_sym; expression=Val{false})[2]
    # ∇ₓG_rows, ∇ₓG_cols, ∇ₓG!(val, x)
    (∇ₓG_rows, ∇ₓG_cols, ∇ₓG_vals_sym) = SparseArrays.findnz(∇ₓG_sym)
    ∇ₓG_vals! = Symbolics.build_function(∇ₓG_vals_sym, x_sym; expression=Val{false})[2]
    # ∇ₓg_rows, ∇ₓg_cols, ∇ₓg!(val, x)
    (∇ₓg_rows, ∇ₓg_cols, ∇ₓg_vals_sym) = SparseArrays.findnz(∇ₓg_sym)
    ∇ₓg_vals! = Symbolics.build_function(∇ₓg_vals_sym, x_sym; expression=Val{false})[2]
    # Second-order derivatives --- ∇ₓ²
    ∇ₓₓF_sym = Symbolics.sparsejacobian(∇ₓF_sym, x_sym)
    ∇ₓₓf_sym = Symbolics.sparsejacobian(∇ₓf_sym, x_sym)
    ∇ₓₓG_sym = Symbolics.sparsejacobian(vec(Matrix(∇ₓG_sym)), x_sym) # passing a full vectorized matrix
    ∇ₓₓg_sym = Symbolics.sparsejacobian(vec(Matrix(∇ₓg_sym)), x_sym)
    @assert((nx, nx) == size(∇ₓₓF_sym))
    @assert((nx, nx) == size(∇ₓₓf_sym))
    @assert((m1 * nx, nx) == size(∇ₓₓG_sym))
    @assert((m2 * nx, nx) == size(∇ₓₓg_sym))
    # ∇ₓ²F_rows, ∇ₓ²F_cols, ∇ₓ²F!(val, x)
    (∇ₓₓF_rows, ∇ₓₓF_cols, ∇ₓₓF_vals_sym) = SparseArrays.findnz(∇ₓₓF_sym)
    ∇ₓₓF_vals! = Symbolics.build_function(∇ₓₓF_vals_sym, x_sym; expression=Val{false})[2]
    # ∇ₓ²f_rows, ∇ₓ²f_cols, ∇ₓ²f!(val, x)
    (∇ₓₓf_rows, ∇ₓₓf_cols, ∇ₓₓf_vals_sym) = SparseArrays.findnz(∇ₓₓf_sym)
    ∇ₓₓf_vals! = Symbolics.build_function(∇ₓₓf_vals_sym, x_sym; expression=Val{false})[2]
    # ∇ₓ²G_rows, ∇ₓ²G_cols, ∇ₓ²G!(val, x)
    (∇ₓₓG_rows, ∇ₓₓG_cols, ∇ₓₓG_vals_sym) = SparseArrays.findnz(∇ₓₓG_sym)
    ∇ₓₓG_vals! = Symbolics.build_function(∇ₓₓG_vals_sym, x_sym; expression=Val{false})[2]
    # ∇ₓ²g_rows, ∇ₓ²g_cols, ∇ₓ²g!(val, x)
    (∇ₓₓg_rows, ∇ₓₓg_cols, ∇ₓₓg_vals_sym) = SparseArrays.findnz(∇ₓₓg_sym)
    ∇ₓₓg_vals! = Symbolics.build_function(∇ₓₓg_vals_sym, x_sym; expression=Val{false})[2]

    funs = (;
        ∇ₓF!,
        ∇ₓf!,
        ∇ₓG_rows,
        ∇ₓG_cols,
        ∇ₓG_vals!,
        ∇ₓg_rows,
        ∇ₓg_cols,
        ∇ₓg_vals!,
        ∇ₓₓF_rows,
        ∇ₓₓF_cols,
        ∇ₓₓF_vals!,
        ∇ₓₓf_rows,
        ∇ₓₓf_cols,
        ∇ₓₓf_vals!,
        ∇ₓₓG_rows,
        ∇ₓₓG_cols,
        ∇ₓₓG_vals!,
        ∇ₓₓg_rows,
        ∇ₓₓg_cols,
        ∇ₓₓg_vals!
    )
    syms = (;
        F=F_sym,
        f=f_sym,
        G=G_sym,
        g=g_sym,
        ∇ₓF=∇ₓF_sym,
        ∇ₓf=∇ₓf_sym,
        ∇ₓG=∇ₓG_sym,
        ∇ₓg=∇ₓg_sym,
        ∇ₓₓF=∇ₓₓF_sym,
        ∇ₓₓf=∇ₓₓf_sym,
        ∇ₓₓG=∇ₓₓG_sym,
        ∇ₓₓg=∇ₓₓg_sym,
    )
    (funs, syms)
end
