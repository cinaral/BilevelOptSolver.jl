struct BilevelOptProb
    n₁::Int # length(x₁)
    n₂::Int # length(x₂)
    m₁::Int # length(G(x))
    m₂::Int # length(g(x))
    F::Function # F(x)
    G::Function # G(x)
    f::Function # f(x)
    g::Function # g(x)
    solve_follower_nlp::Function
    solve_follower_KKT_mcp::Function
    find_bilevel_feas_pt::Function
    nᵥ::Int  # length(v) = n₁ + n₂ + m₂ + m₂ 
    mₕ::Int  # length(h) = n₂ + m₂ + m₂ 
    v_inds::Dict{String,UnitRange{Int64}} # indexing v := [x₁; z], z := [x₂; λ; s]
    v_l₀::Vector{Float64} # default v bounds 
    v_u₀::Vector{Float64}
    Gh_inds::Dict{String,UnitRange{Int64}}
    Gh_l₀::Vector{Float64} # default [G; h] bounds 
    Gh_u₀::Vector{Float64}
    Gh!::Function # Gh(out, v), out := [G(v); h(v)]
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
function construct_bop(n₁, n₂, F, G, f, g; verbosity=0)
    nₓ = n₁ + n₂ # length(x)
    x_dummy = zeros(nₓ)
    m₁ = length(G(x_dummy))
    m₂ = length(g(x_dummy))
    nᵥ = n₁ + n₂ + m₂ + m₂ # length(v)
    mₕ = n₂ + m₂ + m₂ # length(h)
    x₁_sym = Symbolics.@variables(x₁[1:n₁])[1] |> Symbolics.scalarize
    x₂_sym = Symbolics.@variables(x₂[1:n₂])[1] |> Symbolics.scalarize
    λ_sym = Symbolics.@variables(λ[1:m₂])[1] |> Symbolics.scalarize
    s_sym = Symbolics.@variables(s[1:m₂])[1] |> Symbolics.scalarize
    z_sym = [x₂_sym; λ_sym; s_sym] # z = [x₂; λ; s] ∈ R⁽ᵐʰ⁾
    v_sym = [x₁_sym; z_sym] # v = [x₁; x₂; λ; s] ∈ R⁽ⁿᵛ⁾
    @assert(mₕ == length(z_sym))
    @assert(nᵥ == length(v_sym))

    # defined for convenience
    z_inds = Dict([
        ("x₂", 1:n₂),
        ("λ", n₂+1:n₂+m₂),
        ("s", n₂+m₂+1:n₂+m₂+m₂),
    ])
    v_inds = Dict([
        ("x₁", 1:n₁),
        ("z", n₁+1:n₁+mₕ),
        ("x", 1:nₓ),
        ("x₂", n₁+1:nₓ),
        ("λ", nₓ+1:nₓ+m₂),
        ("s", nₓ+m₂+1:nₓ+2*m₂)
    ])

    x_sym = [x₁_sym; x₂_sym]
    F_sym = F(x_sym)
    G_sym = G(x_sym)
    f_sym = f(x_sym)
    g_sym = g(x_sym)

    ∇ₓ₂f = Symbolics.gradient(f_sym, x₂_sym)
    ∇ₓ₂g = Symbolics.sparsejacobian(g_sym, x₂_sym)

    if isempty(λ_sym)
        h_sym = [∇ₓ₂f; g_sym - s_sym; λ_sym]
    else
        h_sym = [∇ₓ₂f - ∇ₓ₂g' * λ_sym; g_sym - s_sym; λ_sym]
    end
    @assert(mₕ == length(h_sym))

    # hᵢ ⟂ z_lᵢ ≤ z ≤ z_uᵢ
    # x₁, x₂, λ are free, but 0 ≤ s ≤ s_ubᵢ
    z_l₀ = fill(-Inf, mₕ) # default z lb
    z_l₀[z_inds["s"]] .= zeros(m₂)
    z_u₀ = fill(Inf, mₕ) # default z ub

    # If there were bounds on leader x₁ it could be added here
    v_l₀ = fill(-Inf, nᵥ) # default l lb
    v_l₀[v_inds["z"]] .= z_l₀
    v_u₀ = fill(Inf, nᵥ) # default u lb
    v_u₀[v_inds["z"]] .= z_u₀

    x₂_l = v_l₀[v_inds["x₂"]]
    x₂_u = v_u₀[v_inds["x₂"]]
    g_l = fill(0.0, m₂)
    g_u = fill(Inf, m₂)
    solve_follower_nlp = setup_follower_nlp(n₁, n₂, m₂, x₂_l, x₂_u, g_l, g_u, x_sym, λ_sym, x₂_sym, f_sym, g_sym)
    solve_follower_KKT_mcp = setup_follower_KKT_mcp(n₁, n₂, m₂, x_sym, x₂_sym, λ_sym, s_sym, f_sym, g_sym)

    find_bilevel_feas_pt = setup_find_bile_feas_pt(n₁, n₂, m₁, m₂, x_sym, G_sym, g_sym)

    m = m₁ + mₕ  # length([G; h])
    Gh_inds = Dict([ # defined for convenience
        ("G", 1:m₁),
        ("h", m₁+1:m),
        ("x₂", m₁+1:m₁+n₂),
        ("λ", m₁+n₂+1:m₁+n₂+m₂),
        ("s", m₁+n₂+m₂+1:m)
    ])
    Gh_l₀ = zeros(m) # default Gh lb (shouldn't change)
    Gh_u₀ = fill(Inf, m) # default Gh ub
    Gh_u₀[Gh_inds["x₂"]] .= 0 # x₂ is free so its complement is zero
    Gh_u₀[Gh_inds["λ"]] .= 0 # λ is free so its complement is zero

    Gh_sym = [G_sym; h_sym] # BOPᵢ constraints
    @assert(m == length(Gh_sym))
    Gh_inds = Dict([("G", 1:m₁), ("h", m₁+1:m₁+mₕ)]) # defined for convenience

    solve_BOPᵢ_nlp, info_BOPᵢ, eval_BOPᵢ!, Gh!, ∇ᵥF!, ∇ᵥGh_rows, ∇ᵥGh_cols, ∇ᵥGh_vals!, ∇ᵥGh_shape = setup_BOPᵢ_nlp(nᵥ, m, v_l₀, v_u₀, Gh_l₀, Gh_u₀, v_sym, F_sym, Gh_sym)

    check_Λ_lp_feas = setup_check_Λ_lp_feas(nᵥ, m₁, mₕ, Gh!, ∇ᵥF!, ∇ᵥGh_rows, ∇ᵥGh_cols, ∇ᵥGh_vals!, ∇ᵥGh_shape, Gh_inds, v_inds)

    # F_path ⟂ θ_l ≤ θ ≤ θ_u
    Λ_sym = Symbolics.@variables(Λ[1:m])[1] |> Symbolics.scalarize
    r_sym = Symbolics.@variables(r[1:m])[1] |> Symbolics.scalarize # slacks for Gh
    θ_sym = [v_sym; Λ_sym; r_sym]
    n_θ = length(θ_sym)
    θ_inds = Dict([ # defined for convenience
        ("v" => 1:nᵥ),
        ("Λ" => nᵥ+1:nᵥ+m),
        ("Λ₁" => nᵥ+1:nᵥ+m₁),
        ("Λₕ" => nᵥ+m₁+1:nᵥ+m),
        ("r" => nᵥ+m+1:nᵥ+2*m),
        ("z" => n₁+1:n₁+mₕ),
        ("r₁" => nᵥ+m+1:nᵥ+m+m₁),
        ("rₕ" => nᵥ+m+m₁+1:nᵥ+2*m),
    ])
    θ_l₀ = fill(-Inf, n_θ)
    θ_l₀[θ_inds["z"]] .= z_l₀
    θ_l₀[θ_inds["Λ₁"]] .= 0
    θ_u₀ = fill(Inf, n_θ)
    θ_u₀[θ_inds["z"]] .= z_u₀
    solve_BOPᵢ_KKT_mcp = setup_BOPᵢ_KKT_mcp(n_θ, θ_l₀, θ_u₀, θ_sym, F_sym, Gh_sym, v_sym, Λ_sym, r_sym)

    # TODO: used only for verification, this should be optional to save construction time in the future
    (deriv_funs, sym_derivs) = generate_derivatives(n₁, n₂, m₁, m₂, x_sym, F_sym, G_sym, f_sym, g_sym)

    BilevelOptProb(
        n₁,
        n₂,
        m₁,
        m₂,
        F,
        G,
        f,
        g,
        solve_follower_nlp,
        solve_follower_KKT_mcp,
        find_bilevel_feas_pt,
        nᵥ,
        mₕ,
        v_inds,
        v_l₀,
        v_u₀,
        Gh_inds,
        Gh_l₀,
        Gh_u₀,
        Gh!,
        solve_BOPᵢ_nlp,
        info_BOPᵢ,
        eval_BOPᵢ!,
        check_Λ_lp_feas,
        n_θ,
        θ_inds,
        θ_l₀,
        θ_u₀,
        solve_BOPᵢ_KKT_mcp,
        deriv_funs,
        sym_derivs
    )
end

function setup_follower_nlp(n₁, n₂, m₂, x₂_l, x₂_u, g_l, g_u, x_sym, λ_sym, x₂_sym, f_sym, g_sym)
    f = Symbolics.build_function(f_sym, x_sym; expression=Val{false})
    g! = Symbolics.build_function(g_sym, x_sym; expression=Val{false})[2]

    ∇ₓ₂f_sym = Symbolics.gradient(f_sym, x₂_sym)
    ∇ₓ₂f! = Symbolics.build_function(∇ₓ₂f_sym, x_sym; expression=Val{false})[2]

    ∇ₓ₂g_sym = Symbolics.sparsejacobian(g_sym, x₂_sym)
    (∇ₓ₂g_rows, ∇ₓ₂g_cols, ∇ₓ₂g_vals) = SparseArrays.findnz(∇ₓ₂g_sym)
    ∇ₓ₂g_vals! = Symbolics.build_function(∇ₓ₂g_vals, x_sym; expression=Val{false})[2]

    obj_factor = Symbolics.@variables(σf)[1]
    if isempty(λ_sym)
        L = obj_factor * f_sym
    else
        L = obj_factor * f_sym + g_sym' * λ_sym # WARN: IPOPT convention: ∇²ₓ₂f(x) + λᵀ ∇²ₓ₂ g(x)
    end

    ∇ₓ₂L = Symbolics.gradient(L, x₂_sym)
    ∇²ₓ₂L = Symbolics.sparsejacobian(∇ₓ₂L, x₂_sym)
    (∇²ₓ₂L_rows, ∇²ₓ₂L_cols, ∇²ₓ₂L_vals_sym) = SparseArrays.findnz(∇²ₓ₂L)
    ∇²ₓ₂L_vals! = Symbolics.build_function(∇²ₓ₂L_vals_sym, x_sym, obj_factor, λ_sym; expression=Val{false})[2]

    x = zeros(n₁ + n₂)

    function solve_follower_nlp(x₁; x₂_init=zeros(n₂), tol=1e-6, max_iter=1000, print_level=0, is_using_HSL=false)
        x[1:n₁] .= x₁

        function eval_f(x₂::Vector{Float64})
            x[n₁+1:end] .= x₂
            f(x)
        end

        function eval_g(g::Vector{Float64}, x₂::Vector{Float64})
            x[n₁+1:end] .= x₂
            g!(g, x)
        end

        function eval_∇ₓ₂f(∇ₓ₂f::Vector{Float64}, x₂::Vector{Float64})
            x[n₁+1:end] .= x₂
            ∇ₓ₂f!(∇ₓ₂f, x)
        end

        function eval_∇ₓ₂g_vals(∇ₓ₂g_vals::Vector{Float64}, x₂::Vector{Float64})
            x[n₁+1:end] .= x₂
            ∇ₓ₂g_vals!(∇ₓ₂g_vals, x)
        end

        function eval_∇²ₓ₂L_vals(∇²ₓ₂L_vals::Vector{Float64}, x₂::Vector{Float64}, obj_factor::Float64, λ::Vector{Float64})
            x[n₁+1:end] .= x₂
            ∇²ₓ₂L_vals!(∇²ₓ₂L_vals, x, obj_factor, λ)
        end

        solve = setup_nlp_solve_IPOPT(n₂, m₂, x₂_l, x₂_u, g_l, g_u, eval_f, eval_g, eval_∇ₓ₂f, ∇ₓ₂g_rows, ∇ₓ₂g_cols, eval_∇ₓ₂g_vals, ∇²ₓ₂L_rows, ∇²ₓ₂L_cols, eval_∇²ₓ₂L_vals)
        solve(; x_init=x₂_init, tol, max_iter, print_level, is_using_HSL)
    end
end


function setup_follower_KKT_mcp(n₁, n₂, m₂, x_sym, x₂_sym, λ_sym, s_sym, f_sym, g_sym)
    θ_sym = [x₂_sym; λ_sym; s_sym]
    n_θ = length(θ_sym)
    θ_inds = Dict([ # defined for convenience
        ("x₂" => 1:n₂),
        ("λ" => n₂+1:n₂+m₂),
        ("s" => n₂+m₂+1:n₂+2*m₂)
    ])
    θ_l = fill(-Inf, n_θ)
    θ_l[θ_inds["λ"]] .= 0.0
    θ_l[θ_inds["s"]] .= 0.0
    θ_u = fill(Inf, n_θ)

    g_w_slack = g_sym .- s_sym

    if isempty(λ_sym)
        L = f_sym
    else
        L = f_sym + g_w_slack' * λ_sym
    end

    ∇ₓ₂L_sym = Symbolics.gradient(L, x₂_sym)
    F_sym = [∇ₓ₂L_sym; g_w_slack; λ_sym]

    θ_w_x₁ = [x_sym; λ_sym; s_sym] # 
    F! = Symbolics.build_function(F_sym, θ_w_x₁; expression=Val(false))[2]
    J = Symbolics.sparsejacobian(F_sym, θ_sym)

    (J_rows, J_cols, J_vals) = SparseArrays.findnz(J)
    J_vals! = Symbolics.build_function(J_vals, θ_w_x₁; expression=Val{false})[2]

    function solve_follower_KKT_mcp(x₁; θ_init=zeros(n_θ), tol=1e-6, max_iter=1000, is_silent=true)
        θ_w_x₁ = zeros(length(θ_w_x₁))
        θ_w_x₁[1:n₁] .= x₁

        function eval_F!(F, θ::Vector{Float64})
            θ_w_x₁[n₁+1:end] .= θ
            F!(F, θ_w_x₁)
        end

        function eval_J_vals!(J_vals, θ::Vector{Float64})
            θ_w_x₁[n₁+1:end] .= θ
            J_vals!(J_vals, θ_w_x₁)
        end

        solve = setup_mcp_solve_PATH(n_θ, θ_l, θ_u, F!, J_rows, J_cols, J_vals!)
        solve(; x_init=θ_init, tol, max_iter, is_silent)
    end
end

function setup_find_bile_feas_pt(n₁, n₂, m₁, m₂, x_sym, G_sym, g_sym)
    x_l = fill(-Inf, n₁ + n₂)
    x_u = fill(Inf, n₁ + n₂)
    Gg_l = fill(0.0, m₁ + m₂)
    Gg_u = fill(Inf, m₁ + m₂)

    function f_zero(x::Vector{Float64})
        0.0
    end

    Gg_sym = [G_sym; g_sym]
    Gg! = Symbolics.build_function(Gg_sym, x_sym; expression=Val{false})[2]

    function ∇ₓf_zero(∇ₓf::Vector{Float64}, x::Vector{Float64})
        ∇ₓf .= 0.0
    end

    ∇ₓGg_sym = Symbolics.sparsejacobian(Gg_sym, x_sym)
    (∇ₓGg_rows, ∇ₓGg_cols, ∇ₓGg_vals) = SparseArrays.findnz(∇ₓGg_sym)
    ∇ₓGg_vals! = Symbolics.build_function(∇ₓGg_vals, x_sym; expression=Val{false})[2]

    Λ_sym = Symbolics.@variables(Λf[1:m₁+m₂])[1] |> Symbolics.scalarize
    if isempty(Λ_sym)
        L_feas = 0.0
    else
        L_feas = Λ_sym' * Gg_sym # WARN: IPOPT convention: λᵀ ∇²ₓGg(x) because zero cost
    end

    ∇ₓL = Symbolics.gradient(L_feas, x_sym)
    ∇ₓ²L = Symbolics.sparsejacobian(∇ₓL, x_sym)
    (∇ₓ²L_rows, ∇ₓ²L_cols, ∇ₓ²L_vals_sym) = SparseArrays.findnz(∇ₓ²L)
    ∇²ₓL_feas_vals! = Symbolics.build_function(∇ₓ²L_vals_sym, x_sym, Λ_sym; expression=Val{false})[2]

    setup_nlp_solve_IPOPT(n₁ + n₂, m₁ + m₂, x_l, x_u, Gg_l, Gg_u, f_zero, Gg!, ∇ₓf_zero, ∇ₓGg_rows, ∇ₓGg_cols, ∇ₓGg_vals!, ∇ₓ²L_rows, ∇ₓ²L_cols, ∇²ₓL_feas_vals!)
end

function setup_BOPᵢ_nlp(nᵥ, m, v_l₀, v_u₀, Gh_l₀, Gh_u₀, v_sym, F_sym, Gh_sym)
    F = Symbolics.build_function(F_sym, v_sym; expression=Val{false})
    Gh! = Symbolics.build_function(Gh_sym, v_sym; expression=Val{false})[2]

    ∇ᵥF_sym = Symbolics.gradient(F_sym, v_sym)
    ∇ᵥF! = Symbolics.build_function(∇ᵥF_sym, v_sym; expression=Val{false})[2]

    ∇ᵥGh_sym = Symbolics.sparsejacobian(Gh_sym, v_sym)
    (∇ᵥGh_rows, ∇ᵥGh_cols, ∇ᵥGh_vals_sym) = SparseArrays.findnz(∇ᵥGh_sym)
    ∇ᵥGh_vals! = Symbolics.build_function(∇ᵥGh_vals_sym, v_sym; expression=Val{false})[2]

    Λ = Symbolics.@variables(Λ[1:m])[1] |> Symbolics.scalarize
    obj_factor = Symbolics.@variables(σf)[1]
    L_sym = obj_factor * F_sym + Gh_sym' * Λ # WARN: IPOPT convention: ∇²ᵥF(v) + Λᵀ ∇²ᵥ[G(v); h(v)]
    ∇ᵥL_sym = Symbolics.gradient(L_sym, v_sym)
    ∇ᵥ²L_sym = Symbolics.sparsejacobian(∇ᵥL_sym, v_sym)
    (∇ᵥ²L_rows, ∇ᵥ²L_cols, ∇ᵥ²L_vals_sym) = SparseArrays.findnz(∇ᵥ²L_sym)
    ∇ᵥ²L_vals! = Symbolics.build_function(∇ᵥ²L_vals_sym, v_sym, obj_factor, Λ; expression=Val{false})[2]

    solve = setup_nlp_solve_IPOPT(nᵥ, m, v_l₀, v_u₀, Gh_l₀, Gh_u₀, F, Gh!, ∇ᵥF!, ∇ᵥGh_rows, ∇ᵥGh_cols, ∇ᵥGh_vals!, ∇ᵥ²L_rows, ∇ᵥ²L_cols, ∇ᵥ²L_vals!)

    # used to verify solutions
    function info()
        (; ∇ᵥGh_rows, ∇ᵥGh_cols, ∇ᵥGh_shape=size(∇ᵥGh_sym), ∇ᵥ²L_rows, ∇ᵥ²L_cols, ∇ᵥ²L_shape=size(∇ᵥ²L_sym))
    end

    function eval!(Gh, ∇ᵥF, ∇ᵥGh_vals, ∇ᵥ²L_vals, v, Λ)
        Gh!(Gh, v)
        ∇ᵥF!(∇ᵥF, v)
        ∇ᵥGh_vals!(∇ᵥGh_vals, v)
        ∇ᵥ²L_vals!(∇ᵥ²L_vals, v, 1.0, Λ)
    end

    (; solve, info, eval!, Gh!, ∇ᵥF!, ∇ᵥGh_rows, ∇ᵥGh_cols, ∇ᵥGh_vals!, ∇ᵥGh_shape=size(∇ᵥGh_sym))
end

"""
```
KKT conditions for BOPᵢ is practically an LP feasibility problem when v is given:
    ∃Λ ∈ R⁽ᵐ¹⁺ᵐʰ⁾, ∃Λ_l ∈ R⁽ⁿᵛ⁾, ∃Λ_u ∈ R⁽ⁿᵛ⁾:
        ∇ᵥF(x) - Λᵀ ∇ᵥGhᵢ(v) - Λ_v_l + Λ_v_u = 0                                            
                Ghᵢ(v) ≥ 0 ⟂   Λ ≥ 0              (KKT conditions of BOPᵢ NLP)
              v - v_lᵢ ≥ 0 ⟂ Λ_v_l ≥ 0                
             -v + v_uᵢ ≥ 0 ⟂ Λ_v_u ≥ 0 
```
"""
function setup_check_Λ_lp_feas(nᵥ, m₁, mₕ, Gh!, ∇ᵥF!, ∇ᵥGh_rows, ∇ᵥGh_cols, ∇ᵥGh_vals!, ∇ᵥGh_shape, Gh_inds, v_inds)
    m = m₁ + mₕ
    check_feas = setup_lp_feas_check_HiGHS(m + 2 * nᵥ)

    function check_Λ_lp_feas(v, z_l, z_u, h_l, h_u)
        Gh = zeros(m)
        Gh!(Gh, v)
        ∇ᵥF = zeros(nᵥ)
        ∇ᵥF!(∇ᵥF, v)
        ∇ᵥGh = sparse(∇ᵥGh_rows, ∇ᵥGh_cols, zeros(Cdouble, length(∇ᵥGh_rows)), ∇ᵥGh_shape[1], ∇ᵥGh_shape[2])
        ∇ᵥGh_vals!(∇ᵥGh.nzval, v)
        ∇ᵥGhb = [∇ᵥGh; LinearAlgebra.I(nᵥ); -LinearAlgebra.I(nᵥ)] # appended with v bounds 

        # Λ_all := [Λ; Λ_v_l; Λ_v_u]
        # Λ_all_l ≤ Λ_all ≤ Λ_all_u
        Λ_all_l = [fill(-Inf, m); zeros(2 * nᵥ)]
        Λ_all_l[Gh_inds["h"]] = z_l
        Λ_all_u = [fill(Inf, m); fill(Inf, 2 * nᵥ)]
        Λ_all_u[Gh_inds["h"]] = z_u

        # Define the row lower bounds and upper bounds A_l = A_u = b
        # A_l ≤ A * Λ_all ≤ A_u
        A_l = [∇ᵥF; 0]
        A_l[v_inds["z"]] = h_l
        A_u = [∇ᵥF; 0]
        A_u[v_inds["z"]] = h_u

        # constraint matrix is column-wise:
        G = @view Gh[Gh_inds["G"]]
        A = vcat(∇ᵥGhb', sparse([G' zeros(mₕ + 2 * nᵥ)'])) # Λ₁ᵀ G = 0 (leader complementarity) added

        #check_feas(Λ_all_l, Λ_all_u, A_l, A_u, A)
        (; check_feas, Λ_all_l, Λ_all_u, A_l, A_u, A)
    end
end


function setup_BOPᵢ_KKT_mcp(n_θ, θ_l₀, θ_u₀, θ_sym, F_sym, Gh_sym, v_sym, Λ_sym, r_sym)
    Gh_w_slack = Gh_sym .- r_sym

    L = F_sym + Gh_sym' * Λ_sym
    ∇ᵥL_sym = Symbolics.gradient(L, v_sym)
    F_sym = [∇ᵥL_sym; Gh_w_slack; Λ_sym]

    F! = Symbolics.build_function(F_sym, θ_sym; expression=Val(false))[2]
    J = Symbolics.sparsejacobian(F_sym, θ_sym)

    (J_rows, J_cols, J_vals) = SparseArrays.findnz(J)
    J_vals! = Symbolics.build_function(J_vals, θ_sym; expression=Val{false})[2]

    setup_mcp_solve_PATH(n_θ, θ_l₀, θ_u₀, F!, J_rows, J_cols, J_vals!)
end

"""
Generate derivative functions (currently unused 2025-06-20)
"""
function generate_derivatives(n₁, n₂, m₁, m₂, x_sym, F_sym, G_sym, f_sym, g_sym)
    nₓ = n₁ + n₂ # length(x)
    # First-order derivatives --- ∇ₓ
    ∇ₓF_sym = Symbolics.gradient(F_sym, x_sym)
    ∇ₓf_sym = Symbolics.gradient(f_sym, x_sym)
    ∇ₓG_sym = Symbolics.sparsejacobian(G_sym, x_sym)
    ∇ₓg_sym = Symbolics.sparsejacobian(g_sym, x_sym)
    @assert(nₓ == length(∇ₓF_sym))
    @assert(nₓ == length(∇ₓf_sym))
    @assert((m₁, nₓ) == size(∇ₓG_sym))
    @assert((m₂, nₓ) == size(∇ₓg_sym))
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
    ∇ₓ²F_sym = Symbolics.sparsejacobian(∇ₓF_sym, x_sym)
    ∇ₓ²f_sym = Symbolics.sparsejacobian(∇ₓf_sym, x_sym)
    ∇ₓ²G_sym = Symbolics.sparsejacobian(vec(Matrix(∇ₓG_sym)), x_sym) # passing a full vectorized matrix
    ∇ₓ²g_sym = Symbolics.sparsejacobian(vec(Matrix(∇ₓg_sym)), x_sym)
    @assert((nₓ, nₓ) == size(∇ₓ²F_sym))
    @assert((nₓ, nₓ) == size(∇ₓ²f_sym))
    @assert((m₁ * nₓ, nₓ) == size(∇ₓ²G_sym))
    @assert((m₂ * nₓ, nₓ) == size(∇ₓ²g_sym))
    # ∇ₓ²F_rows, ∇ₓ²F_cols, ∇ₓ²F!(val, x)
    (∇ₓ²F_rows, ∇ₓ²F_cols, ∇ₓ²F_vals_sym) = SparseArrays.findnz(∇ₓ²F_sym)
    ∇ₓ²F_vals! = Symbolics.build_function(∇ₓ²F_vals_sym, x_sym; expression=Val{false})[2]
    # ∇ₓ²f_rows, ∇ₓ²f_cols, ∇ₓ²f!(val, x)
    (∇ₓ²f_rows, ∇ₓ²f_cols, ∇ₓ²f_vals_sym) = SparseArrays.findnz(∇ₓ²f_sym)
    ∇ₓ²f_vals! = Symbolics.build_function(∇ₓ²f_vals_sym, x_sym; expression=Val{false})[2]
    # ∇ₓ²G_rows, ∇ₓ²G_cols, ∇ₓ²G!(val, x)
    (∇ₓ²G_rows, ∇ₓ²G_cols, ∇ₓ²G_vals_sym) = SparseArrays.findnz(∇ₓ²G_sym)
    ∇ₓ²G_vals! = Symbolics.build_function(∇ₓ²G_vals_sym, x_sym; expression=Val{false})[2]
    # ∇ₓ²g_rows, ∇ₓ²g_cols, ∇ₓ²g!(val, x)
    (∇ₓ²g_rows, ∇ₓ²g_cols, ∇ₓ²g_vals_sym) = SparseArrays.findnz(∇ₓ²g_sym)
    ∇ₓ²g_vals! = Symbolics.build_function(∇ₓ²g_vals_sym, x_sym; expression=Val{false})[2]

    funs = (;
        ∇ₓF!,
        ∇ₓf!,
        ∇ₓG_rows,
        ∇ₓG_cols,
        ∇ₓG_vals!,
        ∇ₓg_rows,
        ∇ₓg_cols,
        ∇ₓg_vals!,
        ∇ₓ²F_rows,
        ∇ₓ²F_cols,
        ∇ₓ²F_vals!,
        ∇ₓ²f_rows,
        ∇ₓ²f_cols,
        ∇ₓ²f_vals!,
        ∇ₓ²G_rows,
        ∇ₓ²G_cols,
        ∇ₓ²G_vals!,
        ∇ₓ²g_rows,
        ∇ₓ²g_cols,
        ∇ₓ²g_vals!
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
        ∇ₓ²F=∇ₓ²F_sym,
        ∇ₓ²f=∇ₓ²f_sym,
        ∇ₓ²G=∇ₓ²G_sym,
        ∇ₓ²g=∇ₓ²g_sym,
    )
    (funs, syms)
end
