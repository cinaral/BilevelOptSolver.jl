struct BilevelOptProb
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

"""
Construct BilevelOptProb
```
Assume, 
    F: Rⁿ¹->R, 
    f: Rⁿ²->R, 
    G: X->Rᵐ¹, 
    g: X->Rᵐ², 
    are all twice continuously differentiable.

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
"""
function construct_bop(n1, n2, F, G, f, g; np=0, verbosity=0)
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

    # follower NLP: min f(x) s.t. g(x) ≥ 0
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

    bop = BilevelOptProb(
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

"""
    x := [x1; x2] and optionally [x1; x2; p] to let F, G, f, g functions accept paramaters
    z := [x2; λ]
    v := [x1; z] and optionally [x1; z; p] to let F, G, f, g functions accept paramaters
    Γ := [G; h; -h; z; -z] 
    θ := [v; Λ] where Λ'Γ=0 (so pars   and optionally [v; Λ; p] to let F, G, f, g functions accept paramaters
"""
function define_index_dicts(n1, n2, m1, m2, nx, np, nz, n, m)
    x_inds = Dict{String,UnitRange{Int64}}([
        ("x1", 1:n1),
        ("x2", n1+1:n1+n2),
        ("x", 1:n1+n2),
        ("p", n1+n2+1:n1+n2+np)
    ])

    z_inds = Dict{String,UnitRange{Int64}}([
        ("x2", 1:n2)
        ("λ", n2+1:n2+m2)
    ])

    v_inds = Dict{String,UnitRange{Int64}}([
        ("x1", 1:n1),
        ("z", n1+1:n1+nz),
        ("v" => 1:n1+nz),
        ("p", n1+nz+1:n1+nz+np),
        ("x", 1:nx),
        ("x2", n1+1:nx),
        ("λ", nx+1:nx+m2)
    ])

    Γ_inds = Dict{String,UnitRange{Int64}}([
        ("G" => 1:m1),
        ("hl" => m1+1:m1+nz),
        ("zl" => m1+nz+1:m1+2*nz),
        ("hu" => m1+2*nz+1:m1+3*nz),
        ("zu" => m1+3*nz+1:m1+4*nz)
    ])

    θ_inds = Dict{String,UnitRange{Int64}}([
        ("v" => 1:n),
        ("Λ" => n+1:n+m),
        ("r" => n+m+1:n+2*m),
        ("p" => n+2*m+1:n+2*m+np),
        ("θ" => 1:n+2*m),
        ("z", n1+1:n1+nz),
        ("ΛG" => n+1:n+m1),
        ("Λhl" => n+m1+1:n+m1+nz),
        ("Λzl" => n+m1+nz+1:n+m1+2*nz),
        ("Λhu" => n+m1+2*nz+1:n+m1+3*nz),
        ("Λzu" => n+m1+3*nz+1:n+m1+4*nz),
        ("rG" => n+m+1:n+m+m1),
        ("rhl" => n+m+m1+1:n+m+m1+nz),
        ("rzl" => n+m+m1+nz+1:n+m+m1+2*nz),
        ("rhu" => n+m+m1+2*nz+1:n+m+m1+3*nz),
        ("rzu" => n+m+m1+3*nz+1:n+m+m1+4*nz),
    ])

    inds = (; x=x_inds, z=z_inds, v=v_inds, Γ=Γ_inds, θ=θ_inds)
end


