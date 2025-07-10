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
    ∇ₓ₂g_size::Tuple{Int, Int}
    ∇ₓ₂g_rows::Vector{Int}
    ∇ₓ₂g_cols::Vector{Int}
    ∇ₓ₂g_vals!::Function # ∇ₓ₂g_vals!(out, x)
    # ∇²ₓ₂L₂, L₂ = of * f(x) - λ' g(x)
    ∇²ₓ₂L₂_size::Tuple{Int, Int}
    ∇²ₓ₂L₂_rows::Vector{Int}
    ∇²ₓ₂L₂_cols::Vector{Int}
    ∇²ₓ₂L₂_vals!::Function # ∇²ₓ₂L_vals!(out, v, λ, of)
    # follower MCP: z := [x₂; λ] s.t. h ⟂ zu ≥ z ≥ zl
    zl₀::Vector{Float64} # default z bounds 
    zu₀::Vector{Float64}
    h!::Function # h!(out, v)
    # ∇_z h
    ∇h_size::Tuple{Int, Int}
    ∇h_rows::Vector{Int}
    ∇h_cols::Vector{Int}
    ∇h_vals!::Function # ∇h_vals!(out, v)
    # SBOPi NLP: min F(v) s.t. Γ := [G; h; -h; z; -z] ≥ Γlᵢ 
    Fv::Function # F(v)
    Γ!::Function # Γ!(out, v)
    ∇ᵥF!::Function # Γ!(out, v) 
    # ∇ᵥΓ
    ∇ᵥΓ_size::Tuple{Int, Int}
    ∇ᵥΓ_rows::Vector{Int}
    ∇ᵥΓ_cols::Vector{Int}
    ∇ᵥΓ_vals!::Function # ∇ᵥΓ_vals!(out, v)
    # ∇²ᵥL₁, L₁ = of * F(v) - Λ' Γ(v)
    ∇²ᵥL₁_size::Tuple{Int, Int}
    ∇²ᵥL₁_rows::Vector{Int}
    ∇²ᵥL₁_cols::Vector{Int}
    ∇²ᵥL₁_vals!::Function # ∇²ᵥL₁_vals!(out, v, Λ, of)
    # SBOPi MCP: θ := [v; Λ] s.t. Φ ⟂ θu ≥ θ ≥ θl
    θl₀::Vector{Float64} # default θ bounds
    θu₀::Vector{Float64}
    Φ!::Function # Φ!(out, v, Λ)
    # ∇_θ Φ
    ∇Φ_size::Tuple{Int, Int}
    ∇Φ_rows::Vector{Int}
    ∇Φ_cols::Vector{Int}
    ∇Φ_vals!::Function # ∇Φ_vals!(out, θ)
    # high-point NLP: min F(x) s.t. [G(x); g(x)] ≥ 0
    Gg!::Function # Gg!(out, x) 
    ∇ₓF!::Function # ∇ₓF!(out, x)
    # ∇ₓGg
    ∇ₓGg_size::Tuple{Int, Int}
    ∇ₓGg_rows::Vector{Int}
    ∇ₓGg_cols::Vector{Int}
    ∇ₓGg_vals!::Function # ∇ₓGg_vals!(out, x)
    # ∇²ₓL₃, L₃ = of * F(x) - Λ₃' Gg(x)
    ∇²ₓL₃_size::Tuple{Int, Int}
    ∇²ₓL₃_rows::Vector{Int}
    ∇²ₓL₃_cols::Vector{Int}
    ∇²ₓL₃_vals!::Function # ∇²ₓ₃L_vals!(out, x, Λ₃, of)
    # for convenience
    x_inds::Dict{String,UnitRange{Int64}} # indexes x := [x₁; x₂]
    z_inds::Dict{String,UnitRange{Int64}} # indexes z := [x₂; λ]
    v_inds::Dict{String,UnitRange{Int64}} # indexes v := [x₁; z]
    Γ_inds::Dict{String,UnitRange{Int64}} # indexes Γ := [G; h; -h; z; -z]
    θ_inds::Dict{String,UnitRange{Int64}} # indexes θ := [v; Λ]
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
    nθ = n + m

    x1_sym = Symbolics.@variables(x1[1:n1])[1] |> Symbolics.scalarize
    x2_sym = Symbolics.@variables(x2[1:n2])[1] |> Symbolics.scalarize
    λ_sym = Symbolics.@variables(λ[1:m2])[1] |> Symbolics.scalarize
    p_sym = Symbolics.@variables(p[1:np])[1] |> Symbolics.scalarize
    x_sym = [x1_sym; x2_sym] # x := [x₁; x₂]
    z_sym = [x2_sym; λ_sym] # z := [x₂; λ]
    v_sym = [x1_sym; z_sym] # z := [x₁; z]
    @assert(nx == length(x_sym))
    @assert(nz == length(z_sym))
    @assert(n == length(v_sym))
    F_sym = F(x_sym)
    G_sym = G(x_sym)
    f_sym = f(x_sym)
    g_sym = g(x_sym)

    xp_sym = [x_sym; p_sym]
    vp_sym = [v_sym; p_sym]

    # we define these for convenience
    x_inds = Dict([
        ("x1", 1:n1),
        ("x2", n1+1:n1+n2),
        ("p", n1+n2+1:n1+n2+np)
    ])
    z_inds = Dict([
        ("x2", 1:n2),
        ("λ", n2+1:n2+m2)
    ])
    v_inds = Dict([
        ("x", 1:nx),
        ("x1", 1:n1),
        ("z", n1+1:n1+nz),
        ("x2", n1+1:nx),
        ("λ", nx+1:nx+m2)
    ])
    Γ_inds = Dict([
        ("G" => 1:m1),
        ("hl" => m1+1:m1+nz),
        ("hu" => m1+nz+1:m1+2*nz),
        ("zl" => m1+2*nz+1:m1+3*nz),
        ("zu" => m1+3*nz+1:m1+4*nz)
    ])
    θ_inds = Dict([
        ("v" => 1:n),
        ("Λ" => n+1:n+m),
        ("ΛG" => n+1:n+m1),
        ("Λhl" => n+m1+1:n+m1+nz),
        ("Λhu" => n+m1+nz+1:n+m1+2*nz),
        ("Λzl" => n+m1+2*nz+1:n+m1+3*nz),
        ("Λzu" => n+m1+3*nz+1:n+m1+4*nz)
    ])

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

    # ∇²ₓ₂L₂, L₂ =  of f(x) - λ' g(x)
    of₂_sym = Symbolics.@variables(of2)[1] # objective factor
    if isempty(λ_sym)
        L₂_sym = of₂_sym * f_sym
    else
        L₂_sym = of₂_sym * f_sym - g_sym' * λ_sym
    end
    ∇ₓ₂L₂_sym = Symbolics.gradient(L₂_sym, x2_sym)
    ∇²ₓ₂L₂_sym = Symbolics.sparsejacobian(∇ₓ₂L₂_sym, x2_sym)
    ∇²ₓ₂L₂_size = size(∇²ₓ₂L₂_sym)
    (∇²ₓ₂L₂_rows, ∇²ₓ₂L₂_cols, ∇²ₓ₂L₂_vals_sym) = SparseArrays.findnz(∇²ₓ₂L₂_sym)
    ∇²ₓ₂L₂_vals! = Symbolics.build_function(∇²ₓ₂L₂_vals_sym, xp_sym, λ_sym, of₂_sym; expression=Val{false})[2]

    # follower MCP: z := [x₂; λ] s.t. h ⟂ zu ≥ z ≥ zl
    # by default x₂ is free and λ ≥ 0, but these will be overwritten later
    ∇ₓ₂L₂_sym = substitute(∇ₓ₂L₂_sym, Dict([of₂_sym=>1.]))
    zl₀ = [fill(-Inf, n2); zeros(m2)] # default z lb
    zu₀ = fill(Inf, nz) # default z ub
    if nz > 0
        h_sym = [∇ₓ₂L₂_sym; g_sym]
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

    # SBOPi NLP: min F(v) s.t. Γ := [G; h; -h; z; -z] ≥ Γlᵢ 
    Γ_sym = [G_sym; h_sym; -h_sym; z_sym; -z_sym]
    @assert(m == length(Γ_sym))
    Fv = Symbolics.build_function(F_sym, vp_sym; expression=Val{false})
    Γ! = Symbolics.build_function(Γ_sym, vp_sym; expression=Val(false))[2]
    ∇ᵥF_sym = Symbolics.gradient(F_sym, vp_sym)
    ∇ᵥF! = Symbolics.build_function(∇ᵥF_sym, vp_sym; expression=Val{false})[2]

    # ∇ᵥΓ
    ∇ᵥΓ_sym = Symbolics.sparsejacobian(Γ_sym, v_sym)
    ∇ᵥΓ_size = size(∇ᵥΓ_sym)
    (∇ᵥΓ_rows, ∇ᵥΓ_cols, ∇ᵥΓ_vals_sym) = SparseArrays.findnz(∇ᵥΓ_sym)
    ∇ᵥΓ_vals! = Symbolics.build_function(∇ᵥΓ_vals_sym, vp_sym; expression=Val{false})[2]

    # ∇²ᵥL₁, L₁ = F(v) - Λ' Γ(v)
    Λ_sym = Symbolics.@variables(Λ[1:m])[1] |> Symbolics.scalarize
    of₁_sym = Symbolics.@variables(of1)[1] # objective factor
    if isempty(λ_sym)
        L₁_sym = of₁_sym * F_sym
    else
        L₁_sym = of₁_sym * F_sym - Γ_sym' * Λ_sym
    end
    ∇ᵥL₁_sym = Symbolics.gradient(L₁_sym, v_sym)
    ∇²ᵥL₁_sym = Symbolics.sparsejacobian(∇ᵥL₁_sym, v_sym)
    ∇²ᵥL₁_size = size(∇²ᵥL₁_sym)
    (∇²ᵥL₁_rows, ∇²ᵥL₁_cols, ∇²ᵥL₁_vals_sym) = SparseArrays.findnz(∇²ᵥL₁_sym)
    ∇²ᵥL₁_vals! = Symbolics.build_function(∇²ᵥL₁_vals_sym, vp_sym, Λ_sym, of₁_sym; expression=Val{false})[2]

    # SBOPi MCP: θ := [v; Λ] s.t. Φ ⟂ θu ≥ θ ≥ θl
    # by default v is free and Λ ≥ 0, but these will be overwritten later
    θ_sym = [v_sym; Λ_sym]
    @assert(nθ == length(θ_sym))
    θp_sym = [v_sym; Λ_sym; p_sym]
    θl₀ = [fill(-Inf, n); zeros(m)] # default θ lb
    θl₀[θ_inds["Λhl"]] .= zl₀
    θu₀ = fill(Inf, nθ) # default θ ub
    θu₀[θ_inds["Λhu"]] .= zu₀
    ∇ᵥL₁_sym = substitute(∇ᵥL₁_sym, Dict([of₁_sym=>1.]))
    Φ_sym = [∇ᵥL₁_sym; Γ_sym]
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

    # ∇²ₓL₃, L₃ = F(x) - Λ₃' Gg(x)
    Λ₃_sym = Symbolics.@variables(Λ[1:m1+m2])[1] |> Symbolics.scalarize
    of₃_sym = Symbolics.@variables(of3)[1] # objective factor
    if isempty(Λ₃_sym)
        L₃_sym = of₃_sym * F_sym
    else
        L₃_sym = of₃_sym * F_sym - Gg_sym' * Λ₃_sym
    end
    ∇ₓL₃_sym = Symbolics.gradient(L₃_sym, x_sym)
    ∇²ₓL_sym = Symbolics.sparsejacobian(∇ₓL₃_sym, x_sym)
    ∇²ₓL₃_size = size(∇²ₓL_sym)
    (∇²ₓL₃_rows, ∇²ₓL₃_cols, ∇²ₓL₃_vals_sym) = SparseArrays.findnz(∇²ₓL_sym)
    ∇²ₓL₃_vals! = Symbolics.build_function(∇²ₓL₃_vals_sym, xp_sym, Λ₃_sym, of₃_sym; expression=Val{false})[2]

    BilevelOptProb(
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
        # ∇²ₓ₂L₂, L₂ = f(x) - λ' g(x)
        ∇²ₓ₂L₂_size,
        ∇²ₓ₂L₂_rows,
        ∇²ₓ₂L₂_cols,
        ∇²ₓ₂L₂_vals!,
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
        # ∇²ᵥL₁, L₁ = F(v) - Λ' Γ(v)
        ∇²ᵥL₁_size,
        ∇²ᵥL₁_rows,
        ∇²ᵥL₁_cols,
        ∇²ᵥL₁_vals!,
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
        # ∇²ₓL₃, L₃ = F(x) - Λ₃' Gg(x)
        ∇²ₓL₃_size,
        ∇²ₓL₃_rows,
        ∇²ₓL₃_cols,
        ∇²ₓL₃_vals!,
        # for convenience
        x_inds,
        z_inds,
        v_inds,
        Γ_inds,
        θ_inds
    )
end
