
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
    fol_nlp
    # follower MCP: z := [x₂; λ] s.t. h ⟂ zu ≥ z ≥ zl
    fol_mcp
    # SBOPi NLP: min F(v) s.t. Γ := [G; h; -h; z; -z] ≥ Γlᵢ 
    sbop_nlp
    # SBOPi MCP: θ := [v; Λ] s.t. Φ ⟂ θu ≥ θ ≥ θl
    sbop_mcp
    # high-point NLP: min F(x) s.t. [G(x); g(x)] ≥ 0
    hp_nlp
    # corresponding MCP to high-point NLP
    hp_mcp
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
    #vp_sym = Num[v_sym; p_sym]
    F_sym = F(xp_sym)
    G_sym = G(xp_sym)
    f_sym = f(xp_sym)
    g_sym = g(xp_sym)

    # we define these index dictionaries for convenience
    inds = define_index_dicts(n1, n2, m1, m2, nx, nz, n, m)

    # follower NLP: min f(x) s.t. g(x) ≥ 0
    fol_nlp, fol_nlp_sym = convert_to_NLP(n2, m2, x2_sym, f_sym, g_sym, λ_sym, x1_sym, p_sym)

    # follower MCP: z := [x₂; λ] s.t. h ⟂ zu ≥ z ≥ zl
    fol_mcp, fol_mcp_sym = convert_NLP_to_MCP(fol_nlp.n, fol_nlp.m, fol_nlp_sym.x, fol_nlp_sym.g, fol_nlp_sym.λ, fol_nlp_sym.of, fol_nlp_sym.∇ₓL, fol_nlp_sym.q, fol_nlp_sym.p)

    # SBOPi NLP: min F(v) s.t. Γ := [G; h; z; -h; -z] ≥ Γlᵢ 
    # we append [-h; -z] to keep the corresponding Λ ≥ 0 without the need for upper bounds
    Γ_sym = [G_sym; fol_mcp_sym.h; fol_mcp_sym.z; -fol_mcp_sym.h; -fol_mcp_sym.z]
    @assert(m == length(Γ_sym))
    Λ_sym = Symbolics.@variables(Λ[1:m])[1] |> Symbolics.scalarize
    sbop_nlp, sbop_nlp_sym = convert_to_NLP(n, m, v_sym, F_sym, Γ_sym, Λ_sym, [], p_sym)

    # SBOPi MCP: θ := [v; Λ; r] s.t. Φ ⟂ θu ≥ θ ≥ θl
    # in order to conveniently specify the bounds the Γ := [G; h; z; -h; -z], we introduce its slack, r ⟂ Λ
    # this way v and Λ is free, and by default r ≥ 0, but r bounds will be overwritten later
    # TODO 2025-07-14: this is a few more problem variables than strictly needed, but we wish to keep Λ identical to the NLP interpretation for the purposes of checking optimality conditions
    sbop_mcp, sbop_mcp_sym = convert_NLP_to_MCP(sbop_nlp.n, sbop_nlp.m, sbop_nlp_sym.x, sbop_nlp_sym.g, sbop_nlp_sym.λ, sbop_nlp_sym.of, sbop_nlp_sym.∇ₓL, sbop_nlp_sym.q, sbop_nlp_sym.p; xl=fill(-Inf, n), xu=fill(Inf, n), do_use_slacks=true)

    # high-point NLP: min F(x) s.t. [G(x); g(x)] ≥ 0
    Gg_sym = [G_sym; g_sym]
    Λ₃_sym = Symbolics.@variables(Λ[1:m1+m2])[1] |> Symbolics.scalarize
    hp_nlp, hp_nlp_sym = convert_to_NLP(n1 + n2, m1 + m2, x_sym, F_sym, Gg_sym, Λ₃_sym, [], p_sym)

    hp_mcp, hp_mcp_sym = convert_NLP_to_MCP(hp_nlp.n, hp_nlp.m, hp_nlp_sym.x, hp_nlp_sym.g, hp_nlp_sym.λ, hp_nlp_sym.of, hp_nlp_sym.∇ₓL, hp_nlp_sym.q, hp_nlp_sym.p)

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
        fol_nlp,
        fol_mcp,
        sbop_nlp,
        sbop_mcp,
        hp_nlp,
        hp_mcp,
        inds
    )

    syms = (; fol_nlp=fol_nlp_sym, fol_mcp=fol_mcp_sym, sbop_nlp=sbop_nlp_sym, sbop_mcp=sbop_mcp_sym, hp_nlp=hp_nlp_sym, hp_mcp=hp_mcp_sym)

    (; bop, syms)
end


"""
Convert Julia functions into an NLP using symbolics.
```
    min f(x) 
    s.t. g(x) ≥ 0
```
- x ∈ Rⁿ, q and p are optional parameters
- f([q;x;p]) ∈ R
- g([q;x;p]) ∈ Rᵐ
- f([q;x;p]) and g([q;x;p]) ≥ 0 are a twice differentiable 
"""
function convert_to_NLP(n, f, g; nq=0, np=0, xl=zeros(n), xu=fill(Inf, n))
    qxp_dummy = zeros(nq + n + np)
    m = length(g(qxp_dummy))
    x_sym = Symbolics.@variables(x[1:n])[1] |> Symbolics.scalarize
    λ_sym = Symbolics.@variables(λ[1:m])[1] |> Symbolics.scalarize
    q_sym = Symbolics.@variables(q[1:nq])[1] |> Symbolics.scalarize
    p_sym = Symbolics.@variables(p[1:np])[1] |> Symbolics.scalarize

    qxp_sym = Num[q_sym; x_sym; p_sym]
    f_sym = f(qxp_sym)
    g_sym = g(qxp_sym)

    nlp, syms = convert_to_NLP(n, m, x_sym, f_sym, g_sym, λ_sym, q_sym, p_sym; xl, xu)
    (; nlp, syms)
end

function convert_to_NLP(n, m, x_sym, f_sym, g_sym, λ_sym, q_sym, p_sym; xl=fill(-Inf, n), xu=fill(Inf, n))
    qxp_sym = Num[q_sym; x_sym; p_sym]
    gl = zeros(m)
    gu = fill(Inf, m)

    f = Symbolics.build_function(f_sym, qxp_sym; expression=Val{false})
    g! = Symbolics.build_function(g_sym, qxp_sym; expression=Val{false})[2]
    ∇ₓf_sym = Symbolics.gradient(f_sym, x_sym)
    ∇ₓf! = Symbolics.build_function(∇ₓf_sym, qxp_sym; expression=Val{false})[2]

    # ∇ₓg
    ∇ₓg_sym = Symbolics.sparsejacobian(g_sym, x_sym)
    ∇ₓg_size = size(∇ₓg_sym)
    (∇ₓg_rows, ∇ₓg_cols, ∇ₓg_vals_sym) = SparseArrays.findnz(∇ₓg_sym)
    ∇ₓg_vals! = Symbolics.build_function(∇ₓg_vals_sym, qxp_sym; expression=Val{false})[2]

    # ∇²ₓL = of * f(x) - λ' g(x)
    of_sym = Symbolics.@variables(of)[1] # objective factor
    if isempty(λ_sym)
        L_sym = of_sym * f_sym
    else
        L_sym = of_sym * f_sym - g_sym' * λ_sym
    end
    ∇ₓL_sym = Symbolics.gradient(L_sym, x_sym)
    ∇²ₓL_sym = Symbolics.sparsejacobian(∇ₓL_sym, x_sym)
    ∇²ₓL_size = size(∇²ₓL_sym)
    (∇²ₓL_rows, ∇²ₓL_cols, ∇²ₓL_vals_sym) = SparseArrays.findnz(∇²ₓL_sym)
    ∇²ₓL_vals! = Symbolics.build_function(∇²ₓL_vals_sym, qxp_sym, λ_sym, of_sym; expression=Val{false})[2]

    nlp = (; n, m, xl, xu, gl, gu, f, g!, ∇ₓf!, ∇ₓg_rows, ∇ₓg_cols, ∇ₓg_vals!, ∇ₓg_size, ∇²ₓL_rows, ∇²ₓL_cols, ∇²ₓL_vals!, ∇²ₓL_size)
    syms = (; x=x_sym, f=f_sym, g=g_sym, ∇ₓf=∇ₓf_sym, ∇ₓg=∇ₓg_sym, λ=λ_sym, of=of_sym, ∇ₓL=∇ₓL_sym, ∇²ₓL=∇²ₓL_sym, q=q_sym, p=p_sym)

    # force compilation
    q = rand(length(q_sym))
    x = rand(n)
    p = rand(length(p_sym))
    λ = rand(m)
    g = zeros(m)
    ∇ₓf = zeros(n)
    ∇ₓg = sparse(∇ₓg_rows, ∇ₓg_cols, zeros(length(∇ₓg_rows)), ∇ₓg_size[1], ∇ₓg_size[2])
    ∇²ₓL = sparse(∇²ₓL_rows, ∇²ₓL_cols, zeros(length(∇²ₓL_rows)), ∇²ₓL_size[1], ∇²ₓL_size[2])
    nlp.f([q; x; p])
    nlp.g!(g, [q; x; p])
    nlp.∇ₓf!(∇ₓf, [q; x; p])
    nlp.∇ₓg_vals!(∇ₓg.nzval, [q; x; p])
    nlp.∇²ₓL_vals!(∇²ₓL.nzval, [q; x; p], λ, 1.0)

    (; nlp, syms)
end

"""
Convert an NLP problem's first-order KKT conditions into an MCP problem.
```
MCP:  h(z) ⟂ zl ≤ z ≤ zu 
From this NLP:
    min f(x) 
    s.t. g(x) ≥ 0   
    KKT conditions:
        L := f(x) - λ' g(x)
        ∇ₓL = 0  ⟂ x free
        g(x) ≥ 0 ⟂ λ ≥ 0 
        Let z := [x; λ] and h := [∇ₓL; g(x)]
```
- z ∈ Rⁿ⁺ᵐ, q and p are optional parameters
- h([q;z;p]) ∈ Rⁿ⁺ᵐ 
"""
function convert_NLP_to_MCP(n, m, x_sym, g_sym, λ_sym, of_sym, ∇ₓL_sym, q_sym, p_sym; xl=fill(-Inf, n), xu=fill(Inf, n), do_use_slacks=false)
    ∇ₓL_sym = substitute(∇ₓL_sym, Dict([of_sym => 1.0]))

    if do_use_slacks
        nz = n + 2 * m
        r_sym = Symbolics.@variables(r[1:m])[1] |> Symbolics.scalarize
        zl = [xl; fill(-Inf, m); zeros(m)]
        zu = [xu; fill(Inf, 2 * m)]
        h_sym = [∇ₓL_sym; g_sym .- r_sym; λ_sym]
    else
        nz = n + m
        zl = [xl; zeros(m)]
        zu = [xu; fill(Inf, m)]
        r_sym = []
        h_sym = [∇ₓL_sym; g_sym]
    end
    z_sym = [x_sym; λ_sym; r_sym] # z := [x₂; λ]

    @assert(nz == length(zl), "nz and zl dim mismatch")
    @assert(nz == length(zu), "nz and zu dim mismatch")
    @assert(nz == length(h_sym), "nz and h dim mismatch")

    qzp_sym = Num[q_sym; z_sym; p_sym]
    h! = Symbolics.build_function(h_sym, qzp_sym; expression=Val(false))[2]

    # ∇_z h,
    ∇h_sym = Symbolics.sparsejacobian(h_sym, z_sym) # Jacobian of h wrt z
    ∇h_size = size(∇h_sym)
    (∇h_rows, ∇h_cols, ∇h_vals_sym) = SparseArrays.findnz(∇h_sym)
    ∇h_vals! = Symbolics.build_function(∇h_vals_sym, qzp_sym; expression=Val{false})[2]

    mcp = (; n=nz, zl, zu, h!, ∇h_rows, ∇h_cols, ∇h_vals!, ∇h_size)
    syms = (; z=z_sym, r=r_sym, h=h_sym, ∇h=∇h_sym)

    # force compilation
    q = rand(length(q_sym))
    z = rand(nz)
    p = rand(length(p_sym))
    h = zeros(nz)
    ∇h = sparse(∇h_rows, ∇h_cols, zeros(length(∇h_rows)), ∇h_size[1], ∇h_size[2])
    mcp.h!(h, [q; z; p])
    mcp.∇h_vals!(∇h.nzval, [q; z; p])

    (; mcp, syms)
end

"""
    x := [x1; x2] and optionally [x1; x2; p] to let F, G, f, g functions accept paramaters
    z := [x2; λ]
    v := [x1; z] and optionally [x1; z; p] to let F, G, f, g functions accept paramaters
    Γ := [G; h; -h; z; -z] 
    θ := [v; Λ] where Λ'Γ=0 (so param and optionally [v; Λ; p] to let F, G, f, g functions accept paramaters
"""
function define_index_dicts(n1, n2, m1, m2, nx, nz, n, m)
    x_inds = Dict{String,UnitRange{Int64}}([
        ("x1", 1:n1),
        ("x2", n1+1:n1+n2),
        ("x", 1:n1+n2)
    ])

    z_inds = Dict{String,UnitRange{Int64}}([
        ("x2", 1:n2)
        ("λ", n2+1:n2+m2)
    ])

    v_inds = Dict{String,UnitRange{Int64}}([
        ("x1", 1:n1),
        ("z", n1+1:n1+nz),
        ("v" => 1:n1+nz),
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
