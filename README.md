# BilevelOptSolver.jl

This is an NLP/LP based solver for the bilevel optimization problem:
```
Assume, 
F: Rⁿ¹->R, 
f: Rⁿ²->R, 
G: X->Rᵐ¹, 
g: X->Rᵐ², 
are all continuous and twice differentiable.

Bilevel Optimization Problem (BOP), Standard (Optimistic), let x=[x₁;x₂]:
min     F(x)
x
s.t.    x₁ ∈ Rⁿ¹
        G(x) ≥ 0                                            (BOP)
        x₂ ∈ S := { x₂ : x₂ ∈ arg min   f(x₁, y) 
                                 y   
                                s.t.    y ∈ Rⁿ²
                                        g(x₁, y) ≥ 0 }           
```
See ```examples/``` folder for usage examples:
```julia
bop = construct_bop(n₁, n₂, F, G, f, g);
sol = solve_bop(bop)
```