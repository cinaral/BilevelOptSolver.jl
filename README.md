# BilevelOptSolver.jl

This package solves bilevel optimization problems using a novel LP and NLP/MPC based approach.
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
```

See ```examples/``` folder for usage examples:
```julia
bop = construct_bop(n₁, n₂, F, G, f, g);
sol = solve_bop(bop)
```

## (Optional) Use HSL solvers with IPOPT
If you would like to use HSL, please [obtain a license and download HSL_jll.jl](https://licences.stfc.ac.uk/products/Software/HSL/LibHSL), and set HSL_PATH environment variable to the extracted location.

## (Optional) Solve larger problems with PATH Solver
If you wish to use PATH Solver with problems larger than 300 variables and 2000 non-zeros, you can [obtain a one-year license here](https://pages.cs.wisc.edu/~ferris/path/julia/LICENSE) and set the environment variable "PATH_LICENSE_STRING".