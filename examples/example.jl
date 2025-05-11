using BilevelOptSolver

# Kleinert (2021)
n₁ = 1
n₂ = 1
n = n₁ + n₂

function F(xx)
    x = @view xx[1:n₁]
    y = @view xx[n₁+1:n₁+n₂]
    x[1] + 6y[1]
end

function G(xx)
    x = @view xx[1:n₁]
    y = @view xx[n₁+1:n₁+n₂]
    [
        x[1] - 5y[1] + 12.5,
    ]
end

function f(xx)
    x = @view xx[1:n₁]
    y = @view xx[n₁+1:n₁+n₂]
    -y[1]
end

function g(xx)
    x = @view xx[1:n₁]
    y = @view xx[n₁+1:n₁+n₂]
    [
        2x[1] - y[1]
        -x[1] - y[1] + 6
        -x[1] + 6y[1] + 3
        x[1] + 3y[1] - 3
    ]
end

bop = construct_bop(n₁, n₂, F, G, f, g);
sol = solve_bop(bop)
#sol = solve_bop(bop; x₁_init=[2.4564338234981746], x₂_init=[0.9845259227566776]) # bilevel feasible init
#sol = solve_bop(bop; x₁_init=[3/7], x₂_init=[6/7]) # optimal init

"""
TODO 2025-05-11:  
1. merge into repo and clear commits
2. fix bop type instability DONE
3. fix feas LP
4. test on BOLIP
5. make modular
"""