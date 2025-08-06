using BilevelOptSolver

# 2025-08-06 TODO: if y_init=-0., we will be stuck at a local maximum. we should randomize to get out of there

""" 
[mb_2007_06](https://basblsolver.github.io/BASBLib/LP-NLP/mb_2007_06)

nonconvex f
stationary points of f are y={0}
"""
function mb_2007_06()
    n1::Int64 = 0
    n2::Int64 = 1

    function F(xy)
        x = @view xy[1:n1]
        y = @view xy[n1+1:n1+n2]
        y[1]
    end

    function G(xy)
        x = @view xy[1:n1]
        y = @view xy[n1+1:n1+n2]
        [
            y[1] + 1.0
            1.0 - y[1]
        ]
    end

    function f(xy)
        x = @view xy[1:n1]
        y = @view xy[n1+1:n1+n2]
        y[1]^3
    end

    function g(xy)
        x = @view xy[1:n1]
        y = @view xy[n1+1:n1+n2]
        [
            y[1] + 1.0
            1.0 - y[1]
        ]
    end

    xy_init = [0.1]
    xy_optimal = [-1.0]
    Ff_optimal = [-1.0; -1.0]
    name = "mb_2007_06"

    (; n1, n2, F, G, f, g, xy_init, xy_optimal, Ff_optimal, name)
end
