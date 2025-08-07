""" 
[mb_2007_05](https://basblsolver.github.io/BASBLib/LP-NLP/mb_2007_05)

nonconvex f
stationary points of f are y={-0.5, -0.09375, 0.5}
"""
function mb_2007_05()
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
        16 * y[1]^4 + 2 * y[1]^3 - 8 * y[1]^2 - 1.5 * y[1] + 0.5
    end

    function g(xy)
        x = @view xy[1:n1]
        y = @view xy[n1+1:n1+n2]
        [
            y[1] + 1.0
            1.0 - y[1]
        ]
    end

    xy_init = [-0.09375]
    xy_optimal = [0.5]
    Ff_optimal = [0.5; -1.0]
    name = "mb_2007_05"

    (; n1, n2, F, G, f, g, xy_init, xy_optimal, Ff_optimal, name)
end