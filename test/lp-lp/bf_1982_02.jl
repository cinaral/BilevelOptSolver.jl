""" 
[bf_1982_02](https://basblsolver.github.io/BASBLib/LP-LP/bf_1982_02)
"""
function bf_1982_02()
    n1::Int64 = 2
    n2::Int64 = 2

    function F(xy)
        x = @view xy[1:n1]
        y = @view xy[n1+1:n1+n2]
        x[2] - 2 * x[1] + y[1] / 2
    end

    function G(xy)
        x = @view xy[1:n1]
        y = @view xy[n1+1:n1+n2]
        [
            x[1]
            10.0 - x[1]
            x[2]
            10.0 - x[2]
            # fails without this:
            y[1]
            10.0 - y[1]
            y[2]
            10.0 - y[2]
        ]
    end

    function f(xy)
        x = @view xy[1:n1]
        y = @view xy[n1+1:n1+n2]
        x[1] + x[2] - 4 * y[1] + y[2]
    end

    function g(xy)
        x = @view xy[1:n1]
        y = @view xy[n1+1:n1+n2]
        [
            2 * x[1] - y[1] + y[2] - 5 / 2
            3 * x[2] - x[1] - y[2] + 2
            2 - x[2] - x[1]
            y[1]
            10.0 - y[1]
            y[2]
            10.0 - y[2]
        ]
    end

    xy_optimal = [2.0; 0.0; 1.5; 0.0]
    Ff_optimal = [-3.25; -4.0]
    name = "bf_1982_02"

    (; n1, n2, F, G, f, g, xy_optimal, Ff_optimal, name)
end