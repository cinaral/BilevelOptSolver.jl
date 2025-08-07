""" 
[ct_1982_01](https://basblsolver.github.io/BASBLib/LP-LP/ct_1982_01)
"""
function ct_1982_01()
    n1::Int64 = 2
    n2::Int64 = 6

    function F(xy)
        x = @view xy[1:n1]
        y = @view xy[n1+1:n1+n2]
        4 * y[1] - 4 * x[2] - 8 * x[1] - 40 * y[2] - 4 * y[3]
    end

    function G(xy)
        x = @view xy[1:n1]
        y = @view xy[n1+1:n1+n2]
        [
            x[1]
            10.0 - x[1]
            x[2]
            10.0 - x[2]
            y[1]
            10.0 - y[1]
            y[2]
            10.0 - y[2]
            y[3]
            10.0 - y[3]
            y[4]
            10.0 - y[4]
            y[5]
            10.0 - y[5]
            y[6]
            10.0 - y[6]
        ]
    end

    function f(xy)
        x = @view xy[1:n1]
        y = @view xy[n1+1:n1+n2]
        x[1] + 2 * x[2] + y[1] + y[2] + 2 * y[3]
    end

    function g(xy)
        x = @view xy[1:n1]
        y = @view xy[n1+1:n1+n2]
        [
            y[1] - y[2] - y[3] - y[4] + 1.0
            y[1] - 2 * x[1] - 2 * y[2] + y[3] / 2 - y[5] + 1.0
            y[2] - 2 * y[1] - 2 * x[2] + y[3] / 2 - y[6] + 1.0
            y[1]
            10.0 - y[1]
            y[2]
            10.0 - y[2]
            y[3]
            10.0 - y[3]
            y[4]
            10.0 - y[4]
            y[5]
            10.0 - y[5]
            y[6]
            10.0 - y[6]
        ]
    end

    xy_optimal = [0.0; 0.9; 0.0; 0.6; 0.4; 0.0; 0.0; 0.0]
    Ff_optimal = [-29.2; 3.2]
    name = "ct_1982_01"

    (; n1, n2, F, G, f, g, xy_optimal, Ff_optimal, name)
end