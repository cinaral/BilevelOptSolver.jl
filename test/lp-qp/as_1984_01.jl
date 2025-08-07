""" 
[as_1984_01](https://basblsolver.github.io/BASBLib/LP-QP/as_1984_01)
"""
function as_1984_01()
    n1::Int64 = 2
    n2::Int64 = 2

    function F(xy)
        x = @view xy[1:n1]
        y = @view xy[n1+1:n1+n2]
        2 * x[1] + 2 * x[2] - 3 * y[1] - 3 * y[2] - 60.0
    end

    function G(xy)
        x = @view xy[1:n1]
        y = @view xy[n1+1:n1+n2]
        [
            2 * y[2] - x[2] - y[1] - x[1] + 40.0
            50.0 - x[1]; 50.0 - x[2];
            x[1]
            50.0 - x[1]
            x[2]
            50.0 - x[2]
        ]
    end

    function f(xy)
        x = @view xy[1:n1]
        y = @view xy[n1+1:n1+n2]
        (y[1] - x[1] + 20.0)^2 + (y[2] - x[2] + 20.0)^2
    end

    function g(xy)
        x = @view xy[1:n1]
        y = @view xy[n1+1:n1+n2]
        [
            x[1] - 2 * y[1] - 10.0
            x[2] - 2 * y[2] - 10.0
            y[1] + 10.0
            20.0 - y[1]
            y[2] + 10.0
            20.0 - y[2]
        ]
    end

    xy_optimal = [25.0, 30.0, 5.0, 10.0]
    Ff_optimal = [5.0; 0.0]
    name = "as_1984_01"

    (; n1, n2, F, G, f, g, xy_optimal, Ff_optimal, name)
end
