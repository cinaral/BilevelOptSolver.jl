""" 
[as_1981_01](https://basblsolver.github.io/BASBLib/QP-QP/as_1981_01)
"""
function as_1981_01()
    n1::Int64 = 4
    n2::Int64 = 4

    function F(xy)
        x = @view xy[1:n1]
        y = @view xy[n1+1:n1+n2]
        -(200.0 - y[1] - y[3]) * (y[1] + y[3]) - (160.0 - y[2] - y[4]) * (y[2] + y[4])
    end

    function G(xy)
        x = @view xy[1:n1]
        y = @view xy[n1+1:n1+n2]
        [
            -x[1] - x[2] - x[3] - x[4] + 40.0
            x[1]
            10.0 - x[1]
            x[2]
            5.0 - x[2]
            x[3]
            15.0 - x[3]
            x[4]
            20.0 - x[4]
            y[1]
            20.0 - y[1]
            y[2]
            20.0 - y[2]
            y[3]
            40.0 - y[3]
            y[4]
            40.0 - y[4]
        ]
    end

    function f(xy)
        x = @view xy[1:n1]
        y = @view xy[n1+1:n1+n2]
        (y[1] - 4.0)^2 + (y[2] - 13.0)^2 + (y[3] - 35.0)^2 + (y[4] - 2.0)^2
    end

    function g(xy)
        x = @view xy[1:n1]
        y = @view xy[n1+1:n1+n2]
        [
            -0.4 * y[1] - 0.7 * y[2] + x[1]
            -0.6 * y[1] - 0.3 * y[2] + x[2]
            -0.4 * y[3] - 0.7 * y[4] + x[3]
            -0.6 * y[3] - 0.3 * y[4] + x[4]
            y[1]
            20.0 - y[1]
            y[2]
            20.0 - y[2]
            y[3]
            40.0 - y[3]
            y[4]
            40.0 - y[4]
        ]
    end

    xy_init = [0.1; 0.1; 0.1; 0.1; 0.1; 0.1; 0.1; 0.1]
    xy_optimal = [7.0; 3.0; 12.0; 18.0; 0.0; 10.0; 30.0; 0.0]
    Ff_optimal = [-6600.0; 54.0]
    name = "as_1981_01"

    (; n1, n2, F, G, f, g, xy_init, xy_optimal, Ff_optimal, name)
end