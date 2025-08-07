""" 
[b_1991_02](https://basblsolver.github.io/BASBLib/LP-QP/b_1991_02)
"""
function b_1991_02()
    n1::Int64 = 1
    n2::Int64 = 2

    function F(xy)
        x = @view xy[1:n1]
        y = @view xy[n1+1:n1+n2]
        x[1] + y[2]
    end

    function G(xy)
        x = @view xy[1:n1]
        y = @view xy[n1+1:n1+n2]
        [
            x[1] - 2.0
            4.0 - x[1]
            #y[1]
            #10.0 - y[1]
            #y[2]
            #10.0 - y[2]
        ]
    end

    function f(xy)
        x = @view xy[1:n1]
        y = @view xy[n1+1:n1+n2]
        2 * y[1] + x[1] * y[2]
    end

    function g(xy)
        x = @view xy[1:n1]
        y = @view xy[n1+1:n1+n2]
        [
            -x[1] + y[1] + y[2] - 4.0
            y[1]
            10.0 - y[1]
            y[2]
            10.0 - y[2]
        ]
    end

    xy_optimal = [2.0; 6.0; 0.0]
    Ff_optimal = [2.0; 12.0]
    name = "b_1991_02"

    (; n1, n2, F, G, f, g, xy_optimal, Ff_optimal, name)
end
