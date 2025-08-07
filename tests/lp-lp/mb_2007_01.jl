""" 
[mb_2007_01](https://basblsolver.github.io/BASBLib/LP-LP/mb_2007_01)

linear f
"""
function mb_2007_01()
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
        -y[1]
    end

    function g(xy)
        x = @view xy[1:n1]
        y = @view xy[n1+1:n1+n2]
        [
            y[1] + 1.0
            1.0 - y[1]
        ]
    end

    xy_optimal = [1.0]
    Ff_optimal = [1.0; -1.0]
    name = "mb_2007_01"

    (; n1, n2, F, G, f, g, xy_optimal, Ff_optimal, name)
end