n1::Int64 = 2
n2::Int64 = 2
n::Int64 = n1 + n2

function F(xy)
	x1 = @view xy[1]
	x2 = @view xy[2]
	y1 = @view xy[n1+1]
	y2 = @view xy[n1+2]
	2*x1 .+ 2*x2 .- 3*y1 .- 3*y2 .- 60
end

function G(xy)
	x1 = @view xy[1]
	x2 = @view xy[2]
	y1 = @view xy[n1+1]
	y2 = @view xy[n1+2]
	[3*y1 .- 2*x2 .- 2*x1 .+ 3*y2 .+ 60,]
end

function f(xy)
	x1 = @view xy[1]
	x2 = @view xy[2]
	y1 = @view xy[n1+1]
	y2 = @view xy[n1+2]
	(y1 .- x1 .+ 20)^2 .+ (y2 .- x2 .+ 20)^2
end

function g(xy)
	x1 = @view xy[1]
	x2 = @view xy[2]
	y1 = @view xy[n1+1]
	y2 = @view xy[n1+2]
	[x1 .- 2*y1 .- 10; x2 .- 2*y2 .- 10; y1 .+ 10; y2 .+ 10; 20 .- y1; 20 .- y2]
end

xy_init = [10  10  20  20]'
Ff_optimal = [5  0  1]
