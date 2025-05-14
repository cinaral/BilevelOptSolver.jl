%cd("projects/BilevelOptSolver.jl/dataset")
addpath(genpath("./Examples"))
addpath(genpath("./Examples/Linear"))
addpath(genpath("./Examples/Nonlinear"))
addpath(genpath("./Examples/Simple"))

num_problems = 1

function [] = print_views(fileID, nx, ny)
for j = 1:nx
	fprintf(fileID, strcat("\tx", string(j), " = @view xy[", string(j), "]\n"));
end
for j = 1:ny
	fprintf(fileID, strcat("\ty", string(j), " = @view xy[n1+", string(j), "]\n"));
end
end

function [] = print_body(fileID, fun_name, expr, nx, ny, is_arr_out)
if nargin < 6
	is_arr_out = false;
end

fprintf(fileID, strcat("function ", fun_name, "(xy)\n"));
print_views(fileID, nx, ny)
fprintf(fileID, "\t"); % tab

char_expr = replace(char(expr), '+', '.+');
char_expr = replace(char_expr, '-', '.-');

if is_arr_out & isscalar(expr) % we expect an array of length 1
	fprintf(fileID, "[");
end
fprintf(fileID, char_expr);
if is_arr_out & isscalar(expr)
	fprintf(fileID, ",]");
end
fprintf(fileID, "\nend\n\n");
end

for i = 1:num_problems
	[probname, dim, xy, Ff] = InfomAllExamp(i);
	x = sym("x", [dim(1) 1], "real");
	y = sym("y", [dim(2) 1], "real");
	fun = str2func(probname);
	F = fun(x, y, "F");
	G = fun(x, y, "F");
	f = fun(x, y, "f");
	g = fun(x, y, "g");

	fileID = fopen(strcat(probname,'.jl'),'w');
	%n?::Int64 = 5
	%n?::Int64 = 5
	%n::Int64 = n? + n?

	fprintf(fileID, strcat("n1::Int64 = ", string(dim(1))));
	fprintf(fileID, "\n");
	fprintf(fileID, strcat("n2::Int64 = ", string(dim(2))));
	fprintf(fileID, "\n");
	fprintf(fileID, "n::Int64 = n1 + n2\n");
	fprintf(fileID, "\n");

	print_body(fileID, "F", F, dim(1), dim(2))
	print_body(fileID, "G", -G, dim(1), dim(2), true)
	print_body(fileID, "f", f, dim(1), dim(2))
	print_body(fileID, "g", -g, dim(1), dim(2), true)

	fprintf(fileID, strcat("xy_init = [", num2str(reshape(xy', 1, [])), "]'\n"));
	fprintf(fileID, strcat("Ff_optimal = [", num2str(Ff), "]\n"));

	fclose(fileID);
end

