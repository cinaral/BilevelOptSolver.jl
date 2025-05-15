%cd("projects/BilevelOptSolver.jl/dataset")
addpath(genpath("./Examples"))
addpath(genpath("./Examples/Linear"))
addpath(genpath("./Examples/Nonlinear"))
addpath(genpath("./Examples/Simple"))

num_problems = 173 % max 173

list_fileID = fopen(strcat("./converted/problems_list.jl"),'w');
fprintf(list_fileID, strcat("problems = [\n"));

for i = 1:num_problems
	if (i == 79 || i == 80 || i == 138 || i == 173)
		continue
	end
	fprintf('converting %d\n',i)
	[probname, dim, xy, Ff] = InfomAllExamp(i);
	x = sym("x", [dim(1) 1], "real");
	y = sym("y", [dim(2) 1], "real");
	fun = str2func(probname);
	F = fun(x, y, "F");
	G = fun(x, y, "G");
	f = fun(x, y, "f");
	g = fun(x, y, "g");
	
	fprintf(list_fileID, strcat("\t""", probname, """\n"));
	
	prob_fileID = fopen(strcat("./converted/", string(i), "_", probname,'.jl'),'w');
	
	%n?::Int64 = 5
	%n?::Int64 = 5
	%n::Int64 = n? + n?
	
	fprintf(prob_fileID, strcat("function ", probname, "()\n\n"));
	fprintf(prob_fileID, strcat("n1::Int64 = ", string(dim(1))));
	fprintf(prob_fileID, "\n");
	fprintf(prob_fileID, strcat("n2::Int64 = ", string(dim(2))));
	fprintf(prob_fileID, "\n");
	
	print_body(prob_fileID, "F", F, dim(1), dim(2))
	print_body(prob_fileID, "G", -G, dim(1), dim(2), true)
	print_body(prob_fileID, "f", f, dim(1), dim(2))
	print_body(prob_fileID, "g", -g, dim(1), dim(2), true)
	
	fprintf(prob_fileID, strcat("xy_init = Float64", string(mat2str(xy)), "\n"));
	fprintf(prob_fileID, strcat("Ff_optimal = Float64", string(mat2str(Ff')), "\n\n"));
	
	fprintf(prob_fileID, strcat("(; n1, n2, F, G, f, g, xy_init, Ff_optimal)\n\n"));
	fprintf(prob_fileID, strcat("end ", "\n\n"));
	
	fclose(prob_fileID);
end

fprintf(list_fileID, strcat("]\n"));
fclose(list_fileID);

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

char_expr = char(expr);
char_expr = replace(char_expr, '*', ' .* ');
char_expr = replace(char_expr, '+', '.+');
char_expr = replace(char_expr, '-', '.-');
char_expr = replace(char_expr, '^', ' .^');
char_expr = replace(char_expr, 'exp(', 'exp.(');

if is_arr_out & length(expr) < 2 % we expect an array of length 1
	fprintf(fileID, "[");
end
fprintf(fileID, char_expr);
if is_arr_out & length(expr) == 1
	fprintf(fileID, ";");
end
if is_arr_out & length(expr) < 2
	fprintf(fileID, "]");
end
if ~is_arr_out
	fprintf(fileID, "[1]");
end
fprintf(fileID, "\nend\n\n");
end