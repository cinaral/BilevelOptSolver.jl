%cd("projects/BilevelOptSolver.jl/dataset")
addpath(genpath("./Examples"))
addpath(genpath("./Examples/Linear"))
addpath(genpath("./Examples/Nonlinear"))
addpath(genpath("./Examples/Simple"))

num_problems = 173; % max 173

list_fileID = fopen(strcat("./converted/working_problems_list.jl"),'w');
fprintf(list_fileID, strcat("problems = [\n"));

for i = 1:num_problems
	% 79-80 has if else statements, 138 requires partial differential equation toolbox, 173 weird probname
	if i == 79 || i == 80 || i == 138 || i == 173
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
	
	% 36, 49, 50, 51  piecewise
	% 114, 115 does not compile
	% 126, 160 not a valid solution or failed to solve follower NLP
	if ~(i == 36 || i == 49 || i == 50 || i == 51 || i == 114 || i == 115 || i == 126 || i == 160)
		fprintf(list_fileID, strcat("\t""", probname, """\n"));
	end
	
	prob_fileID = fopen(strcat("./converted/", string(i), "_", probname,'.jl'),'w');
	
	fprintf(prob_fileID, strcat("function ", probname, "()\n"));
	fprintf(prob_fileID, strcat("\tn1::Int64 = ", string(dim(1))));
	fprintf(prob_fileID, "\n");
	fprintf(prob_fileID, strcat("\tn2::Int64 = ", string(dim(2))));
	fprintf(prob_fileID, "\n\n");
	
	print_body(prob_fileID, "F", F, dim(1), dim(2))
	print_body(prob_fileID, "G", -G, dim(1), dim(2), true)
	print_body(prob_fileID, "f", f, dim(1), dim(2))
	print_body(prob_fileID, "g", -g, dim(1), dim(2), true)
	
	fprintf(prob_fileID, strcat("\txy_init = Float64", string(mat2str(xy)), "\n"));
	fprintf(prob_fileID, strcat("\tFf_optimal = Float64", string(mat2str(Ff')), "\n\n"));
	
	fprintf(prob_fileID, strcat("\t(; n1, n2, F, G, f, g, xy_init, Ff_optimal)\n"));
	fprintf(prob_fileID, strcat("end ", "\n\n"));
	
	fclose(prob_fileID);
end

fprintf(list_fileID, strcat("]\n"));
fclose(list_fileID);

function [] = print_views(fileID, nx, ny)
fprintf(fileID, strcat("\t\tx = @view xy[1:n1]\n"));
fprintf(fileID, strcat("\t\ty = @view xy[n1+1:n1+n2]\n"));
end


function [] = print_body(fileID, fun_name, expr, nx, ny, is_arr_out)
if nargin < 6
	is_arr_out = false;
end

fprintf(fileID, strcat("\tfunction ", fun_name, "(xy)\n"));
print_views(fileID, nx, ny)
fprintf(fileID, "\t"); % tab

% replace x1, y2... etc with x[1], y[2]...
char_expr = char(expr);
for j = nx:-1:1
	char_expr = replace(char_expr, strcat("x", string(j)), strcat("x[", string(j), "]"));
end

for j = ny:-1:1
	char_expr = replace(char_expr, strcat("y", string(j)), strcat("y[", string(j), "]"));
end

if is_arr_out && length(expr) < 2 % we expect an array of length 1
	fprintf(fileID, "\t[");
end
fprintf(fileID, strcat("\t", char_expr));

if is_arr_out && length(expr) == 1
	fprintf(fileID, ";");
end
if is_arr_out && length(expr) < 2
	fprintf(fileID, "]");
end

fprintf(fileID, "\n\tend\n\n");
end