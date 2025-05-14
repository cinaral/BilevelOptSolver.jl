function w=TuyEtal2007Ex3(x,y,keyf,keyxy)
% This file provides all functions defining TuyEtal2007Ex3 problem 
% and their first and second order derivatives.
 
% [dim_x dim_y dim_G dim_g] = [10 6 12 13]  

[A,B,C,D]=data();

if nargin<4 || isempty(keyxy)

    switch keyf
    case 'F'; w = [12 -1 -12 13 0 2 0 -5  6 -11]*x -[5 6 4 7 0 0]*y;
    case 'G'; w = [-A*x-B*y + [-30; 134]; -x];    
    case 'f'; w = [3 -2 -3 -3 1 6]*y;  
    case 'g'; w = [-C*x-D*y+[-83; -92; -168; 96; 133; -89; 192]; -y];
    end    
else
    switch keyf
    case 'F'
        switch keyxy
        case 'x' ; w = [12; -1; -12; 13; 0; 2; 0; -5; 6; -11];      
        case 'y' ; w = [-5; -6; -4; -7; 0; 0];     
        case 'xx'; w = zeros(10,10);
        case 'xy'; w = zeros(6,10);
        case 'yy'; w = zeros(6,6);
        end 
    case 'G'  
        switch keyxy
        case 'x' ; w = [-A; -eye(10)];    
        case 'y' ; w = [-B; zeros(10,6)];         
        case 'xx'; w = zeros(120,10);
        case 'xy'; w = zeros(72,10);
        case 'yy'; w = zeros(72,6);
        end           
	case 'f'   
        switch keyxy
        case 'x' ; w = zeros(10,1);    
        case 'y' ; w = [3; -2; -3; -3; 1; 6];         
        case 'xx'; w = zeros(10,10);
        case 'xy'; w = zeros(6,10);
        case 'yy'; w = zeros(6,6);
        case 'yxx'; w = [];
        case 'yxy'; w = [];
        case 'yyy'; w = []; 
        end           
	case 'g'   
        switch keyxy
        case 'x' ; w = [-C; zeros(6,10)];
        case 'y' ; w = [-D; -eye(6)];          
        case 'xx'; w = zeros(130,10);
        case 'xy'; w = zeros(78,10);
        case 'yy'; w = zeros(78,6);
        case 'yxx'; w = [];
        case 'yxy'; w = [];
        case 'yyy'; w = [];            
        end        
   end   
end

end

function [A,B,C,D]=data()
 A = [2 3 -14 2 9 -2 -1 4 0 -2;
    -1 7 -13 0 15 -2 8 4 -4 7];
B = [3 -9 2 8 -1 8; 
    6 2 -6 -2 -8 4];
C = [5 -7 -4 2 -3 9 -9 1 3 -11; 
    -6 5 3 2 -8 -5 -8 3 -7 -3; 
    6 4 -2 0 2 -3 3 -2 -2 -4; 
    -5 -6 0 4 -3 8 -1 0 -2 3;
    -11 11 -4 -5 10 6 -14 7 11 3; 
    -9 12 4 10 -2 -8 -5 11 4 -1; 
    -7 2 6 0 11 -1 2 2 1 2];
D = [-10 9 6 -4 -6 3; 
    5 7 -1 -1 6 -4; 
    -10 -5 -6 4 -3 1; 
    4 3 4 4 -1 -1; 
    10 7 -7 -7 -2 -7;
    -2 5 -10 -1 -4 -5; 
    5 5 6 5 -1 12];
end



