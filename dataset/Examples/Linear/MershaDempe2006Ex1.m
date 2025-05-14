function w=MershaDempe2006Ex1(x,y,keyf,keyxy)
% This file provides all functions defining MershaDempe2006Ex1 problem 
% and their first and second order derivatives.
% [dim_x dim_y dim_G dim_g] = [1 1 1 5]  

if nargin<4 || isempty(keyxy)
    switch keyf
    case 'F'; w = x-8*y;
    case 'G'; w = -x;    
    case 'f'; w = y;  
    case 'g'; w = [5*x-2*y-33; -x-2*y+9; -7*x+3*y-5; x+y-15; -y];
    end    
else
    switch keyf
    case 'F'
        switch keyxy
        case 'x' ; w = 1;      
        case 'y' ; w = -8;     
        case 'xx'; w = 0;
        case 'xy'; w = 0;
        case 'yy'; w = 0;
        end 
    case 'G'  
        switch keyxy
        case 'x' ; w = -1;    
        case 'y' ; w = 0;         
        case 'xx'; w = 0;
        case 'xy'; w = 0;
        case 'yy'; w = 0;
        end           
	case 'f'   
        switch keyxy
        case 'x' ; w = 0;    
        case 'y' ; w = 1;         
        case 'xx'; w = 0;
        case 'xy'; w = 0;
        case 'yy'; w = 0;
        end           
	case 'g'   
        switch keyxy
        case 'x' ; w = [5; -1;-7; 1;0];
        case 'y' ; w = [-2; -2; 3; 1; -1];             
        case 'xx'; w = zeros(5,1);
        case 'xy'; w = zeros(5,1);
        case 'yy'; w = zeros(5,1);
        end        
   end   
end

end




