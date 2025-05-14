function w=ClarkWesterberg1988(x,y,keyf,keyxy)
% This file provides all functions defining ClarkWesterberg88 problem 
% and their first and second order derivatives.
% [dim_x dim_y dim_G dim_g] = [1 1 0 3]  

if nargin<4 || isempty(keyxy)
    switch keyf
    case 'F'; w = x-4*y;   
    case 'G'; w = [];
    case 'f'; w = y;  
    case 'g'; w = [-2*x+y; 2*x+5*y-108; 2*x-3*y+4];
    end    
else
    switch keyf
    case 'F'
        switch keyxy
        case 'x' ; w = 1;      
        case 'y' ; w = -4;     
        case 'xx'; w = 0;
        case 'xy'; w = 0;
        case 'yy'; w = 0;
        end 
    case 'G'   
        switch keyxy
        case 'x' ; w = [];    
        case 'y' ; w = [];         
        case 'xx'; w = [];
        case 'xy'; w = [];
        case 'yy'; w = [];
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
        case 'x' ; w = [-2; 2; 2];
        case 'y' ; w = [1; 5; -3];             
        case 'xx'; w = zeros(3,1);
        case 'xy'; w = zeros(3,1);
        case 'yy'; w = zeros(3,1);
        end        
   end   
end

end




