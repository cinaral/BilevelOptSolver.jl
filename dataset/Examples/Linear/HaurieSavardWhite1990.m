function w=HaurieSavardWhite1990(x,y,keyf,keyxy)
% This file provides all functions defining HaurieSavardWhite1990 problem 
% and their first and second order derivatives.
% [dim_x dim_y dim_G dim_g] = [1 1 0 4]  

if nargin<4 || isempty(keyxy)
    switch keyf
    case 'F'; w = x+5*y; 
    case 'G'; w = [];    
    case 'f'; w = -y;  
    case 'g'; w = [-3*x+2*y-6; 3*x+4*y-48; 2*x-5*y-9; -x-y+8];
    end    
else
    switch keyf
    case 'F'
        switch keyxy
        case 'x' ; w = 1;      
        case 'y' ; w = 5;     
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
        case 'y' ; w = -1;         
        case 'xx'; w = 0;
        case 'xy'; w = 0;
        case 'yy'; w = 0;
        end           
	case 'g'   
        switch keyxy
        case 'x' ; w = [-3; 3; 2; -1];
        case 'y' ; w = [2; 4; -5; -1];             
        case 'xx'; w = zeros(4,1);
        case 'xy'; w = zeros(4,1);
        case 'yy'; w = zeros(4,1);
        end        
   end   
end

end




