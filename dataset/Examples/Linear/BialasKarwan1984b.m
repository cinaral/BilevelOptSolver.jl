function w=BialasKarwan1984b(x,y,keyf,keyxy)
% This file provides all functions defining BialasKarwan84b problem 
% and their first and second order derivatives.
% [dim_x dim_y dim_G dim_g] = [1 1 1 6]  

if nargin<4 || isempty(keyxy)
    switch keyf
    case 'F'; w = -y;
    case 'G'; w = -x;    
    case 'f'; w = y;  
    case 'g'; w = [-1; 1; -1; 1; 2; 0]*x+[-2; 2; 2; -2; -1; -1]*y+...
                  [10; -38; -18; -6; -21; 0];
    end    
else
    switch keyf
    case 'F'
        switch keyxy
        case 'x' ; w = 0;      
        case 'y' ; w = -1;     
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
        case 'x' ; w = [-1; 1; -1; 1; 2; 0];
        case 'y' ; w = [-2; 2; 2; -2; -1; -1];          
        case 'xx'; w = zeros(6,1);
        case 'xy'; w = zeros(6,1);
        case 'yy'; w = zeros(6,1);
        end        
   end   
end

end




