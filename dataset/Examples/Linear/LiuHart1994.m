function w=LiuHart1994(x,y,keyf,keyxy)
% This file provides all functions defining LiuHart1994 problem 
% and their first and second order derivatives.
% [dim_x dim_y dim_G dim_g] = [1 1 1 4]  

if nargin<4 || isempty(keyxy)
    switch keyf
    case 'F'; w = -x-3*y; 
    case 'G'; w = -x;        
    case 'f'; w = y;  
    case 'g'; w = [-x+y-3; x+2*y-12; 4*x-y-12; -y];
    end    
else
    switch keyf
    case 'F'
        switch keyxy
        case 'x' ; w = -1;      
        case 'y' ; w = -3;     
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
        case 'x' ; w = [-1; 1; 4;  0];
        case 'y' ; w = [1; 2; -1; -1];             
        case 'xx'; w = zeros(4,1);
        case 'xy'; w = zeros(4,1);
        case 'yy'; w = zeros(4,1);
        end        
   end   
end

end




