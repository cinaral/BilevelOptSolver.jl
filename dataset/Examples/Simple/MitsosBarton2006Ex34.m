function w=MitsosBarton2006Ex34(x,y,keyf,keyxy)
% This file provides all functions defining ex34 problem 
% and their first and second order derivatives.
% [dim_x dim_y dim_G dim_g] = [1 1 2 2]   
 
if nargin<4 || isempty(keyxy)
    switch keyf
    case 'F'; w = y;
    case 'G'; w = [-y-0.5; y-1];    
    case 'f'; w = -y^2;  
    case 'g'; w = [-y-0.5; y-1];
    end    
else
    switch keyf
    case 'F'
        switch keyxy
        case 'x' ; w = 0;         
        case 'y' ; w = 1;  
        case 'xx'; w = 0;
        case 'xy'; w = 0;
        case 'yy'; w = 0;
        end 
    case 'G'  
        switch keyxy
        case 'x' ; w = [0; 0];    
        case 'y' ; w = [-1; 1];         
        case 'xx'; w = [0; 0];
        case 'xy'; w = [0; 0];
        case 'yy'; w = [0; 0];      
        end        
	case 'f'   
        switch keyxy
        case 'x' ; w = 0;         
        case 'y' ; w = -2*y;  
        case 'xx'; w = 0;
        case 'xy'; w = 0;
        case 'yy'; w = -2;
        end 
	case 'g'   
        switch keyxy
        case 'x' ; w = [0; 0];    
        case 'y' ; w = [-1; 1];         
        case 'xx'; w = [0; 0]; 
        case 'xy'; w = [0; 0]; 
        case 'yy'; w = [0; 0];       
        end        
   end   
end

end




