function w=AnandalinghamWhite1990(x,y,keyf,keyxy)
% This file provides all functions defining AnandalinghamWhite1990 problem 
% and their first and second order derivatives.
% [dim_x dim_y dim_G dim_g] = [1 1 1 6] 

if nargin<4 || isempty(keyxy)
    switch keyf
    case 'F'; w = -x-3*y;
    case 'G'; w = -x;    
    case 'f'; w = -x+3*y;  
    case 'g'; w = [-x-2*y+10; x-2*y-6; 2*x-y-21; x+2*y-38; -x+2*y-18; -y];
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
        case 'x' ; w = -1;    
        case 'y' ; w = 3;         
        case 'xx'; w = 0;
        case 'xy'; w = 0;
        case 'yy'; w = 0;
        end           
	case 'g'   
        switch keyxy
        case 'x' ; w = [-1; 1; 2; 1; -1; 0];
        case 'y' ; w = [-2; -2; -1; 2; 2; -1];             
        case 'xx'; w = zeros(6,1);
        case 'xy'; w = zeros(6,1);
        case 'yy'; w = zeros(6,1);
        end        
   end   
end

end




