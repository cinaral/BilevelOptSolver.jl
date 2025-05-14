function w=Bard1984b(x,y,keyf,keyxy)
% This file provides all functions defining Bard1984b problem 
% and their first and second order derivatives.
% [dim_x dim_y dim_G dim_g] = [1 1 1 5] 

if nargin<4 || isempty(keyxy)
    switch keyf
    case 'F'; w = -5*x-y;
    case 'G'; w = -x;    
    case 'f'; w = y;  
    case 'g'; w = [-x-0.5*y+2; -0.25*x+y-2; x+0.5*y-8; x-2*y-4; -y];
    end    
else
    switch keyf
    case 'F'
        switch keyxy
        case 'x' ; w = -5;      
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
        case 'x' ; w = [-1; -0.25; 1; 1; 0];
        case 'y' ; w = [-0.5; 1; 0.5; -2; -1];             
        case 'xx'; w = zeros(5,1);
        case 'xy'; w = zeros(5,1);
        case 'yy'; w = zeros(5,1);
        end        
   end   
end

end




