function w=MershaDempe2006Ex2(x,y,keyf,keyxy)
% This file provides all functions defining MershaDempe2006Ex2 problem 
% and their first and second order derivatives.
% [dim_x dim_y dim_G dim_g] = [1 1 2 2]  

if nargin<4 || isempty(keyxy)
    switch keyf
    case 'F'; w = -x-2*y;
    case 'G'; w = [-2*x+3*y-12; x+y-14];    
    case 'f'; w = -y;  
    case 'g'; w = [-3*x+y+3; 3*x+y-30];
    end    
else
    switch keyf
    case 'F'
        switch keyxy
        case 'x' ; w = -1;      
        case 'y' ; w = -2;     
        case 'xx'; w = 0;
        case 'xy'; w = 0;
        case 'yy'; w = 0;
        end 
    case 'G'  
        switch keyxy
        case 'x' ; w = [-2; 1];    
        case 'y' ; w = [3; 1];         
        case 'xx'; w = zeros(2,1);
        case 'xy'; w = zeros(2,1);
        case 'yy'; w = zeros(2,1);
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
        case 'x' ; w = [-3; 3];
        case 'y' ; w = [1; 1];             
        case 'xx'; w = zeros(2,1);
        case 'xy'; w = zeros(2,1);
        case 'yy'; w = zeros(2,1);
        end        
   end   
end

end




