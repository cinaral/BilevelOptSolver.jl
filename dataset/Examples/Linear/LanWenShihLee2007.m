function w=LanWenShihLee2007(x,y,keyf,keyxy)
% This file provides all functions defining LanWenShihLee2007 problem 
% and their first and second order derivatives.
% [dim_x dim_y dim_G dim_g] = [1 1 1 7]  

if nargin<4 || isempty(keyxy)
    switch keyf
    case 'F'; w = 2*x-11*y;   
    case 'G'; w = -x;
    case 'f'; w = x+3*y;  
    case 'g'; w = [1; 2; 3; 1; -4; -1; 0]*x+[-2; -1; 4; 7; 5; -4; -1]*y+...
                  [-4; -24; -96; -126; -65; 8;  0];
    end    
else
    switch keyf
    case 'F'
        switch keyxy
        case 'x' ; w = 2;      
        case 'y' ; w = -11;     
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
        case 'x' ; w = 1;    
        case 'y' ; w = 3;         
        case 'xx'; w = 0;
        case 'xy'; w = 0;
        case 'yy'; w = 0;
        end           
	case 'g'   
        switch keyxy
        case 'x' ; w = [1; 2; 3; 1; -4; -1;  0];
        case 'y' ; w = [-2; -1; 4; 7; 5; -4; -1];             
        case 'xx'; w = zeros(7,1);
        case 'xy'; w = zeros(7,1);
        case 'yy'; w = zeros(7,1);
        end        
   end   
end

end




