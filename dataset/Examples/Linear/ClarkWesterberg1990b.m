function w=ClarkWesterberg1990b(x,y,keyf,keyxy)
% This file provides all functions defining ClarkWesterberg90b problem 
% and their first and second order derivatives.
% [dim_x dim_y dim_G dim_g] = [1 2 2 5]  

if nargin<4 || isempty(keyxy)
    switch keyf
    case 'F'; w = -x-3*y(1)+2*y(2);
    case 'G'; w = [-x; x-8];    
    case 'f'; w = -y(1);  
    case 'g'; w = [0; 0; -2; 8; -2]*x+...
                  [-1 0; 1 0; 1 4; 3 -2; 1 -3]*y+[0; -4; -16; -48; 12];
    end    
else
    switch keyf
    case 'F'
        switch keyxy
        case 'x' ; w = -1;      
        case 'y' ; w = [-3; 2];     
        case 'xx'; w = 0;
        case 'xy'; w = zeros(2,1);
        case 'yy'; w = zeros(2,2);
        end 
    case 'G'  
        switch keyxy
        case 'x' ; w = [-1; 1];    
        case 'y' ; w = zeros(2,2);         
        case 'xx'; w = zeros(2,1);
        case 'xy'; w = zeros(4,1);
        case 'yy'; w = zeros(4,2);
        end           
	case 'f'   
        switch keyxy
        case 'x' ; w = 0;    
        case 'y' ; w = [-1; 0];         
        case 'xx'; w = 0;
        case 'xy'; w = zeros(2,1);
        case 'yy'; w = zeros(2,2);
        end           
	case 'g'   
        switch keyxy
        case 'x' ; w = [0; 0; -2; 8; -2];
        case 'y' ; w = [-1 0; 1 0; 1 4; 3 -2; 1 -3];          
        case 'xx'; w = zeros(5,1);
        case 'xy'; w = zeros(10,1);
        case 'yy'; w = zeros(10,2);
        end        
   end   
end

end




