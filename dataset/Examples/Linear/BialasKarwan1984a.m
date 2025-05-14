function w=BialasKarwan1984a(x,y,keyf,keyxy)
% This file provides all functions defining BialasKarwan84a problem 
% and their first and second order derivatives.
% [dim_x dim_y dim_G dim_g] = [1 2 1 7]  

if nargin<4 || isempty(keyxy)
    switch keyf
    case 'F'; w = -x-y(2);
    case 'G'; w = -x;    
    case 'f'; w = -y(2);  
    case 'g'; w = [1; -1; -1; 1; 0; 0; 0]*x+[1 1; -1 1; 1 1; -1 1; 0 1; -eye(2)]*y+...
                  [-3; 1; -1; -1; -0.5; 0; 0];
    end    
else
    switch keyf
    case 'F'
        switch keyxy
        case 'x' ; w = -1;      
        case 'y' ; w = [0; -1];     
        case 'xx'; w = 0;
        case 'xy'; w = zeros(2,1);
        case 'yy'; w = zeros(2,2);
        end 
    case 'G'  
        switch keyxy
        case 'x' ; w = -1;    
        case 'y' ; w = [0 0];         
        case 'xx'; w = 0;
        case 'xy'; w = zeros(2,1);
        case 'yy'; w = zeros(2,2);
        end           
	case 'f'   
        switch keyxy
        case 'x' ; w = 0;    
        case 'y' ; w = [0; -1];         
        case 'xx'; w = 0;
        case 'xy'; w = zeros(2,1);
        case 'yy'; w = zeros(2,2);
        end           
	case 'g'   
        switch keyxy
        case 'x' ; w = [1; -1; -1; 1; 0; 0; 0];
        case 'y' ; w = [1 1; -1 1; 1 1; -1 1; 0 1; -eye(2)];          
        case 'xx'; w = zeros(7,1);
        case 'xy'; w = zeros(14,1);
        case 'yy'; w = zeros(14,2);
        end        
   end   
end

end




