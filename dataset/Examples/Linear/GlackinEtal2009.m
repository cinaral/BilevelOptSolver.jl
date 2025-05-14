function w=GlackinEtal2009(x,y,keyf,keyxy)
% This file provides all functions defining GlackinEtal09 problem 
% and their first and second order derivatives.
% [dim_x dim_y dim_G dim_g] = [2 1 3 3]  

if nargin<4 || isempty(keyxy)
    switch keyf
    case 'F'; w = -2*x(1)+4*x(2)+3*y;
    case 'G'; w = [x(1)-x(2)+1; -x];    
    case 'f'; w = -y;  
    case 'g'; w = [x(1)+x(2)+y-4; 2*x(1)+2*x(2)+y-6; -y];
    end    
else
    switch keyf
    case 'F'
        switch keyxy
        case 'x' ; w = [-2; 4];      
        case 'y' ; w = 3;     
        case 'xx'; w = zeros(2,2);
        case 'xy'; w = zeros(1,2);
        case 'yy'; w = 0;
        end 
    case 'G'  
        switch keyxy
        case 'x' ; w = [1 -1; -1 0; 0 -1];    
        case 'y' ; w = zeros(3,1);         
        case 'xx'; w = zeros(6,2);
        case 'xy'; w = zeros(3,2);
        case 'yy'; w = zeros(3,1);
        end           
	case 'f'   
        switch keyxy
        case 'x' ; w = [0; 0];    
        case 'y' ; w = -1;         
        case 'xx'; w = zeros(2,2);
        case 'xy'; w = zeros(1,2);
        case 'yy'; w = 0;
        end           
	case 'g'   
        switch keyxy
        case 'x' ; w = [1 1; 2 2; 0 0];
        case 'y' ; w = [1; 1; -1];          
        case 'xx'; w = zeros(6,2);
        case 'xy'; w = zeros(3,2);
        case 'yy'; w = zeros(3,1);
        end        
   end   
end

end




