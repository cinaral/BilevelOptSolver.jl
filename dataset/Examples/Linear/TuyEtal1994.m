function w=TuyEtal1994(x,y,keyf,keyxy)
% This file provides all functions defining TuyEtal94 problem 
% and their first and second order derivatives.
% [dim_x dim_y dim_G dim_g] = [2 2 3 3]  

if nargin<4 || isempty(keyxy)
    switch keyf
    case 'F'; w = 3*x(1)+2*x(2)+y(1)+y(2);
    case 'G'; w = [x(1)+x(2)+y(1)+y(2)-4; -x];    
    case 'f'; w = 4*y(1)+y(2);  
    case 'g'; w = [-3*x(1)-5*x(2)-6*y(1)-2*y(2)+15; -y];
    end    
else
    switch keyf
    case 'F'
        switch keyxy
        case 'x' ; w = [3; 2];      
        case 'y' ; w = [1; 1];     
        case 'xx'; w = zeros(2,2);
        case 'xy'; w = zeros(2,2);
        case 'yy'; w = zeros(2,2);
        end 
    case 'G'  
        switch keyxy
        case 'x' ; w = [1 1; -eye(2)];    
        case 'y' ; w = zeros(3,2);         
        case 'xx'; w = zeros(6,2);
        case 'xy'; w = zeros(6,2);
        case 'yy'; w = zeros(6,2);
        end           
	case 'f'   
        switch keyxy
        case 'x' ; w = [0; 0];    
        case 'y' ; w = [4; 1];         
        case 'xx'; w = zeros(2,2);
        case 'xy'; w = zeros(2,2);
        case 'yy'; w = zeros(2,2);
        end           
	case 'g'   
        switch keyxy
        case 'x' ; w = [-3 -5; 0  0; 0  0];
        case 'y' ; w = [-6 -2; -1 0; 0 -1];          
        case 'xx'; w = zeros(6,2);
        case 'xy'; w = zeros(6,2);
        case 'yy'; w = zeros(6,2);
        end        
   end   
end

end




