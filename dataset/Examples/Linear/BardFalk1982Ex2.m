function w=BardFalk1982Ex2(x,y,keyf,keyxy)
% This file provides all functions defining BardFalk1982Ex2 problem 
% and their first and second order derivatives.
% [dim_x dim_y dim_G dim_g] = [2 2 2 5]  

if nargin<4 || isempty(keyxy)
    switch keyf
    case 'F'; w = -2*x(1)+x(2)+0.5*y(1);
    case 'G'; w = -x;    
    case 'f'; w = x(1)+x(2)-4*y(1)+y(2);  
    case 'g'; w = [-2*x(1)+y(1)-y(2)+2.5; x(1)-3*x(2)+y(2)-2; x(1)+x(2)-2; -y];
    end    
else
    switch keyf
    case 'F'
        switch keyxy
        case 'x' ; w = [-2; 1];      
        case 'y' ; w = [0.5; 0];     
        case 'xx'; w = zeros(2,2);
        case 'xy'; w = zeros(2,2);
        case 'yy'; w = zeros(2,2);
        end 
    case 'G'  
        switch keyxy
        case 'x' ; w = [-1 0; 0 -1];    
        case 'y' ; w = [0 0; 0 0];         
        case 'xx'; w = zeros(4,2);
        case 'xy'; w = zeros(4,2);
        case 'yy'; w = zeros(4,2);
        end           
	case 'f'   
        switch keyxy
        case 'x' ; w = [1; 1];    
        case 'y' ; w = [-4; 1];         
        case 'xx'; w = zeros(2,2);
        case 'xy'; w = zeros(2,2);
        case 'yy'; w = zeros(2,2);
        end           
	case 'g'   
        switch keyxy
        case 'x' ; w = [-2 0; 1 -3; 1 1; 0 0; 0 0];
        case 'y' ; w = [1 -1; 0 1; 0 0; -1 0; 0 -1];             
        case 'xx'; w = zeros(10,2);
        case 'xy'; w = zeros(10,2);
        case 'yy'; w = zeros(10,2);
        end        
   end   
end

end




