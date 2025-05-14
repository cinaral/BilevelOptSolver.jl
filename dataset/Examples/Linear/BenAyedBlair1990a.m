function w=BenAyedBlair1990a(x,y,keyf,keyxy)
% This file provides all functions defining BenAyedBlair1990a problem 
% and their first and second order derivatives.
% [dim_x dim_y dim_G dim_g] = [1 2 2 4]  

if nargin<4 || isempty(keyxy)
    switch keyf
    case 'F'; w = -1.5*x-6*y(1)-y(2);
    case 'G'; w = [-x; x-1];    
    case 'f'; w = -y(1)-5*y(2);  
    case 'g'; w = [x+3*y(1)+y(2)-5; 2*x+y(1)+3*y(2)-5; -y];
    end    
else
    switch keyf
    case 'F'
        switch keyxy
        case 'x' ; w = 1.5;      
        case 'y' ; w = [-6; -1];     
        case 'xx'; w = 0;
        case 'xy'; w = zeros(2,1);
        case 'yy'; w = zeros(2,2);
        end 
    case 'G'  
        switch keyxy
        case 'x' ; w = [-1; 1];    
        case 'y' ; w = [0 0; 0 0];         
        case 'xx'; w = zeros(2,1);
        case 'xy'; w = zeros(4,1);
        case 'yy'; w = zeros(4,2);
        end           
	case 'f'   
        switch keyxy
        case 'x' ; w = 0;    
        case 'y' ; w = [-1; -5];         
        case 'xx'; w = 0;
        case 'xy'; w = zeros(2,1);
        case 'yy'; w = zeros(2,2);
        end           
	case 'g'   
        switch keyxy
        case 'x' ; w = [1; 2; 0; 0];
        case 'y' ; w = [3 1; 1 3; -1 0; 0 -1];             
        case 'xx'; w = zeros(4,1);
        case 'xy'; w = zeros(8,1);
        case 'yy'; w = zeros(8,2);
        end        
   end   
end

end




