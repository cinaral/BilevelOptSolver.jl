function w=FrankeEtal2018Ex521(x,y,keyf,keyxy)
% This file provides all functions defining FrankeEtal18d problem 
% and their first and second order derivatives.
% [dim_x dim_y dim_G dim_g] = [1 2 0 3]   
 
if nargin<4 || isempty(keyxy)
    switch keyf
    case 'F'; w = -y(2);    
    case 'G'; w = [];     
    case 'f'; w = y(1);  
    case 'g'; w = [(y(1)-1)^2-(y(2)-0.5)^2-1.25; y(1)+y(2)-1; -y(1)];
    end    
else
    switch keyf
    case 'F'
        switch keyxy
        case 'x' ; w = 0;    
        case 'y' ; w = [0; -1];          
        case 'xx'; w = 0;
        case 'xy'; w = zeros(2,1);
        case 'yy'; w = zeros(2,2);  
        end        
    case 'G'
        switch keyxy
        case 'x' ; w = []; 
        case 'y' ; w = [];
        case 'xx'; w = [];
        case 'xy'; w = [];
        case 'yy'; w = [];    
        end        
	case 'f'   
        switch keyxy
        case 'x' ; w = 0;    
        case 'y' ; w = [1; 0];          
        case 'xx'; w = 0;
        case 'xy'; w = zeros(2,1);
        case 'yy'; w = zeros(2,2); 
        end           
	case 'g'   
        switch keyxy
        case 'x' ; w = [0;0;0];    
        case 'y' ; w = [2*y(1)-2 -2*y(2)+1; 1 1; -1 0];         
        case 'xx'; w = zeros(3,1);
        case 'xy'; w = zeros(6,1);           
        case 'yy'; w = [2 0; -2 0; 0 0; 0 0; 0 0; 0 0];
        end        
   end   
end

end




