function w=FrankeEtal2018Ex511(x,y,keyf,keyxy)
% This file provides all functions defining FrankeEtal18b problem 
% and their first and second order derivatives.
% [dim_x dim_y dim_G dim_g] = [1 3 0 4]   
 
if nargin<4 || isempty(keyxy)
    switch keyf
    case 'F'; w = 0.5*(y(1)-2)^2+0.5*y(2)^2+0.5*(y(3)-2)^2;   
    case 'G'; w = [];    
    case 'f'; w = y(1)+y(2)+y(3);  
    case 'g'; w = [-y(1)-y(2); -y(1)+y(2); -y(1); -y(3)];
    end    
else
    switch keyf
    case 'F'
        switch keyxy
        case 'x' ; w = 0; 
        case 'y' ; w = [y(1)-2; y(2); y(3)-2];
        case 'xx'; w = 0;
        case 'xy'; w = zeros(3,1);
        case 'yy'; w = [1 0 0; 0 1 0; 0 0 1];    
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
        case 'y' ; w = [1; 1; 1];
        case 'xx'; w = 0;
        case 'xy'; w = zeros(3,1);
        case 'yy'; w = zeros(3,3);    
        end           
	case 'g'   
        switch keyxy
        case 'x' ; w = zeros(4,1);             
        case 'y' ; w = [-1 -1 0; -1 1 0; -1 0 0; 0 0 -1];
        case 'xx'; w = zeros(4,1);
        case 'xy'; w = zeros(12,1);
        case 'yy'; w = zeros(12,3);     
        end        
   end   
end

end




