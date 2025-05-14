function w=CandlerTownsley1982(x,y,keyf,keyxy)
% This file provides all functions defining CandlerTownsley1982 problem 
% and their first and second order derivatives.
% [dim_x dim_y dim_G dim_g] = [2 3 2 6]  

if nargin<4 || isempty(keyxy)
    switch keyf
    case 'F'; w = -[8 4]*x + [4 -40 -4]*y;
    case 'G'; w = -x;    
    case 'f'; w = [1 2]*x + [1 1 2]*y;  
    case 'g'; w = [0 0; 2 0; 0 2; zeros(3,2)]*x+...
                  [-1 1 1; -1 2 -0.5; 2 -1 -0.5; -eye(3)]*y-...
                  [1; 1; 1; 0; 0; 0;];
    end    
else
    switch keyf
    case 'F'
        switch keyxy
        case 'x' ; w = [-8; -4];      
        case 'y' ; w = [4; -40; -4];     
        case 'xx'; w = zeros(2,2);
        case 'xy'; w = zeros(3,2);
        case 'yy'; w = zeros(3,3);
        end 
    case 'G'  
        switch keyxy
        case 'x' ; w = [-1 0; 0 -1];    
        case 'y' ; w = zeros(2,3);         
        case 'xx'; w = zeros(4,2);
        case 'xy'; w = zeros(6,2);
        case 'yy'; w = zeros(6,3);
        end           
	case 'f'   
        switch keyxy
        case 'x' ; w = [1; 2];    
        case 'y' ; w = [1; 1; 2];         
        case 'xx'; w = zeros(2,2);
        case 'xy'; w = zeros(3,2);
        case 'yy'; w = zeros(3,3);
        end           
	case 'g'   
        switch keyxy
        case 'x' ; w = [0 0; 2 0; 0 2; zeros(3,2)];
        case 'y' ; w = [-1 1 1; -1 2 -0.5; 2 -1 -0.5; -eye(3)];          
        case 'xx'; w = zeros(12,2);
        case 'xy'; w = zeros(18,2);
        case 'yy'; w = zeros(18,3);
        end        
   end   
end

end




