function w=DesignCentringP4(x,y,keyf,keyxy)
% This file provides all functions defining DesignCentringP4 problem 
% and their first and second order derivatives.
% [dim_x dim_y dim_G dim_g] = [4 6 3 3]   
 
if nargin<4 || isempty(keyxy)
    switch keyf
    case 'F'; w = -(x(1)-x(3))*(x(2)-x(4));
    case 'G'; w = [-y(1)-y(2)^2; y(3)/4+y(4)-3/4; -y(6)-1];
    case 'f'; w = y(1)+y(2)^2-y(3)/4-y(4)+y(6);
    case 'g'; A = [eye(2);eye(2);eye(2)]; 
              w = [y-A*x(1:2); -y+A*x(3:4)];
    end    
else
    switch keyf
    case 'F'
        switch keyxy
        case 'x' ; w = [x(4)-x(2); x(3)-x(1); x(2)-x(4); x(1)-x(3)];         
        case 'y' ; w = zeros(6,1); 
        case 'xx'; w = [0 -1 0 1; -1 0 1 0; 0 1 0 -1; 1 0 -1 0]; 
        case 'xy'; w = zeros(6,4);
        case 'yy'; w = zeros(6,6);
        end 
    case 'G'  
        switch keyxy
        case 'x' ; w = zeros(3,4);   
        case 'y' ; w = [-1 -2*y(2) 0 0 0 0; 0 0 1/4 1 0 0; 0 0 0 0 0 -1];         
        case 'xx'; w = zeros(12,4);
        case 'xy'; w = zeros(18,4);
        case 'yy'; w = zeros(18,6); w(2,2)=-2;
        end           
	case 'f'   
        switch keyxy
        case 'x' ;  w = zeros(4,1);    
        case 'y' ;  w = [1;2*y(2);-1/4;-1;0;1];           
        case 'xx';  w = zeros(4,4);  
        case 'xy';  w = zeros(6,4);  
        case 'yy';  w = zeros(6,6); w(2,2)=2;       
        case 'yxx'; w = [];  
        case 'yxy'; w = [];  
        case 'yyy'; w = [];   
        end                
	case 'g'  
        switch keyxy
        case 'x' ; A = [eye(2);eye(2);eye(2)]; 
                   w = [-A zeros(6,2); zeros(6,2) A]; 
        case 'y' ; w = [eye(6);-eye(6)];         
        case 'xx'; w = zeros(48,4);
        case 'xy'; w = zeros(72,4);
        case 'yy'; w = zeros(72,6);
        end
    end
end
end



