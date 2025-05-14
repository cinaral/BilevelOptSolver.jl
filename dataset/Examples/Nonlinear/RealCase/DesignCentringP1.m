function w=DesignCentringP1(x,y,keyf,keyxy)
% This file provides all functions defining DesignCentringP1 problem 
% and their first and second order derivatives.
% [dim_x dim_y dim_G dim_g] = [3 6 3 3]   
 
if nargin<4 || isempty(keyxy)
    switch keyf
    case 'F'; w = -pi*x(3)^3;
    case 'G'; w = [-y(1)-y(2)^2; y(3)/4+y(4)-3/4;-y(6)-1];
    case 'f'; w = y(1)+y(2)^2-y(3)/4-y(4)+y(6);
    case 'g'; w = [ (y(1)-x(1))^2+(y(2)-x(2))^2-x(3)^2;
                    (y(3)-x(1))^2+(y(4)-x(2))^2-x(3)^2;
                    (y(5)-x(1))^2+(y(6)-x(2))^2-x(3)^2;];
    end    
else
    switch keyf
    case 'F'
        switch keyxy
        case 'x' ; w = [0; 0; -3*pi*x(3)^2];         
        case 'y' ; w = zeros(6,1); 
        case 'xx'; w = zeros(3,3); w(3,3)= -6*pi*x(3);  
        case 'xy'; w = zeros(6,3);
        case 'yy'; w = zeros(6,6);
        end 
    case 'G'  
        switch keyxy
        case 'x' ; w = zeros(3,3);   
        case 'y' ; w = [-1 -2*y(2) 0 0 0 0; 0 0 1/4 1 0 0; 0 0 0 0 0 -1];         
        case 'xx'; w = zeros(9,3);
        case 'xy'; w = zeros(18,3);
        case 'yy'; w = zeros(18,6); w(2,2)=-2;
        end           
	case 'f'   
        switch keyxy
        case 'x' ;  w = zeros(3,1);    
        case 'y' ;  w = [1;2*y(2);-1/4;-1;0;1];           
        case 'xx';  w = zeros(3,3);  
        case 'xy';  w = zeros(6,3);  
        case 'yy';  w = zeros(6,6);  w(2,2) =2;         
        case 'yxx';  w = [];  
        case 'yxy';  w = [];  
        case 'yyy';  w = [];   
        end                
	case 'g'   
        switch keyxy
        case 'x' ;  w = 2*[x(1)-y(1) x(2)-y(2) -x(3);
                           x(1)-y(3) x(2)-y(4) -x(3);
                           x(1)-y(5) x(2)-y(6) -x(3)];    
        case 'y' ;  w = 2*[y(1)-x(1) y(2)-x(2) 0 0 0 0;
                           0 0 y(3)-x(1) y(4)-x(2) 0 0;
                           0 0 0 0 y(5)-x(1) y(6)-x(2)];           
         case 'xx';  w = [diag([2 2 -2]);diag([2 2 -2]);diag([2 2 -2])];  
         case 'xy';  w = [diag([-2 -2 0]);zeros(5,3);
                          diag([-2 -2 0]);zeros(5,3);-2 0 0; 0 -2 0];    
         case 'yy';  w = zeros(18,6);
                     w(1,1) =2;  w(2,2) =2; 
                     w(9,3) =2;  w(10,4)=2;
                     w(17,5)=2;  w(18,6)=2;           
        case 'yxx';  w = [];  
        case 'yxy';  w = [];  
        case 'yyy';  w = [];  
        end
    end
end
end



