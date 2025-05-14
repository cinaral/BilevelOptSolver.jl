function w=DesignCentringP2(x,y,keyf,keyxy)
% This file provides all functions defining DesignCentringP2 problem 
% and their first and second order derivatives.
% [dim_x dim_y dim_G dim_g] = [4 6 5 3]   
 
if nargin<4 || isempty(keyxy)
    switch keyf
    case 'F'; w = -pi*x(3)*x(4);
    case 'G'; w = [-y(1)-y(2)^2; y(3)/4+y(4)-3/4; -y(6)-1; x(3)-1; x(4)-1];
    case 'f'; w = y(1)+y(2)^2-y(3)/4-y(4)+y(6);
    case 'g'; w = [ (x(1)-y(1))^2/x(3)^2+(x(2)-y(2))^2/x(4)^2;
                    (x(1)-y(3))^2/x(3)^2+(x(2)-y(4))^2/x(4)^2;
                    (x(1)-y(5))^2/x(3)^2+(x(2)-y(6))^2/x(4)^2];
    end    
else
    switch keyf
    case 'F'
        switch keyxy
        case 'x' ; w = [0; 0; -pi*x(4); -pi*x(3)];         
        case 'y' ; w = zeros(6,1); 
        case 'xx'; w = zeros(4,4); w(3,4)= -pi; w(4,3)= -pi; 
        case 'xy'; w = zeros(6,4);
        case 'yy'; w = zeros(6,6);
        end 
    case 'G'  
        switch keyxy
        case 'x' ; w = [zeros(3,4); 0 0 1 0; 0 0 0 1];   
        case 'y' ; w = [-1 -2*y(2) 0 0 0 0; 0 0 1/4 1 0 0; 0 0 0 0 0 -1; zeros(2,6)];         
        case 'xx'; w = zeros(20,4);
        case 'xy'; w = zeros(30,4);
        case 'yy'; w = zeros(30,6); w(2,2)=-2;
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
        xy1=x(1)-y(1); xy3=x(1)-y(3); xy5=x(1)-y(5);
        xy2=x(2)-y(2); xy4=x(2)-y(4); xy6=x(2)-y(6);
        switch keyxy
        case 'x' ;  w = 2*[xy1/x(3)^2, xy2/x(4)^2, -xy1^2/x(3)^3, -xy2^2/x(4)^3;
                           xy3/x(3)^2, xy4/x(4)^2, -xy3^2/x(3)^3, -xy4^2/x(4)^3;
                           xy5/x(3)^2, xy6/x(4)^2, -xy5^2/x(3)^3, -xy6^2/x(4)^3 ];    
        case 'y' ;  w = 2*[-xy1/x(3)^2, -xy2/x(4)^2, 0, 0, 0, 0;
                           0, 0, -xy3/x(3)^2, -xy4/x(4)^2, 0, 0;
                           0, 0, 0, 0, -xy5/x(3)^2, -xy6/x(4)^2];           
         case 'xx'; w1 =  [1/x(3)^2, 0, -2*xy1/x(3)^3, 0;
                           0 1/x(4)^2, 0, -2*xy2/x(4)^3;
                           -2*xy1/x(3)^3, 0, 3*xy1^2/x(3)^4 0;
                           0, -2*xy2/x(4)^3, 0,  3*xy2^2/x(4)^4];
                    w2 =  [1/x(3)^2, 0, -2*xy3/x(3)^3, 0;
                           0 1/x(4)^2, 0, -2*xy4/x(4)^3;
                           -2*xy3/x(3)^3, 0, 3*xy3^2/x(3)^4 0;
                           0, -2*xy4/x(4)^3, 0,  3*xy4^2/x(4)^4];
                    w3 =  [1/x(3)^2, 0, -2*xy5/x(3)^3, 0;
                           0 1/x(4)^2, 0, -2*xy6/x(4)^3;
                           -2*xy5/x(3)^3, 0, 3*xy5^2/x(3)^4 0;
                           0, -2*xy6/x(4)^3, 0,  3*xy6^2/x(4)^4];
                    w  = [w1;w2;w3]*2;  
         case 'xy'; w1 = [-1/x(3)^2 0 2*xy1/x(3)^3 0; 0 -1/x(4)^2 0 2*xy2/x(4)^3];
                    w2 = [-1/x(3)^2 0 2*xy3/x(3)^3 0; 0 -1/x(4)^2 0 2*xy4/x(4)^3];
                    w3 = [-1/x(3)^2 0 2*xy5/x(3)^3 0; 0 -1/x(4)^2 0 2*xy6/x(4)^3];
                    w  = [w1;zeros(6,4); w2 ;zeros(6,4);w3]*2; 
         case 'yy';  w = zeros(18,6);
                     w(1,1) =2/x(3)^2;  w(2,2) =2/x(4)^2; 
                     w(9,3) =2/x(3)^2;  w(10,4)=2/x(4)^2;
                     w(17,5)=2/x(3)^2;  w(18,6)=2/x(4)^2;   
        end
    end
end
end



