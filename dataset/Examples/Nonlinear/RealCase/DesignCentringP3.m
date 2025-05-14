function w=DesignCentringP3(x,y,keyf,keyxy)
% This file provides all functions defining DesignCentringP3 problem 
% and their first and second order derivatives.
% [dim_x dim_y dim_G dim_g] = [6 6 3 3]   
 
if nargin<4 || isempty(keyxy)
    switch keyf
    case 'F'; w = -pi*abs(x(3)*x(6)-x(4)*x(5));
    case 'G'; w = [-y(1)-y(2)^2; y(3)/4+y(4)-3/4;-y(6)-1];
    case 'f'; w = y(1)+y(2)^2-y(3)/4-y(4)+y(6);
    case 'g'; z = [x(1);x(2)]; 
              A = inv([x(3) x(4); x(5) x(6)]*[x(3) x(5); x(4) x(6)]);%
              w = [ ([y(1);y(2)]-z)'*A*([y(1);y(2)]-z)-1;
                    ([y(3);y(4)]-z)'*A*([y(3);y(4)]-z)-1;
                    ([y(5);y(6)]-z)'*A*([y(5);y(6)]-z)-1];
    end    
else
    switch keyf
    case 'F'
        t  = x(3)*x(6)-x(4)*x(5);
        st = sign(t);
        switch keyxy
        case 'x' ; w = pi*st*[0 0 -x(6) x(5) x(4) -x(3)]';         
        case 'y' ; w = zeros(6,1);  
        case 'xy'; w = zeros(6,6);
        case 'yy'; w = zeros(6,6);
        end 
    case 'G'  
        switch keyxy
        case 'x' ; w = zeros(3,6);   
        case 'y' ; w = [-1 -2*y(2) 0 0 0 0; 0 0 1/4 1 0 0; 0 0 0 0 0 -1];         
        case 'xx'; w = zeros(18,6);
        case 'xy'; w = zeros(18,6);
        case 'yy'; w = zeros(18,6); w(2,2)=-2;
        end           
	case 'f'   
        switch keyxy
        case 'x' ;  w = zeros(6,1);    
        case 'y' ;  w = [1;2*y(2);-1/4;-1;0;1];           
        case 'xx';  w = zeros(6,6);  
        case 'xy';  w = zeros(6,6);  
        case 'yy';  w = zeros(6,6);  w(2,2)=2;
        end                
   end   
end
end



