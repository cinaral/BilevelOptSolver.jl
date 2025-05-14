function w=HuHuangZhang2009(x,y,keyf,keyxy)
% This file provides all functions defining HuHuangZhang2009 problem 
% and their first and second order derivatives.
% [dim_x dim_y dim_G dim_g] = [1 2 1 5]  

if nargin<4 || isempty(keyxy)
    switch keyf
    case 'F'; w = -4*x-y(1)-y(2);
    case 'G'; w = -x;    
    case 'f'; w = -x-3*y(1);  
    case 'g'; w = [x+y(1)+y(2)-25/9; x+y(2)-2; y(1)+y(2)-8/9; -y];
    end    
else
    switch keyf
    case 'F'
        switch keyxy
        case 'x' ; w = -4;      
        case 'y' ; w = [-1; -1];     
        case 'xx'; w = 0;
        case 'xy'; w = zeros(2,1);
        case 'yy'; w = zeros(2,2);
        end 
    case 'G'  
        switch keyxy
        case 'x' ; w = -1;    
        case 'y' ; w = zeros(1,2);         
        case 'xx'; w = 0;
        case 'xy'; w = zeros(2,1);
        case 'yy'; w = zeros(2,2);
        end           
	case 'f'   
        switch keyxy
        case 'x' ; w = -1;    
        case 'y' ; w = [-3; 0];         
        case 'xx'; w = 0;
        case 'xy'; w = zeros(2,1);
        case 'yy'; w = zeros(2,2);
        end           
	case 'g'   
        switch keyxy
        case 'x' ; w = [1; 1; 0; 0; 0];
        case 'y' ; w = [1 1; 0 1; 1 1; -1 0; 0 -1];          
        case 'xx'; w = zeros(5,1);
        case 'xy'; w = zeros(10,1);
        case 'yy'; w = zeros(10,2);
        end        
   end   
end

end




