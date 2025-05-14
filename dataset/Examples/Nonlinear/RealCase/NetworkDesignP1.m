function w=NetworkDesignP1(x,y,keyf,keyxy)
% This file provides all functions defining DesignCentringP4 problem 
% and their first and second order derivatives.
% [dim_x dim_y dim_G dim_g] = [5 5 5 11]   
 
if nargin<4 || isempty(keyxy)
    switch keyf
    case 'F'; xy = y./(1+x); 
              w  = [50+xy(1) 10*xy(2) 10+xy(3) 10*xy(4) 50+xy(5)]*y+100*sum(x) ;
    case 'G'; w  = -x-1;
    case 'f'; xy = y./(1+x)/2; 
              w  = [50+xy(1) 10*xy(2) 10+xy(3) 10*xy(4) 50+xy(5)]*y;
    case 'g'; A = [1 0 1 0 1; 0 1 -1 0 -1; -1 0 -1 1 0]; 
              w = [[A;-A]*y+[-6;0;0;6;0;0];-y];   
    end    
else
    switch keyf
    case 'F'
        a=[2;20;2;20;2];
        switch keyxy
        case 'x' ; w = -(a/2).*(y.*y./(1+x).^2)+100*ones(5,1);         
        case 'y' ; w =  a.*(y./(1+x))+[50;0;10;0;50];
        case 'xx'; w =  diag(a.*(y.*y./(1+x).^3)); 
        case 'xy'; w = -diag(a.*(y./(1+x).^2)); 
        case 'yy'; w =  diag(a./(1+x)); 
        end 
    case 'G'  
        switch keyxy
        case 'x' ; w = -eye(5,5);   
        case 'y' ; w = zeros(5,5);         
        case 'xx'; w = zeros(25,5);
        case 'xy'; w = zeros(25,5);
        case 'yy'; w = zeros(25,5); 
        end           
	case 'f'   
        a=[1;10;1;10;1];
        switch keyxy
        case 'x' ;  w = -(a/2).*(y.*y./(1+x).^2);    
        case 'y' ;  w =  a.*(y./(1+x))+[50;0;10;0;50];        
        case 'xx';  w =  diag(a.*(y.*y./(1+x).^3));
        case 'xy';  w = -diag(a.*(y./(1+x).^2));   
        case 'yy';  w =  diag(a./(1+x));     
        end                
	case 'g'  
        switch keyxy
        case 'x' ; w = zeros(11,5); 
        case 'y' ; A = [1 0 1 0 1; 0 1 -1 0 -1; -1 0 -1 1 0]; 
                   w = [A;-A; -eye(5)];         
        case 'xx'; w = zeros(55,5);
        case 'xy'; w = zeros(55,5);
        case 'yy'; w = zeros(55,5);
        end
    end
end
end



