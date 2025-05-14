function w=ShehuEtal2019Ex42(x,y,keyf,keyxy,data)
% This file provides all functions defining IOC problem 
% and their first and second order derivatives.
% [dim_x dim_y dim_G dim_g] = [2 2*ni 3 4*ni]   

 
m  = length(y);

if nargin<4 || isempty(keyxy)
    switch keyf
    case 'F'; w  = y'*data.Q*y /2;
    case 'G'; w  = []; 
    case 'f'; w  = sum( (data.A*y-data.b).^2 ) /2 + data.mu* sum(abs(y)) ;
    case 'g'; w  = []; 
    end    
else
    switch keyf
    case 'F'
        switch keyxy
        case 'x' ; w  = 0;         
        case 'y' ; w  = data.Q*y;
        case 'xx'; w  = 0;
        case 'xy'; w  = zeros(m,1);  
        case 'yy'; w  = data.Q;
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
        case 'y' ; w = data.At*(data.A*y-data.b)+ data.mu*sign(y);      
        case 'xx'; w = 0;  
        case 'xy'; w = zeros(m,1);  
        case 'yy'; w = data.AA;
        case 'yxx'; w = [];
        case 'yxy'; w = [];
        case 'yyy'; w = [];  
        end           
	case 'g'   
        switch keyxy
        case 'x' ; w = [];    
        case 'y' ; w = [];         
        case 'xx'; w = []; 
        case 'xy'; w = []; 
        case 'yy'; w = [];
        case 'yxx'; w = [];
        case 'yxy'; w = [];
        case 'yyy'; w = [];              
        end       
   end   
end
end