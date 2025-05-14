function w=TollSettingP5(x,y,keyf,keyxy)
% This file provides all functions defining TollSettingP5 problem 
% and their first and second order derivatives.
% [dim_x dim_y dim_G dim_g] = [1 4 0 8]   
 
if nargin<4 || isempty(keyxy)
    switch keyf
    case 'F'; w  = -(y(2)+y(3))*x ;
    case 'G'; w  = [];
    case 'f'; w  = [8 3+2*x 4+2*x 6]*y;
    case 'g'; A  = [ 1 1 0 0; 0 0 1 1]; 
              w  = [[A;-A]*y+[-1 -1 1 1]'; -y];   
    end    
else 
    switch keyf
    case 'F'
        switch keyxy
        case 'x' ; w = -(y(2)+y(3));       
        case 'y' ; w = -[0; x; x; 0];
        case 'xx'; w = 0; 
        case 'xy'; w = -[0; 1; 1; 0];
        case 'yy'; w = sparse(zeros(4,4));
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
        case 'x' ; w = 2*(y(2)+y(3));   
        case 'y' ; w = [8; 3+2*x; 4+2*x; 6];        
        case 'xx'; w = 0; 
        case 'xy'; w = [0; 2; 2; 0];
        case 'yy'; w = sparse(zeros(4));    
        case 'yxx'; w = []; 
        case 'yxy'; w = [];
        case 'yyy'; w = [];
        end                
	case 'g'  
        switch keyxy
        case 'x' ; w = sparse(zeros(8,1)); 
        case 'y' ; A  = [ 1 1 0 0; 0 0 1 1]; 
                   w = sparse([A;-A; -eye(4)]);         
        case 'xx'; w = sparse(zeros(8,1));
        case 'xy'; w = sparse(zeros(32,1));
        case 'yy'; w = sparse(zeros(32,4));      
        case 'yxx'; w = []; 
        case 'yxy'; w = [];
        case 'yyy'; w = [];
        end
    end

end
end



