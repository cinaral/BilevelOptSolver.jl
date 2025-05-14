function w=TollSettingP1(x,y,keyf,keyxy)
% This file provides all functions defining TollSettingP1 problem 
% and their first and second order derivatives.
% [dim_x dim_y dim_G dim_g] = [3 8 3 18]   
 
if nargin<4 || isempty(keyxy)
    switch keyf
    case 'F'; w  = -[y(3) y(4) y(8)]*x ;
    case 'G'; w  = -x;
    case 'f'; w  = [2 6 5+x(1) x(2) 4 2 6 x(3)]*y;
    case 'g'; A  = sparse([1 1 1 0 0 0 0 0; -1 0 0 1 1 0 0 0;
                   0 -1 0 -1 0 1 1 0; 0 0 0 0 -1 -1 0 1; 0 0 1 0 0 0 1 1]); 
              w  = [[A;-A]*y+[-1;0;0;0;-1;1;0;0;0;1]; -y];   
    end    
else
    switch keyf
    case 'F'
        switch keyxy
        case 'x' ; w = -[y(3);y(4);y(8)];       
        case 'y' ; w = -[0; 0; x(1); x(2); 0; 0; 0; x(3)];
        case 'xx'; w = zeros(3,3); 
        case 'xy'; w = -[zeros(2,3); 1 0 0; 0 1 0; zeros(3,3); 0 0 1];
        case 'yy'; w = zeros(8,8);
        end 
    case 'G'  
        switch keyxy
        case 'x' ; w = sparse(-eye(3));   
        case 'y' ; w = sparse(zeros(3,8));         
        case 'xx'; w = sparse(zeros(9,3));
        case 'xy'; w = sparse(zeros(24,3));
        case 'yy'; w = sparse(zeros(24,8)); 
        end           
	case 'f'   
        switch keyxy
        case 'x' ; w = [y(3);y(4);y(8)];    
        case 'y' ; w = [2; 6; 5+x(1); x(2); 4; 2; 6; x(3)];        
        case 'xx'; w = sparse(zeros(3,3)); 
        case 'xy'; w = sparse([zeros(2,3); 1 0 0; 0 1 0; zeros(3,3); 0 0 1]);
        case 'yy'; w = sparse(zeros(8,8));         
        case 'yxx'; w = []; 
        case 'yxy'; w = [];
        case 'yyy'; w = [];
        end                
	case 'g'  
        switch keyxy
        case 'x' ; w = sparse(zeros(18,3)); 
        case 'y' ; A = [1 1 1 0 0 0 0 0; -1 0 0 1 1 0 0 0;
                   0 -1 0 -1 0 1 1 0; 0 0 0 0 -1 -1 0 1; 0 0 1 0 0 0 1 1]; 
                   w = sparse([A;-A; -eye(8)]);         
        case 'xx'; w = sparse(zeros(54,3));
        case 'xy'; w = sparse(zeros(144,3));
        case 'yy'; w = sparse(zeros(144,8));      
        case 'yxx'; w = []; 
        case 'yxy'; w = [];
        case 'yyy'; w = [];
        end
    end
    w = sparse(w);
end
end



