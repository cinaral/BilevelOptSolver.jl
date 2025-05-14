function w=TollSettingP2(x,y,keyf,keyxy)
% This file provides all functions defining TollSettingP2 problem 
% and their first and second order derivatives.
% [dim_x dim_y dim_G dim_g] = [3 18 3 38]   
 
if nargin<4 || isempty(keyxy)
    switch keyf
    case 'F'; w  = -[y(1)+y(2) y(3)+y(4) y(5)+y(6)]*x ;
    case 'G'; w  = -x;
    case 'f'; w  = [2*x(1) 2*x(1) 2*x(2) 2*x(2) 2*x(3) 2*x(3) 5 7 14 7 2 4 29 20 12 8 5 2]*y;
    case 'g'; A  = sparse([zeros(1,6) 1 1 1 zeros(1,9); zeros(1,9) 1 1 1 zeros(1,6);
                           zeros(1,12) 1 1 1 zeros(1,3); zeros(1,15) 1 1 1;
                           1 0 0 0 1 0 -1 0 0 0 0 0 1 0 0 0 0 0;
                           0 1 0 0 0 1 0 0 0 -1 0 0 0 0 0 1 0 0;
                           -1 0 1 0 0 0 0 -1 0 0 0 0 0 1 0 0 0 0;
                           0 -1 0 1 0 0 0 0 0 0 -1 0 0 0 0 0 1 0;
                           0 0 -1 0 -1 0 0 0 -1 0 0 0 0 0 1 0 0 0;
                           0 0 0 -1 0 -1 0 0 0 0 0 -1 0 0 0 0 0 1]); 
              w  = [[A;-A]*y+[-ones(4,1);zeros(6,1);ones(4,1);zeros(6,1)]; -y];   
    end    
else 
    switch keyf
    case 'F'
        switch keyxy
        case 'x' ; w = -[y(1)+y(2); y(3)+y(4); y(5)+y(6)];       
        case 'y' ; w = -[x(1); x(1); x(2); x(2); x(3); x(3); zeros(12,1)];
        case 'xx'; w = sparse(zeros(3,3)); 
        case 'xy'; w = -sparse([1 0 0; 1 0 0; 0 1 0; 0 1 0; 0 0 1; 0 0 1; zeros(12,3)]);
        case 'yy'; w = sparse(zeros(18,18));
        end 
    case 'G'  
        switch keyxy
        case 'x' ; w = sparse(-eye(3));   
        case 'y' ; w = sparse(zeros(3,18));         
        case 'xx'; w = sparse(zeros(9,3));
        case 'xy'; w = sparse(zeros(54,3));
        case 'yy'; w = sparse(zeros(54,18)); 
        end           
	case 'f'   
        switch keyxy
        case 'x' ; w = 2*[y(1)+y(2); y(3)+y(4); y(5)+y(6)];   
        case 'y' ; w = [2*x(1);2*x(1);2*x(2);2*x(2);2*x(3);2*x(3);5;7;14;7;2;4;29;20;12;8;5;2];        
        case 'xx'; w = sparse(zeros(3,3)); 
        case 'xy'; w = sparse([2 0 0; 2 0 0; 0 2 0; 0 2 0; 0 0 2; 0 0 2; zeros(12,3)]);
        case 'yy'; w = sparse(zeros(18,18));      
        case 'yxx'; w = []; 
        case 'yxy'; w = [];
        case 'yyy'; w = [];
        end                
	case 'g'  
        switch keyxy
        case 'x' ; w = sparse(zeros(38,3)); 
        case 'y' ; A = [zeros(1,6) 1 1 1 zeros(1,9); zeros(1,9) 1 1 1 zeros(1,6);
                           zeros(1,12) 1 1 1 zeros(1,3); zeros(1,15) 1 1 1;
                           1 0 0 0 1 0 -1 0 0 0 0 0 1 0 0 0 0 0;
                           0 1 0 0 0 1 0 0 0 -1 0 0 0 0 0 1 0 0;
                           -1 0 1 0 0 0 0 -1 0 0 0 0 0 1 0 0 0 0;
                           0 -1 0 1 0 0 0 0 0 0 -1 0 0 0 0 0 1 0;
                           0 0 -1 0 -1 0 0 0 -1 0 0 0 0 0 1 0 0 0;
                           0 0 0 -1 0 -1 0 0 0 0 0 -1 0 0 0 0 0 1];
                   w = sparse([A;-A; -eye(18)]);         
        case 'xx'; w = sparse(zeros(114,3));
        case 'xy'; w = sparse(zeros(684,3));
        case 'yy'; w = sparse(zeros(684,18));      
        case 'yxx'; w = []; 
        case 'yxy'; w = [];
        case 'yyy'; w = [];
        end
    end

end
end



