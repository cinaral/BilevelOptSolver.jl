function w=OptimalControl(x,y,keyf,keyxy,data)
% This file provides all functions defining OptimalControl problem 
% and their first and second order derivatives.
% [dim_x dim_y dim_G dim_g] = [2 2*ni 3 4*ni]   

 
ip   = data.ip;
ni   = data.ni;
n    = data.n ;
ni2  = 2*ni;
ni4  = 4*ni;

if keyf=='f'
d   = data.C*y(1:ni) - data.P*x;
e   = y(ni+1:2*ni)   - data.Q*x;
end

if nargin<4 || isempty(keyxy)
    switch keyf
    case 'F'; yy = zeros(n,1); yy(ip) = y(1:ni);
              w  = 1/2 * (yy - data.target)' * data.M_full * (yy - data.target) - x'*data.k;
    case 'G'; w  = [x(1)+x(2)-1; -x]; 
    case 'f'; w  = (d'*data.W*d + data.sigma*(e'*data.U*e)) / 2;
    case 'g'; eq = [data.A -data.B]*y;
              w  = [y(ni+1:2*ni)-data.ub; data.ua-y(ni+1:2*ni); eq; -eq]; 
    end    
else
    switch keyf
    case 'F'
        switch keyxy
        case 'x' ; w  = -data.k;         
        case 'y' ; yy = zeros(n,1); yy(ip) = y(1:ni);
                   w  = [data.M_full(ip,:) * (yy - data.target); zeros(ni,1)];
        case 'xx'; w  = zeros(2);
        case 'xy'; w  = sparse(ni2,2);  
        case 'yy'; w  = sparse(ni2,ni2); w(1:ni,1:ni)=data.M_full(ip,ip);
        end 
    case 'G'  
        switch keyxy
        case 'x' ; w = [1  1; -1 0; 0 -1];    
        case 'y' ; w = sparse(3,ni2);         
        case 'xx'; w = []; 
        case 'xy'; w = []; 
        case 'yy'; w = []; 
        end           
	case 'f'   
        switch keyxy
        case 'x' ; w = -data.Pt*data.W*d - data.sigma*data.Qt*data.U*e ;  
        case 'y' ; w = [data.Ct*data.W*d; data.sigma*data.U*e];      
        case 'xx'; w = data.fxx;  
        case 'xy'; w = data.fxy;
        case 'yy'; w = data.fyy;
        case 'yxx'; w = [];
        case 'yxy'; w = [];
        case 'yyy'; w = [];  
        end           
	case 'g'   
        switch keyxy
        case 'x' ; w = sparse(ni4,2);
        case 'y' ; w = sparse([zeros(ni) eye(ni); zeros(ni) -eye(ni);...
                               data.A -data.B; -data.A data.B]);              
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