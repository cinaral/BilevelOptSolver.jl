function w=RobustPortfolioP2(x,y,keyf,keyxy)
% This file provides all functions defining RobustPortfolioP2 problem 
% and their first and second order derivatives.
% [dim_x dim_y dim_G dim_g] = [N+1 N N+3 N+1]   

N     = length(y);  
 
if nargin<4 || isempty(keyxy)
    switch keyf
    case 'F'; w  = -x(end);  
    case 'G'; w  = [x(end)-y'*x(1:N); -x(1:N); sum(x(1:N))-1; 1-sum(x(1:N))];
    case 'f'; w  = y'*x(1:N)-x(end);
    case 'g'; I  = (1:N)'; 
              si = ((0.05/3/N)^2*(2*N*(N+1)*I));
              yi = 1.15+(0.05/N)*I;
              sx = 1.5*(1+sum((x(1:N)-1/N).^2));
              w  = [sum((abs(y-yi)).^2./si)-sx^2;-y];   
    end    
else
    switch keyf
    case 'F'
        switch keyxy
        case 'x' ; w = sparse([zeros(N,1);-1]);         
        case 'y' ; w = sparse( zeros(N,1));
        case 'xx'; w = sparse( zeros(N+1,N+1)); 
        case 'xy'; w = sparse( zeros(N,N+1)); 
        case 'yy'; w = sparse( zeros(N,N)); 
        end 
    case 'G'  
        switch keyxy
        case 'x' ; w = sparse([-y' 1; -eye(N) zeros(N,1); ones(1,N) 0; -ones(1,N) 0]); 
        case 'y' ; w = sparse([-x(1:N)'; zeros(N+2,N)]);       
        case 'xx'; w = sparse(zeros((N+3)*(N+1),N+1));
        case 'xy'; w = sparse([-eye(N) zeros(N,1); zeros((N+2)*N,N+1)]);
        case 'yy'; w = sparse(zeros((N+3)*N,N));
        end           
	case 'f'   
        switch keyxy
        case 'x' ; w = [y;-1];         
        case 'y' ; w =  x(1:N);
        case 'xx'; w = sparse(zeros(N+1,N+1)); 
        case 'xy'; w = sparse([eye(N) zeros(N,1)]); 
        case 'yy'; w = sparse(zeros(N,N));   
        case 'yxx'; w = []; 
        case 'yxy'; w = []; 
        case 'yyy'; w = [];
        end                
	case 'g'  
        switch keyxy
        case 'x' ; sx = 1+sum((x(1:N)-1/N).^2);
                   w  = sparse([-9*sx*(x(1:N)-1/N)' 0; -zeros(N,N+1)]); 
        case 'y' ; I  = (1:N)'; 
                   si = ((0.05/3/N)*sqrt(2*N*(N+1)*I)).^2;
                   yi = 1.15+(0.05/N)*I;
                   w  = sparse([( 2*(y-yi)./si )'; -eye(N,N)]);         
        case 'xx'; sx = 1+sum((x(1:N)-1/N).^2);
                   A  = -9*(sx*diag(ones(N,1))+2*(x(1:N)-1/N)*(x(1:N)-1/N)');
                   w  = sparse([A zeros(N,1); zeros(1,N+1); zeros(N*(N+1),N+1)]);
        case 'xy'; w  = sparse(zeros((N+1)*N,N+1));
        case 'yy'; I  = (1:N)'; 
                   si = ((0.05/3/N)^2*(2*N*(N+1)*I));
                   w  = sparse([diag(2./si); zeros(N*N,N)]);
        end
    end
 
end

end



