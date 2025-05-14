function w=MorganPatrone2006b(x,y,keyf,keyxy)
% This file provides all functions defining MorganPatrone2006b problem 
% and their first and second order derivatives.
% [dim_x dim_y dim_G dim_g] = [1 1 0 4]

  
if nargin<4 || isempty(keyxy)
    switch keyf
    case 'F'; w = -(x+y);
    case 'G'; w = [];
    case 'f'; w = fxy(x,y); 
    case 'g'; w = [-x-0.5;x-0.5;-y-1; y-1];  
    end    
else
    switch keyf
    case 'F'
        switch keyxy
        case 'x' ; w = -1;         
        case 'y' ; w = -1;    
        case 'xx'; w =  0;
        case 'xy'; w =  0;
        case 'yy'; w =  0;
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
        case 'x' ; w = fxy(x,y,'x');    
        case 'y' ; w = fxy(x,y,'y');         
        case 'xx'; w = fxy(x,y,'xx');
        case 'xy'; w = fxy(x,y,'xy');
        case 'yy'; w = fxy(x,y,'yy');
        case 'yxx';  w = [];  
        case 'yxy';  w = []; 
        case 'yyy';  w = []; 
        end           
	case 'g'   
        switch keyxy
        case 'x' ; w = [-1; 1; 0; 0];
        case 'y' ; w = [ 0; 0; -1; 1];             
        case 'xx'; w = [ 0; 0; 0;0];
        case 'xy'; w = [ 0; 0; 0;0];
        case 'yy'; w = [ 0; 0; 0;0];
        case 'yxx';  w = [];  
        case 'yxy';  w = []; 
        case 'yyy';  w = [];             
        end        
   end   
end
end


function w = fxy(x,y,keyxy)
if nargin<3
    if     x>-0.51 && x<=-0.25;      w=(x+0.25)*y;
    elseif x>-0.25 && x<0.25;         w=0;
    elseif x>=0.25 && x<0.51;        w=(x-0.25)*y;
    else;                             w=10^10;
    end
else
    if     x>-0.51 && x<=-0.25  
           switch keyxy 
           case 'x' ; w = y;
           case 'y' ; w = x+0.25;             
           case 'xx'; w = 0;
           case 'xy'; w = 1;
           case 'yy'; w = 0;
           end  
    elseif x>=0.25 && x < 0.51  
           switch keyxy
           case 'x' ; w = y;
           case 'y' ; w = x-0.25;             
           case 'xx'; w = 0;
           case 'xy'; w = 1;
           case 'yy'; w = 0;
           end          
    else                      
           switch keyxy
           case 'x' ; w = 0;
           case 'y' ; w = 0;             
           case 'xx'; w = 0;
           case 'xy'; w = 0;
           case 'yy'; w = 0;
           end          
    end
    
end    
end

% function w=MorganPatrone2006b(x,y,keyf,keyxy)
% % This file provides all functions defining MorganPatrone2006b problem 
% % and their first and second order derivatives.
% % [dim_x dim_y dim_G dim_g] = [1 1 0 4]
% 
%   
% if nargin<4 || isempty(keyxy)
%     switch keyf
%     case 'F'; w = -(x+y);
%     case 'G'; w = [];
%     case 'f'; w = (x+0.25)*y*(-0.5<=x && x<=-0.25)+(x-0.25)*y*(0.25<=x && x<=0.5)+10^10*(-0.5>x || x>0.5); 
%     case 'g'; w = [-x-0.5;x-0.5;-y-1; y-1];  
%     end    
% else
%     switch keyf
%     case 'F'
%         switch keyxy
%         case 'x' ; w = -1;         
%         case 'y' ; w = -1;    
%         case 'xx'; w =  0;
%         case 'xy'; w =  0;
%         case 'yy'; w =  0;
%         end 
%     case 'G'  
%         switch keyxy            
%         case 'x' ; w = [];    
%         case 'y' ; w = [];      
%         case 'xx'; w = [];
%         case 'xy'; w = [];
%         case 'yy'; w = [];
%         end           
% 	case 'f'    
%         switch keyxy
%         case 'x' ; w = y*((-0.5<=x && x<=-0.25) || (0.25<=x && x<=0.5));    
%         case 'y' ; w = (x+0.25)*(-0.5<=x && x<=-0.25)+(x-0.25)*(0.25<=x && x<=0.5);         
%         case 'xx'; w = 0;
%         case 'xy'; w = ((-0.5<=x && x<=-0.25) || (0.25<=x && x<=0.5));
%         case 'yy'; w = 0;
%         case 'yxx';  w = [];  
%         case 'yxy';  w = []; 
%         case 'yyy';  w = []; 
%         end           
% 	case 'g'   
%         switch keyxy
%         case 'x' ; w = [-1; 1; 0; 0];
%         case 'y' ; w = [ 0; 0; -1; 1];             
%         case 'xx'; w = [ 0; 0; 0;0];
%         case 'xy'; w = [ 0; 0; 0;0];
%         case 'yy'; w = [ 0; 0; 0;0];
%         case 'yxx';  w = [];  
%         case 'yxy';  w = []; 
%         case 'yyy';  w = [];             
%         end        
%    end   
% end
% end

% function w = fxy(x,y,keyxy)
% 
% syms x1  x2;
% a=-0.5001;
% b=0.50001;
% 
% if nargin<3
%     fw = piecewise(          x1 <= a,     10^10,...
%                        a  <  x1 <-0.25, (x1+0.25)*x2,...
%                     -0.25 <= x1 <= 0.25,  0,...
%                      0.25 <  x1 < b,    (x1-0.25)*x2,...
%                        b  <= x1 ,         10^10);
% else
%     switch keyxy
%         case 'x' ; fw = piecewise(          x1 <= a,     0,...
%                                       a  <  x1 <-0.25,   x2,...
%                                    -0.25 <= x1 <= 0.25,  0,...
%                                     0.25 <  x1 < b,      x2,...
%                                       b  <= x1 ,         0);
%         case 'y' ; fw = piecewise(          x1 <= a,     0,...
%                                       a  <  x1 <-0.25,   (x1+0.25),...
%                                    -0.25 <= x1 <= 0.25,  0,...
%                                     0.25 <  x1 < b,      (x1-0.25),...
%                                       b  <= x1 ,         0); 
%         case 'xx'; fw = 0;  
%         case 'xy'; fw = piecewise(          x1 <= a,     0,...
%                                       a  <  x1 <-0.25,   1,...
%                                    -0.25 <= x1 <= 0.25,  0,...
%                                     0.25 <  x1 < b,      1,...
%                                       b  <= x1 ,         0);  
%         case 'yy'; fw = 0;
%      end   
% end
% 
% w =  subs(fw, [x1; x2], [x;y]) ;
% 
% end

 