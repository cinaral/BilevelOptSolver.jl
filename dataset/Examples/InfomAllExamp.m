function [probname, dim, xy, Ff]=InfomAllExamp(ProbNoName)
% This file provides information of all examples including:  
%  (a) Dimensions (n_x,n_y,n_G,n_g)
%  (b) Best known objectives F and f:  F = Ff(1) f = Ff(2)  
%  (c) Status F(3) of best known objectives: 
%             Ff(3)==1 global optimal 
%             Ff(3)==2 best known (may not be global optimal) 
%             Ff(3)==0 unknown. In such case  Ff = f = NaN
%  (d) Initial points xy=(x^0;y^0). Part of them are taken from 
%      current publications, and the rest are set as xy=(1,1,...,1)';
%
%
% Input:
%      ProbNoName: the number or the name of an input example.
%                  For example, ProbNoName=3 or ProbNoName='AnEtal2009'
%                  represents the example AnEtal2009.
% Outputs:
%      probname:   the name of the problem for all example except for 2 examples, 
%                        'OptimalControl' and 'ShehuEtal2019Ex42'
%                  For these two examples,  probname is their function handle
%                  probname= @(x,y,keyf,keyxy)OptimalControl(x,y,keyf,keyxy,data);
%                  probname= @(x,y,keyf,keyxy)ShehuEtal2019Ex42(x,y,keyf,keyxy,data);
%      dim:        dimensions of the input example namely, dim=[n_x,n_y,n_G,n_g]
%      xy:         Starting points xy=(x^0;y^0)\in\R^{n_x+n_y,1}
%      Ff:         =[F,f,Statu], where F and f are best known upper and lower level objectives
%                  function values. Statu is descibed as above in (c).

switch ProbNoName
%%%%%%              Nonlinear bilevel examples                       %%%%%%
    case {1, 'AiyoshiShimizu1984Ex2'}
        probname    = 'AiyoshiShimizu1984Ex2';
        dim        = [2 2 5 6];                  
        xy         = [10 10 20 20]' ; 
        Ff         = [5 0 1];
    case {2,'AllendeStill2013'}       
        probname    = 'AllendeStill2013';
        dim        = [2 2 5 2];                  
        xy         = [0 0 0 0]'; 
        Ff         = [1 -0.5 1 ];
    case {3,'AnEtal2009'}
        probname    = 'AnEtal2009';
        dim        = [2 2 6 4];                  
        xy         = [1 1 1 1]'; 
        Ff         = [2251.55 565.78 1];   
    case {4,'Bard1988Ex1'}
        probname    = 'Bard1988Ex1';  
        dim        = [1 1 1 4];                  
        xy         = [4 0]';
        Ff         = [17 1 1];
    case {5,'Bard1988Ex2'}
        probname    = 'Bard1988Ex2' ;
        dim        = [4 4 9 12];  
        xy         = [5 5 15 15 0 0 0 0]';
        Ff         = [-6600 54 1];
    case {6, 'Bard1988Ex3'}
        probname    = 'Bard1988Ex3'; 
        dim        = [2 2 3 4];                  
        xy         = [0 2 4 1]';
        Ff         = [-12.68 -1.02 1];
    case {7,'Bard1991Ex1' }     
        probname    = 'Bard1991Ex1';
        dim        = [1 2 2 3];                  
        xy         = [1 1 1]';
        Ff         = [2 12 1];
    case {8,'BardBook1998'}
        probname    = 'BardBook1998';
        dim        = [2 2 4 7];                  
         xy         = [1 1 1 1]'; 
        Ff         = [0 5 1];    
    case {9,'CalamaiVicente1994a'}
        probname    = 'CalamaiVicente1994a' ;
        dim        = [1 1 0 3];  
        xy         = [1 1]';
        Ff         = [0 0 1];       
    case {10,'CalamaiVicente1994b'}
        probname    = 'CalamaiVicente1994b' ;
        dim        = [4 2 0 6];  
        xy         = [1 1 1 1 1 1 1]';
        Ff         = [0.3125 -0.4063 1];    
    case {11,'CalamaiVicente1994c'}
        probname    = 'CalamaiVicente1994c' ;
        dim        = [4 2 0 6];  
        xy         = [1 1 1 1 1 1]';  
        Ff         = [0.3125 -0.4063 1];    
    case {12,'CalveteGale1999P1'}
        probname    = 'CalveteGale1999P1' ;
        dim        = [2 3 2 6];  
        xy         = [0 0.5 0 0.5 0]';
        Ff         = [-29.2 0.31 1];          
    case {13,'ClarkWesterberg1990a'}
        probname    = 'ClarkWesterberg1990a' ;
        dim        = [1 1 2 3];  
        xy         = [1 1]';
        Ff         = [5 4 1];
   case {14,'Colson2002BIPA1'}
        probname    = 'Colson2002BIPA1'; 
        dim        = [1 1 3 3 ];                  
        xy         = [10 10]';
        Ff         = [250 0 1];
    case {15,'Colson2002BIPA2'}
        probname    = 'Colson2002BIPA2' ;
        dim        = [1 1 1 4];  
        xy         = [4 0]';
        Ff         = [17 2 2];
    case {16,'Colson2002BIPA3'}
        probname    = 'Colson2002BIPA3' ;
        dim        = [1 1 2 2];  
        xy         = [4 0]';
        Ff         = [2 24.02 2];
    case {17,'Colson2002BIPA4'}
        probname    = 'Colson2002BIPA4' ;
        dim        = [1 1 2 2];  
        xy         = [1.5 2.25]';
        Ff         = [88.79 -0.77 2];
    case {18,'Colson2002BIPA5'}
        probname    = 'Colson2002BIPA5' ;
        dim        = [1 2 1 6];  
        xy         = [1 1 1]';  
        Ff         = [2.75 0.57 2];
    case {19,'Dempe1992a'}
        probname    = 'Dempe1992a' ;
        dim        = [2 2 1 2];  
        xy         = [1 1 1 1]'; 
        Ff         = [NaN NaN 0];           
    case {20,'Dempe1992b'}
        probname    = 'Dempe1992b' ;
        dim        = [1 1 0 1];  
        xy         = -[1 1]'; 
        Ff         = [31.25 4 1];   
    case {21,'DempeDutta2012Ex24'}
        probname    = 'DempeDutta2012Ex24' ;
        dim        = [1 1 0 1];  
        xy         = [1 1]';
        Ff         = [0 0 1];
    case {22,'DempeDutta2012Ex31'}
        probname    = 'DempeDutta2012Ex31' ;
        dim        = [2 2 4 2];  
        xy         = [1 1 1 1]'; 
        Ff         = [-1 4 1];
    case {23,'DempeEtal2012'}
        probname    = 'DempeEtal2012' ;
        dim        = [1 1 2 2];  
        xy         = [1 1]';  
        Ff         = [-1 -1 1];        
    case {24,'DempeFranke2011Ex41'}
        probname    = 'DempeFranke2011Ex41' ;
        dim        = [2 2 4 4 ];  
        xy         = [1 1 1 1]';  
        Ff         = [5  -2 1];         
    case {25,'DempeFranke2011Ex42'}
        probname    = 'DempeFranke2011Ex42' ;
        dim        = [2 2 4 3 ];  
        xy         = [1 1 1 1]'; 
        Ff         = [2.13 -3.5 1];         
    case {26,'DempeFranke2014Ex38'}
        probname    = 'DempeFranke2014Ex38' ;
        dim        = [2 2 4 4 ];  
        xy         = [1 1 1 1]';  
        Ff         = [-1 -4 1];               
    case {27,'DempeLohse2011Ex31a'} 
        probname    = 'DempeLohse2011Ex31a' ;
        dim        = [2 2 0 4 ];  
        xy         = [1 1 1 1]';  
        Ff         = [-5.5 0 1];                 
    case {28,'DempeLohse2011Ex31b'}
        probname    = 'DempeLohse2011Ex31b' ;
        dim        = [3 3 0 5];  
        xy         = [1 1 1 1 1 1]'; 
        Ff         = [-12 0 1];          
    case {29,'DeSilva1978'}
        probname    = 'DeSilva1978' ;
        dim        = [2 2 0 4];  
        xy         = [0 0 1 1]';
        Ff         = [-1 0 1];                  
    case {30,'FalkLiu1995'}
        probname    = 'FalkLiu1995' ;
        dim        = [2 2 0 4];  
        xy         = [0 0 1 1]';  
        Ff         = [ -2.1962  0 1];
    case {31,'FloudasEtal2013'}
        probname    = 'FloudasEtal2013' ;
        dim        = [2 2 4 7];  
        xy         = [10 10 20 20]';  
        Ff         = [0 200 1];                
    case {32,'FloudasZlobec1998'}
        probname    = 'FloudasZlobec1998' ;
        dim        = [1 2 2 6];  
        xy         = [1 1 1]'; 
        Ff         = [1 -1 1];          
    case {33,'GumusFloudas2001Ex1'}
        probname    = 'GumusFloudas2001Ex1' ;
        dim        = [1 1 3 3];  
        xy         = [1 1]'; 
        Ff         = [2250 197.75 1];         
    case {34,'GumusFloudas2001Ex3'}
        probname    = 'GumusFloudas2001Ex3' ;
        dim        = [2 3 4 9];  
        xy         = [0 1/2 0 1/2 0]';
        Ff         = [-29.2 0.31 1];          
    case {35,'GumusFloudas2001Ex4'}
        probname    = 'GumusFloudas2001Ex4' ;
        dim        = [1 1 5 2];  
        xy         = [1 1]';
        Ff         = [9 0 1];          
    case {36,'GumusFloudas2001Ex5'}
        probname    = 'GumusFloudas2001Ex5' ;
        dim        = [1 2 2 6];  
        xy         = [1 1 1]';
        Ff         = [0.194 -7.23 1];           
    case {37,'HatzEtal2013'}
        probname    = 'HatzEtal2013' ;
        dim        = [1 2 0 2];  
        xy         = [1 1 1]';  
        Ff         = [0 0 1];          
    case {38,'HendersonQuandt1958'}
        probname    = 'HendersonQuandt1958' ;
        dim        = [1 1 2 1];  
        xy         = [0 0]';
        Ff         = [-3266.67 -711.11 2];          
    case {39,'HenrionSurowiec2011'}
        probname    = 'HenrionSurowiec2011' ;
        dim        = [1 1 0 0];  
        xy         = [1 1]'; 
        Ff         = [0 0 1];         
    case {40,'IshizukaAiyoshi1992a'}
        probname    = 'IshizukaAiyoshi1992a' ;
        dim        = [1 2 1 5];  
        xy         = [1 1 1]' ;
        Ff         = [0 -1.5 1];                  
    case {41,'KleniatiAdjiman2014Ex3'}
        probname    = 'KleniatiAdjiman2014Ex3' ;
        dim        = [1 1 2 2];  
        xy         = [1 1]'; 
        Ff         = [-1 0 1];           
    case {42,'KleniatiAdjiman2014Ex4'}
        probname    = 'KleniatiAdjiman2014Ex4' ;
        dim        = [5 5 13 11];  
        xy         = ones(10,1); 
        Ff         = [-10  -3.1 2];   
    case {43, 'LamparielloSagratella2017Ex23'}
        probname    = 'LamparielloSagratella2017Ex23' ;
        dim        = [1 2 2 2];  
        xy         = [1 1 1]';
        Ff         = [-1  1 1]; 
    case {44,'LamparielloSagratella2017Ex31'}
        probname    = 'LamparielloSagratella2017Ex31' ;
        dim        = [1 1 1 1];  
        xy         = [1 1]'; 
        Ff         = [1 0 1];           
    case {45,'LamparielloSagratella2017Ex32' }
        probname    = 'LamparielloSagratella2017Ex32' ;
        dim        = [1 1 0 0];  
        xy         = [1 1]';  
        Ff         = [0.5  0 1];          
    case {46,'LamparielloSagratella2017Ex33'}
        probname    = 'LamparielloSagratella2017Ex33' ;
        dim        = [1 2 1 3];  
        xy         = [1 1 1]';
        Ff         = [0.5  0 1];           
    case {47,'LamparielloSagratella2017Ex35'}
        probname    = 'LamparielloSagratella2017Ex35' ;
        dim        = [1 1 2 3];  
        xy         = [1 1]'; 
        Ff         = [0.80 -0.40 1];                   
    case {48,'LucchettiEtal1987'}
        probname    = 'LucchettiEtal1987' ;
        dim        = [1 1 2 2];  
        xy         = [1 1]'; 
        Ff         = [0  0 1];           
    case {49,'LuDebSinha2016a'}
        probname    = 'LuDebSinha2016a' ;
        dim        = [1 1 4 0];  
        xy         = [0 2]'; 
        Ff         = [1.1360 1.1838 2];       %0.8395 1    
    case {50,'LuDebSinha2016b'}
        probname    = 'LuDebSinha2016b' ;
        dim        = [1 1 4 0];  
        xy         = [1 1]'; 
        Ff         = [0 1.6645 2];            %0.3077 0.9945      
    case {51,'LuDebSinha2016c'}
        probname    = 'LuDebSinha2016c' ;
        dim        = [1 1 4 0];  
        xy         = [1 1]';  
        Ff         = [1.12 0.06 2];           
    case {52,'LuDebSinha2016d'}
        probname    = 'LuDebSinha2016d' ;
        dim        = [2 2 11 3];  
        xy         = [1 1 1 1]'; 
        Ff         = [NaN NaN 0 ];          
    case {53,'LuDebSinha2016e'}
        probname    = 'LuDebSinha2016e' ;
        dim        = [1 2 6 3];  
        xy         = [1 1 1]';  
        Ff         = [NaN NaN 0];          
    case {54,'LuDebSinha2016f'}
        probname    = 'LuDebSinha2016f' ;
        dim        = [2 1 9 0];  
        xy         = [1 1 1]'; 
        Ff         = [NaN NaN 0];          
    case {55, 'MacalHurter1997'}
        probname    = 'MacalHurter1997' ;
        dim        = [1 1 0 0];  
        xy         = [1 1]'; 
        Ff         = [81.327 -0.333 1];          
    case {56,'Mirrlees1999'}
        probname    = 'Mirrlees1999' ;
        dim        = [1 1 0 2];  
        xy         = [1 1]';  
        Ff         = [1.002 -1.02 1];            
    case {57,'MitsosBarton2006Ex38'}
        probname    = 'MitsosBarton2006Ex38' ;
        dim        = [1 1 4 2];  
        xy         = [1 1]';
        Ff         = [0 0 1];
    case {58,'MitsosBarton2006Ex39' }
        probname    = 'MitsosBarton2006Ex39' ;
        dim        = [1 1 3 2];  
        xy         = [1 1]';
        Ff         = [-1 -1 1];         
    case {59,'MitsosBarton2006Ex310'}
        probname    = 'MitsosBarton2006Ex310' ;
        dim        = [1 1 2 2];  
        xy         = [1 1]';  
        Ff         = [0.5 -0.1 1];         
    case {60,'MitsosBarton2006Ex311'}
        probname    = 'MitsosBarton2006Ex311' ;
        dim        = [1 1 2 2];  
        xy         = [1 1]';  
        Ff         = [-0.80 0.00 1];          
    case {61,'MitsosBarton2006Ex312'}
        probname    = 'MitsosBarton2006Ex312' ;
        dim        = [1 1 2 2];  
        xy         = [1 1]'; 
        Ff         = [0 0 1];          
     case {62,'MitsosBarton2006Ex313'}
        probname    = 'MitsosBarton2006Ex313' ;
        dim        = [1 1 2 2];  
        xy         = [1 1]';  
        Ff         = [-1 0 1];        
    case {63,'MitsosBarton2006Ex314'}
        probname    = 'MitsosBarton2006Ex314' ;
        dim        = [1 1 2 2];  
        xy         = [1 1]'; 
        Ff         = [1/4 -1/12 1];                 
    case {64,'MitsosBarton2006Ex315'}
        probname    = 'MitsosBarton2006Ex315' ;
        dim        = [1 1 2 2];  
        xy         = [1 1]';  
        Ff         = [0 -0.83 1];          
    case {65,'MitsosBarton2006Ex316'}
        probname    = 'MitsosBarton2006Ex316' ;
        dim        = [1 1 2 2];  
        xy         = [1 1]';  
        Ff         = [-2 0 1];         
    case {66,'MitsosBarton2006Ex317'}
        probname    = 'MitsosBarton2006Ex317' ;
        dim        = [1 1 2 2];  
        xy         = [1 1]';   
        Ff         = [3/16 -1/64 1];         
    case {67,'MitsosBarton2006Ex318'}
        probname    = 'MitsosBarton2006Ex318' ;
        dim        = [1 1 2 2];  
        xy         = [1 1]'; 
        Ff         = [-1/4 0 1]; 
    case {68,'MitsosBarton2006Ex319'}
        probname    = 'MitsosBarton2006Ex319' ;
        dim        = [1 1 2 2];  
        xy         = [-1 1]';
        Ff         = [-0.258 -0.0178 1];
    case {69,'MitsosBarton2006Ex320'}
        probname    = 'MitsosBarton2006Ex320' ;
        dim        = [1 1 2 2];  
        xy         = [1 1]'; 
        Ff         = [0.3125 -0.0833 1];        
    case {70,'MitsosBarton2006Ex321'}
        probname    = 'MitsosBarton2006Ex321' ;
        dim        = [1 1 2 2];  
        xy         = [1 1]';  
        Ff         = [0.2095 -0.0656 1];        
    case {71,'MitsosBarton2006Ex322'}
        probname    = 'MitsosBarton2006Ex322' ;
        dim        = [1 1 2 3];  
        xy         = [1 1]';   
        Ff         = [0.2095 -0.0656 1];        
    case {72,'MitsosBarton2006Ex323'}
        probname    = 'MitsosBarton2006Ex323' ;
        dim        = [1 1 3 3];  
        xy         = [0 1]';   
        Ff         = [0.176 -1 1];        
    case {73,'MitsosBarton2006Ex324'}
        probname    = 'MitsosBarton2006Ex324' ;
        dim        = [1 1 2 2];  
        xy         = [1 1]'; 
        Ff         = [-1.755 0 1];         
    case {74,'MitsosBarton2006Ex325'}
        probname    = 'MitsosBarton2006Ex325' ;
        dim        = [2 3 6 9];  
        xy         = [1 1 1 1 1]'; 
        Ff         = [-1 -2 2];        
    case {75,'MitsosBarton2006Ex326'}
        probname    = 'MitsosBarton2006Ex326' ;
        dim        = [2 3 7 6];  
        xy         = [1 1 1 1 1]';  
        Ff         = [-2.354 -2 1];         
    case {76,'MitsosBarton2006Ex327'}
        probname    = 'MitsosBarton2006Ex327' ;
        dim        = [5 5 13 13];  
        xy         = ones(10,1);  
        Ff         = [2 -1.1 2];        
    case {77,'MitsosBarton2006Ex328'}
        probname    = 'MitsosBarton2006Ex328' ;
        dim        = [5 5 13 13];  
        xy         = ones(10,1); 
        Ff         = [-10 -3.1 2];           
    case {78,'MorganPatrone2006a'}
        probname    = 'MorganPatrone2006a' ;
        dim        = [1 1 2 2];  
        xy         = [1 1]';  
        Ff         = [-1 0 1];         
    case {79,'MorganPatrone2006b'}
        probname    = 'MorganPatrone2006b' ;
        dim        = [1 1 0 4];  
        xy         = [1 1]';
        Ff         = [-1.25 0 1];          
     case {80,'MorganPatrone2006c'}
        probname    = 'MorganPatrone2006c' ;
        dim        = [1 1 0 4];  
        xy         = [1 1]';  
        Ff         = [-1  -0.25 1];         
    case {81,'MuuQuy2003Ex1'}
        probname    = 'MuuQuy2003Ex1' ;
        dim        = [1 2 2 3];  
        xy         = [2 0 0]'; 
        Ff         = [-2.0769  -0.5868 2];           
    case {82,'MuuQuy2003Ex2'}
        probname    = 'MuuQuy2003Ex2' ;
        dim        = [2 3 3 4 ];  
        xy         = [0 1 1 1 1 ]'; 
        Ff         = [0.6426 1.6708 2];         
    case {83,'NieWangYe2017Ex34'}
        probname    = 'NieWangYe2017Ex34' ;
        dim        = [1 2 2 2 ];  
        xy         = [1 1 1]'; 
        Ff         = [2 0 1];            
    case {84,'NieWangYe2017Ex52'}
        probname    = 'NieWangYe2017Ex52' ;
        dim        = [2 3 5 2];  
        xy         = [1 1 1 1 1]';
        Ff         = [-1.710 -2.232 1]; 
    case {85,'NieWangYe2017Ex54'}
        probname    = 'NieWangYe2017Ex54' ;
        dim        = [4 4 3 2];  
        xy         = [1 1 1 1 1 1 1 1]';  
        Ff         = [-0.437 -1.190 1];        
    case {86,'NieWangYe2017Ex57'}
        probname    = 'NieWangYe2017Ex57' ;
        dim        = [2 3 5 2 ];  
        xy         = [1 1 1 1 1]'; 
        Ff         = [-2 -1 2];           
    case {87,'NieWangYe2017Ex58'}
        probname    = 'NieWangYe2017Ex58' ;
        dim        = [4 4 3 2];  
        xy         = [0.5442  0.4682  0.4904  0.4942 -0.7792 -0.5034 -0.2871 -0.1855]' ;
        Ff         = [-3.488 -0.862 2];           
    case {88,'NieWangYe2017Ex61'}
        probname    = 'NieWangYe2017Ex61' ;
        dim        = [2 2 5 1 ];  
        xy         = [1 -1 -1 1]';
        Ff         = [-1.022 -1.084 2];          
    case {89,'Outrata1990Ex1a'}
        probname    = 'Outrata1990Ex1a' ;
        dim        = [2 2 0 4 ];  
        xy         = [0 0 3 3]'; 
        Ff         = [-8.92 -6.05 2];          
    case {90,'Outrata1990Ex1b'}
        probname    = 'Outrata1990Ex1b' ;
        dim        = [2 2 0 4 ];  
        xy         = [0 0 3 3]';
        Ff         = [-7.56 -0.58 2];  
    case {91,'Outrata1990Ex1c'}
        probname    = 'Outrata1990Ex1c' ;
        dim        = [2 2 0 4 ];  
        xy         = [0 0 3 3]';
        Ff         = [-12  -112.71 2];         
    case {92,'Outrata1990Ex1d'}
        probname    = 'Outrata1990Ex1d' ;
        dim        = [2 2 0 4 ];  
        xy         = [6 -3  3 3]';
        Ff         = [-3.60 -2 2];        
    case {93,'Outrata1990Ex1e'}
        probname    = 'Outrata1990Ex1e' ;
        dim        = [2 2 0 4 ];  
        xy         = [6 -3  3 3]';
        Ff         = [-3.15 -16.29 2];         
     case {94,'Outrata1990Ex2a'}
        probname    = 'Outrata1990Ex2a' ;
        dim        = [1 2 1 4];  
        xy         = [1 1 1]';
        Ff         = [0.50 -14.53 2];          
      case {95,'Outrata1990Ex2b'}
        probname    = 'Outrata1990Ex2b' ;
        dim        = [1 2 1 4];  
        xy         = [5 0 0]';
        Ff         = [0.50 -4.50 2];   
     case {96,'Outrata1990Ex2c' }
        probname    = 'Outrata1990Ex2c' ;
        dim        = [1 2 1 4];  
        xy         = [0 0 0]'; 
        Ff         = [1.860  -10.931 2];           
     case {97,'Outrata1990Ex2d'}
        probname    = 'Outrata1990Ex2d' ;
        dim        = [1 2 1 4];  
        xy         = [0 0 0]';
        Ff         = [0.92  -19.47 2];            
     case {98,'Outrata1990Ex2e'}
        probname    = 'Outrata1990Ex2e' ;
        dim        = [1 2 1 4];  
        xy         = [0 0 0]';
        Ff         = [0.90 -14.94 2];           
     case {99,'Outrata1993Ex31'}
        probname    = 'Outrata1993Ex31' ;
        dim        = [1 2 1 4];  
        xy         = [0 0 0]';
        Ff         = [1.56 -11.68 2];         
     case {100,'Outrata1993Ex32'} 
        probname    = 'Outrata1993Ex32' ;
        dim        = [1 2 1 4];  
        xy         = [0 0 0]';
        Ff         = [3.208  -20.531 2];                
     case {101,'Outrata1994Ex31'}
        probname    = 'Outrata1994Ex31' ;
        dim        = [1 2 2 4];  
        xy         = [0 0 0]';
        Ff         = [3.208  -20.531 2];            
     case {102,'OutrataCervinka2009'}
        probname    = 'OutrataCervinka2009' ;
        dim        = [2 2 1 3];  
        xy         = [1 1 1 1]';
        Ff         = [0 0 1];               
     case {103,'PaulaviciusAdjiman2017a'}
        probname    = 'PaulaviciusAdjiman2017a' ;
        dim        = [1 1 4 2];  
        xy         = [1 1]';
        Ff         = [0.25 0 1];            
     case {104,'PaulaviciusAdjiman2017b'} 
        probname    = 'PaulaviciusAdjiman2017b' ;
        dim        = [1 1 4 2];  
        xy         = [1 1]';
        Ff         = [-2 -1.5 1];    
     case {105,'SahinCiric1998Ex2' }
        probname    = 'SahinCiric1998Ex2' ;
        dim        = [1 1 2 3];  
        xy         = [1  1]'; 
        Ff         = [5 4 1];    
     case {106,'ShimizuAiyoshi1981Ex1'}
        probname    = 'ShimizuAiyoshi1981Ex1' ;
        dim        = [1 1 3 3];  
        xy         = [1 1]';
        Ff         = [100 0 1];            
     case {107,'ShimizuAiyoshi1981Ex2'} 
        probname    = 'ShimizuAiyoshi1981Ex2' ;
        dim        = [2 2 3 4];  
        xy         = [1 1 1 1]';
        Ff         = [225 100 1];                    
      case {108,'ShimizuEtal1997a'}
        probname    = 'ShimizuEtal1997a' ;
        dim        = [1 1 0 3];  
        xy         = [1 1]';
        Ff         = [NaN NaN 0];                    
     case {109,'ShimizuEtal1997b'}
        probname    = 'ShimizuEtal1997b' ;
        dim        = [1 1 2 2];  
        xy         = [1 1]'; 
        Ff         = [2250  197.75 1];                  
     case {110,'SinhaMaloDeb2014TP3'}
        probname    = 'SinhaMaloDeb2014TP3' ;
        dim        = [2 2 3 4];  
        xy         = [1 1 1 1]' ; 
        Ff         = [-18.679  -1.016 2];           
     case {111,'SinhaMaloDeb2014TP6'}
        probname    = 'SinhaMaloDeb2014TP6' ;
        dim        = [1 2 1 6];  
        xy         = [1 1 1]';
        Ff         = [-1.209 7.615 2];            
     case {112,'SinhaMaloDeb2014TP7'}
        probname    = 'SinhaMaloDeb2014TP7' ;
        dim        = [2 2 4 4];  
        xy         = [1 1 0 1]';
        Ff         = [-1.961 1.961 2];           
     case {113,'SinhaMaloDeb2014TP8'}
        probname    = 'SinhaMaloDeb2014TP8' ;
        dim        = [2 2 5 6];  
        xy         = [1 1 1 1]';
        Ff         = [0  100 1];                              
     case {114,'SinhaMaloDeb2014TP9'}
        probname    = 'SinhaMaloDeb2014TP9' ;
        dim        = [10 10 0 20];  
        xy         = 0*ones(20,1);
        Ff         = [0 1 2];            
     case {115,'SinhaMaloDeb2014TP10'}
        probname    = 'SinhaMaloDeb2014TP10' ;
        dim        = [10 10 0 20];  
        xy         = 0*ones(20,1);
        Ff         = [0 1 2];            
     case {116,'TuyEtal2007' }
        probname    = 'TuyEtal2007' ;
        dim        = [1 1 2 3];  
        xy         = [5 0]';
        Ff         = [22.5  -1.5 1];           
     case {117,'Vogel2012' }
        probname    = 'Vogel2012' ;
        dim        = [1 1 2 1];  
        xy         = [1 1]'; 
        Ff         = [1 -2 1];            
     case {118,'WanWangLv2011'}
        probname    = 'WanWangLv2011' ;
        dim        = [2 3 0 8];  
        xy         = [0 0.5 0 0 0 ]';
        Ff         = [10.62 -0.50 1];     
     case {119,'YeZhu2010Ex42'}
        probname    = 'YeZhu2010Ex42' ;
        dim        = [1 1 2 1];  
        xy         = -[1 1]';
        Ff         = [1 -2 1];     
	case {120,'YeZhu2010Ex43'}
        probname    = 'YeZhu2010Ex43' ;
        dim        = [1 1 2 1];  
        xy         = -[1 1]';
        Ff         = [1.25 -2 1];            
	case {121,'Yezza1996Ex31'}
        probname    = 'Yezza1996Ex31' ;
        dim        = [1 1 2 2];  
        xy         = [1 1]';
        Ff         = [1.5 -2.5 1];           
	case {122,'Yezza1996Ex41'}
        probname    = 'Yezza1996Ex41' ;
        dim        = [1 1 0 2];  
        xy         = -[1 1]';  
        Ff         = [0.5 2.5 1];            
	case {123,'Zlobec2001a'}
        probname    = 'Zlobec2001a' ;
        dim        = [1 2 0 3];  
        xy         = [1 1 1]';  
        Ff         = [-1 -1 1];          
	case {124,'Zlobec2001b'}
        probname    = 'Zlobec2001b' ;
        dim        = [1 1 2 4];  
        xy         = [1 1]'; 
        Ff         = [NaN NaN 0];          
	case {125 ,'DesignCentringP1'}       
        probname    = 'DesignCentringP1' ;
        dim        = [3 6 3 3];  
        xy         = ones(9,1); 
        Ff         = [NaN NaN 0];  
	case {126, 'DesignCentringP2'}       
        probname    = 'DesignCentringP2' ;
        dim        = [4 6 5 3];  
        xy         = 10*ones(10,1); 
        Ff         = [NaN NaN 0]; 
	case {127,'DesignCentringP3'}
        probname    = 'DesignCentringP3' ;
        dim        = [6 6 3 3];  
        xy         = [1 1 1 0 0 1 1 1 1 1 1 1]'; 
        Ff         = [NaN NaN 0]; 
	case {128,'DesignCentringP4'}
        probname    = 'DesignCentringP4' ;
        dim        = [4 6 3 12];  
        xy         = ones(10,1); 
        Ff         = [NaN NaN 0];   
	case {129,'NetworkDesignP1'}
        probname    = 'NetworkDesignP1' ;
        dim        = [5 5 5 11];  
        xy         = zeros(10,1);  
        Ff         = [300.5 419.8 2];                
     case {130,'NetworkDesignP2'}
        probname    = 'NetworkDesignP2' ;
        dim        = [5 5 5 11];  
        xy         = zeros(10,1);
        Ff         = [142.9 81.95 2];  
     case {131,'RobustPortfolioP1'}
        probname    = 'RobustPortfolioP1' ;
        dim        = [11 10 13 11];  %N=10,50,100..., dim=[N+1 N N+3 N+1]
        xy         = ones(dim(1)+dim(2),1); 
        Ff         = [1.15  0 2];     
     case {132,'RobustPortfolioP2'}
        probname    = 'RobustPortfolioP2' ;
        dim        = [11 10 13 11];  %N=10,50,100..., dim=[N+1 N N+3 N+1]
        xy         = ones(dim(1)+dim(2),1); 
        Ff         = [1.15  0 2];     
     case {133,'TollSettingP1'}
        probname    = 'TollSettingP1' ;
        dim        = [3 8 3 18];   
        xy         = [0; 5; 5; 0; 1; 0; 0; 0; 1; 0; 0]; 
        Ff         = [-7 12 2];      
     case {134,'TollSettingP2'}
        probname    = 'TollSettingP2' ;
        dim        = [3 18 3 38];   
        xy         = [0;0;0;0; 0; 0; 10; 1; 0; 1; 0; 0; 0; 10; 0; 0; 0; 1; 0; 0; 10];
        Ff         = [-4.5 32 2];     
     case {135,'TollSettingP3'}
        probname    = 'TollSettingP3' ;
        dim        = [3 18 3 38];   
        xy         = [0;0;0;0; 0; 0; 10; 1; 0; 1; 0; 0; 0; 10; 0; 0; 0; 1; 0; 0; 10];
        Ff         = [-3.5 32 2];     
     case {136,'TollSettingP4'}
        probname    = 'TollSettingP4' ;
        dim        = [2 4 0 8];   
        xy         = ones(6,1); 
        Ff         = [-4 14 2];    
     case {137,'TollSettingP5'}
        probname    = 'TollSettingP5' ;
        dim        = [1 4 0 8];   
        xy         = ones(5,1); 
        Ff         = [-2.5 14 2];   
      case {138,'OptimalControl'}
        rate       = 0.3;  %smaller 'rate' is, larger size of 'dim'
        data       = OptimalControl_data(rate);
        dim        = [length(data.k) 2*data.ni 3 4*data.ni]; 
        xy         = -ones(dim(1)+dim(2),1); 
        Ff         = [NaN NaN 0];
        disp('Example ''OptimalControl'' with extra input data');
        probname    = @(x,y,keyf,keyxy)OptimalControl(x,y,keyf,keyxy,data);
        
%%%%%%              Linear bilevel examples                          %%%%%%
    case {139,  'AnandalinghamWhite1990'}
        probname    = 'AnandalinghamWhite1990';
        dim        = [1 1 1 6];                   
        xy         = [1 1]';   
        Ff         = [-49 15 1]; 
    case {140,  'Bard1984a'  }   
        probname    = 'Bard1984a';
        dim        = [1 1 1 5] ;                  
        xy         = [1 1]'; 
        Ff         = [28/9 -60/9 1];
    case {141,'Bard1984b'}
        probname    = 'Bard1984b';
        dim        = [1 1 1 5] ;                  
        xy         = [1 1]'; 
        Ff         = [-37.6 1.6 1];   
    case {142,'Bard1991Ex2'}
        probname    = 'Bard1991Ex2';  
        dim        = [1 2 1 5];                  
        xy         = [0 0 0]';
        Ff         = [-1 -1 1];
    case {143,'BardFalk1982Ex2'}
        probname    = 'BardFalk1982Ex2' ;
        dim        = [2 2 2 5];  
        xy         = [1 1 1 1]';
        Ff         = [-3.25 -4 1];
    case {144, 'BenAyedBlair1990a'}
        probname    = 'BenAyedBlair1990a'; 
        dim        = [1 2 2 4];                  
        xy         = [1 1 1]';
        Ff         = [-2.5 -5 1];
    case {145, 'BenAyedBlair1990b'}      
        probname    = 'BenAyedBlair1990b';
        dim        = [1 1 1 4];                  
        xy         = [3 4]';
        Ff         = [-6 5 1];
    case {146,'BialasKarwan1984a'}
        probname    = 'BialasKarwan1984a';
        dim        = [1 2 1 7];                  
        xy         = [1 1 1]'; 
        Ff         = [-2 -0.5 1];    
    case {147, 'BialasKarwan1984b'}
        probname    = 'BialasKarwan1984b';
        dim        = [1 1 1 6];                  
        xy         = [10 14]'; 
        Ff         = [-11 11 1];       
    case {148, 'CandlerTownsley1982'} 
        probname    = 'CandlerTownsley1982' ;
        dim        = [2 3 2 6];  
        xy         = [1 1 1 1 1]';
        Ff         = [-29.2 3.2 1];    
    case {149, 'ClarkWesterberg1988'} 
        probname    = 'ClarkWesterberg1988' ;
        dim        = [1 1 0 3];  
        xy         = [20 0]';  
        Ff         = [-37 14 1];    
    case {150,'ClarkWesterberg1990b' }
        probname    = 'ClarkWesterberg1990b' ;
        dim        = [1 2 2 5];  
        xy         = [1 1 1 ]';
        Ff         = [-13 -4 1];          
    case {151,'GlackinEtal2009'}
        probname    = 'GlackinEtal2009' ;
        dim        = [2 1 3 3];  
        xy         = [0 1 0]';
        Ff         = [6 0 1];
   case {152,'HaurieSavardWhite1990'}
        probname    = 'HaurieSavardWhite1990'; 
        dim        = [1 1 0 4 ];                  
        xy         = [1 1]';
        Ff         = [27 -3 1];
    case {153,'HuHuangZhang2009' }
        probname    = 'HuHuangZhang2009' ;
        dim        = [1 2 1 5];  
        xy         = [1 1 1]';
        Ff         = [-76/9 -41/9 1];
    case {154,'LanWenShihLee2007'}
        probname    = 'LanWenShihLee2007' ;
        dim        = [1 1 1 7];  
        xy         = [1 1]';
        Ff         = [-85.0855 50.174 2];
    case {155,'LiuHart1994' }
        probname    = 'LiuHart1994' ;
        dim        = [1 1 1 4];  
        xy         = [1 1]';
        Ff         = [-16 4 1];
    case {156,'MershaDempe2006Ex1'}
        probname    = 'MershaDempe2006Ex1' ;
        dim        = [1 1 1 5];  
        xy         = [1 1]'; 
        Ff         = [NaN NaN 0];
    case {157,'MershaDempe2006Ex2'}
        probname    = 'MershaDempe2006Ex2' ;
        dim        = [1 1 2 2];  
        xy         = [1 1 ]'; 
        Ff         = [-20 -6 1];           
    case {158,'TuyEtal1993' }
        probname    = 'TuyEtal1993' ;
        dim        = [2 2 3 4];  
        xy         = [0 0 0 0]'; 
        Ff         = [-3.25 -6 1];   
    case {159,'TuyEtal1994'}
        probname    = 'TuyEtal1994' ;
        dim        = [2 2 3 3];  
        xy         = [1 1 1 1]';
        Ff         = [6 0 1];
    case {160,'TuyEtal2007Ex3'}
        probname    = 'TuyEtal2007Ex3' ;
        dim        = [10 6 12 13]  ;  
        xy         = 10*ones(16,1); 
        Ff         = [-467.4613 -11.6194 2]; 
    case {161,'VisweswaranEtal1996'}
        probname    = 'VisweswaranEtal1996' ;
        dim        = [1 1 1 5];  
        xy         = [1 1]'; 
        Ff         = [28/9 -60/9 1];   
    case {162,'WangJiaoLi2005' }
        probname    = 'WangJiaoLi2005' ;
        dim        = [1 2 2 2];  
        xy         = [1 1 1]';  
        Ff         = [-1000 -1 1];        
        
%%%%%               Simple bilevel examples                          %%%%%%
    case {163,'FrankeEtal2018Ex53'}  
        probname    = 'FrankeEtal2018Ex53';
        dim        = [1 2 4 4];                   
        xy         = [1 1 1]';   
        Ff         = [1 1 1]; 
    case {164, 'FrankeEtal2018Ex511'}  
        probname    = 'FrankeEtal2018Ex511';
        dim        = [1 3 0 4] ;                  
        xy         = [1 1 1 1]'; 
        Ff         = [3 0 1];
    case {165,'FrankeEtal2018Ex513'}
        probname    = 'FrankeEtal2018Ex513';
        dim        = [1 3 0 3] ;                  
        xy         = [1 1 1 1]'; 
        Ff         = [-1 0 1];
    case {166,'FrankeEtal2018Ex521'}
        probname    = 'FrankeEtal2018Ex521';
        dim        = [1 2 0 3] ;                  
        xy         = [0 0 0]'; 
        Ff         = [-1 0 1];
    case {167,'MitsosBarton2006Ex31'} 
        probname    = 'MitsosBarton2006Ex31' ;
        dim        = [1 1 2 2];  
        xy         = [-1 -1]';
        Ff         = [1 -1 1];
    case {168,'MitsosBarton2006Ex32'}
        probname    = 'MitsosBarton2006Ex32' ;
        dim        = [1 1 3 2];  
        xy         = [-1 -1]';
        Ff         = [NaN NaN 0];
    case {169, 'MitsosBarton2006Ex33'}
        probname    = 'MitsosBarton2006Ex33' ;
        dim        = [1 1 2 3];  
        xy         = [-10 -10]';
        Ff         = [-1 1 1];
    case {170,'MitsosBarton2006Ex34'}
        probname    = 'MitsosBarton2006Ex34' ;
        dim        = [1 1 2 2];  
        xy         = [-0.5 -0.5]';
        Ff         = [1 -1 1];
    case {171,'MitsosBarton2006Ex35'}
        probname    = 'MitsosBarton2006Ex35' ;
        dim        = [1 1 2 2];  
        xy         = [-1 -1]';
        Ff         = [0.5 -1 1];      
    case {172,'MitsosBarton2006Ex36'}
        probname    = 'MitsosBarton2006Ex36' ;
        dim        = [1 1 2 2];  
        xy         = [1 1]';
        Ff         = [-1 -1 1];        
    case {173,'ShehuEtal2019Ex42'}
        ny         = 10;
        dim        = [1 ny 0 0]; 
        xy         = ones(1+ny,1);
        Ff         = [NaN NaN 0];  
        disp('Example ''ShehuEtal2019Ex42'' with extra input data');
        row        = ceil(ny/2);
        data       = ShehuEtal2019Ex42data(row,ny,1e-3);
        probname    = @(x,y,keyf,keyxy)ShehuEtal2019Ex42(x,y,keyf,keyxy,data);
end
end