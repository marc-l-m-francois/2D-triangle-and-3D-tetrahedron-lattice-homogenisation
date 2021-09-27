clear; clc; close all; 
%% Parameters Entry

% The Constructing parameters                                                            
coord_SNude=[0,0,0]; % Coordinate of the starting node    
Num_msh=9; % Number of elements per beam                      
L=10; %the length of the beam in (meter)                                   


% The three imposed tensors 
%         [ eps_11 0  0 ]         [ 0    0      0 ]         [ 0   0    0   ]
%    H_1= |  0    0  0 |  , H_2= | 0  eps_22   0 | ,  H_3= | 0   0    0   |
%         [  0     0  0 ]         [ 0   0       0 ]         [ 0   0 eps_33 ]
%
%         [ 0     eps_12  0 ]         [ 0    0       0  ]         [  0      0  eps_13 ]
%    H_4= |eps_12    0    0 |  , H_5= | 0    0   eps_23 | ,  H_6= |  0      0   0     |
%         [  0       0    0 ]         [ 0  eps_23    0  ]         [ eps_13  0   0     ]
%
%  H_1 -->  Tension or Comprision (+/-) in X direction
%  H_2 -->  Tension or Comprision (+/-) in Y direction
%  H_3 -->  Tension or Comprision (+/-) in Z direction
%  H_4 -->  Shear(Counter clock or clock wise +/-) in the XY direction
%  H_5 -->  Shear(Counter clock or clock wise +/-) in the YZ direction
%  H_6 -->  Shear(Counter clock or clock wise +/-) in the XZ direction

H=[1 0 0 ; 0 0 0 ; 0 0 0 ];    % Filling as choices from above

% Beam Parameters
E1_4=2*10^11;   E2_5=2*10^11;   E3_6=2*10^11;   % Modulus of Elasticity for each beam in (Pascale)
v_PR = 0.3;   % Poisson’s ratio

% Second moment of area for each beam
b_I=0.2; %Width of the beam in (meter)
h_I=0.2; %Hight of the beam in (meter)

% Multi_Cell plot
Num_hor = 2;        % Number of the unit cell in the X direction
Num_ver = 2;        % Number of the unit cell in the Y direction
Num_up = 2;         % Number of the unit cell in the Z direction

PostProcessing = true;   % Enable (true) or Disable (false)
Plot_Multi_Cell = true;   % Enable (true) or Disable (false)
Elasticity_Tensor = false;   % Enable (true) or Disable (false)

%% Creating geometry
% Calculation of Geometry Parameters
beta=pi/3;  % Angle (1 0 2) in (rad)                           
%gama=pi/3; % Angle (0 1 2) in (rad)
%alpha=pi-(beta+gama);
%a=(sin(gama)*L)/(sin(alpha));
%b=(sin(beta)*L)/(sin(alpha));

angles=[0 beta 2*beta pi pi+beta pi+(2*beta) ];

j=1;
nodes_coord=coord_SNude; % matrix will contain nudes in the rows and their coordinates in the columns
for length=[L L L L L L] 
    for i=0: Num_msh
      node_coord=[coord_SNude(1)+((length/Num_msh)*i*cos(angles(j))) ,coord_SNude(2)+((length/Num_msh)*i*sin(angles(j))),0];  
      nodes_coord=[nodes_coord;node_coord];  
    end 
    j=j+1;
end
nodes_coord(1,:) = []; % Remove the first duplicate row 

lu=(length/Num_msh);
u_nodes_coord=[0, 0, 0];
for ango=[pi/2 (7*pi)/6 (11*pi)/6] 
    for i=0: Num_msh
      node_coord=[coord_SNude(1)+(lu*i*cos(ango)*sqrt(3)/3) ,coord_SNude(2)+(lu*i*sin(ango)*sqrt(3)/3),coord_SNude(3)+((length/Num_msh)*i*sqrt(2/3))];  
      u_nodes_coord=[u_nodes_coord;node_coord];  
    end 
end
u_nodes_coord(1,:) = []; % Remove the first duplicate row 


ld=(length/Num_msh);
d_nodes_coord=[0, 0, 0];
for ango=[pi/6 (5*pi)/6 (3*pi)/2] 
    for i=0: Num_msh
      %node_coord=[coord_SNude(1)+(lu*i*cos(ango)) ,coord_SNude(2)+(lu*i*sin(ango)),coord_SNude(3)+((length/Num_msh)*i*-sin(beta))];  
      node_coord=[coord_SNude(1)+(lu*i*cos(ango)*sqrt(3)/3) ,coord_SNude(2)+(lu*i*sin(ango)*sqrt(3)/3),coord_SNude(3)+((length/Num_msh)*i*-sqrt(2/3))];
      d_nodes_coord=[d_nodes_coord;node_coord];  
    end 
end
d_nodes_coord(1,:) = []; % Remove the first duplicate row 

nodes_coord=[nodes_coord ; u_nodes_coord ; d_nodes_coord];

%% Imposed Boundary Condition
P_1_dis=[ H(1,:)*[nodes_coord(Num_msh+1,1) nodes_coord(Num_msh+1,2) nodes_coord(Num_msh+1,3) ]';H(2,:)*[nodes_coord(Num_msh+1,1) nodes_coord(Num_msh+1,2) nodes_coord(Num_msh+1,3)]';H(3,:)*[nodes_coord(Num_msh+1,1) nodes_coord(Num_msh+1,2) nodes_coord(Num_msh+1,3)]'; 0 ; 0 ; 0];
P_2_dis=[ H(1,:)*[nodes_coord((Num_msh+1)*2,1) nodes_coord((Num_msh+1)*2,2) nodes_coord((Num_msh+1)*2,3) ]';H(2,:)*[nodes_coord((Num_msh+1)*2,1) nodes_coord((Num_msh+1)*2,2) nodes_coord((Num_msh+1)*2,3) ]'; H(3,:)*[nodes_coord((Num_msh+1)*2,1) nodes_coord((Num_msh+1)*2,2) nodes_coord((Num_msh+1)*2,3)]' ; 0 ; 0 ; 0];                            
P_3_dis=[ H(1,:)*[nodes_coord((Num_msh+1)*3,1) nodes_coord((Num_msh+1)*3,2) nodes_coord((Num_msh+1)*3,3) ]';H(2,:)*[nodes_coord((Num_msh+1)*3,1) nodes_coord((Num_msh+1)*3,2) nodes_coord((Num_msh+1)*3,3) ]'; H(3,:)*[nodes_coord((Num_msh+1)*3,1) nodes_coord((Num_msh+1)*3,2) nodes_coord((Num_msh+1)*3,3)]' ; 0 ; 0 ; 0];                            
P_4_dis=[ H(1,:)*[nodes_coord((Num_msh+1)*4,1) nodes_coord((Num_msh+1)*4,2) nodes_coord((Num_msh+1)*4,3) ]';H(2,:)*[nodes_coord((Num_msh+1)*4,1) nodes_coord((Num_msh+1)*4,2) nodes_coord((Num_msh+1)*4,3) ]'; H(3,:)*[nodes_coord((Num_msh+1)*4,1) nodes_coord((Num_msh+1)*4,2) nodes_coord((Num_msh+1)*4,3)]' ; 0 ; 0 ; 0];                            
P_5_dis=[ H(1,:)*[nodes_coord((Num_msh+1)*5,1) nodes_coord((Num_msh+1)*5,2) nodes_coord((Num_msh+1)*5,3) ]';H(2,:)*[nodes_coord((Num_msh+1)*5,1) nodes_coord((Num_msh+1)*5,2) nodes_coord((Num_msh+1)*5,3) ]'; H(3,:)*[nodes_coord((Num_msh+1)*5,1) nodes_coord((Num_msh+1)*5,2) nodes_coord((Num_msh+1)*5,3)]' ; 0 ; 0 ; 0];                            
P_6_dis=[ H(1,:)*[nodes_coord((Num_msh+1)*6,1) nodes_coord((Num_msh+1)*6,2) nodes_coord((Num_msh+1)*6,3) ]';H(2,:)*[nodes_coord((Num_msh+1)*6,1) nodes_coord((Num_msh+1)*6,2) nodes_coord((Num_msh+1)*6,3) ]'; H(3,:)*[nodes_coord((Num_msh+1)*6,1) nodes_coord((Num_msh+1)*6,2) nodes_coord((Num_msh+1)*6,3)]' ; 0 ; 0 ; 0];                            
P_7_dis=[ H(1,:)*[nodes_coord((Num_msh+1)*7,1) nodes_coord((Num_msh+1)*7,2) nodes_coord((Num_msh+1)*7,3) ]';H(2,:)*[nodes_coord((Num_msh+1)*7,1) nodes_coord((Num_msh+1)*7,2) nodes_coord((Num_msh+1)*7,3) ]'; H(3,:)*[nodes_coord((Num_msh+1)*7,1) nodes_coord((Num_msh+1)*7,2) nodes_coord((Num_msh+1)*7,3)]' ; 0 ; 0 ; 0];                            
P_8_dis=[ H(1,:)*[nodes_coord((Num_msh+1)*8,1) nodes_coord((Num_msh+1)*8,2) nodes_coord((Num_msh+1)*8,3) ]';H(2,:)*[nodes_coord((Num_msh+1)*8,1) nodes_coord((Num_msh+1)*8,2) nodes_coord((Num_msh+1)*8,3) ]'; H(3,:)*[nodes_coord((Num_msh+1)*8,1) nodes_coord((Num_msh+1)*8,2) nodes_coord((Num_msh+1)*8,3)]' ; 0 ; 0 ; 0];                            
P_9_dis=[ H(1,:)*[nodes_coord((Num_msh+1)*9,1) nodes_coord((Num_msh+1)*9,2) nodes_coord((Num_msh+1)*9,3) ]';H(2,:)*[nodes_coord((Num_msh+1)*9,1) nodes_coord((Num_msh+1)*9,2) nodes_coord((Num_msh+1)*9,3) ]'; H(3,:)*[nodes_coord((Num_msh+1)*9,1) nodes_coord((Num_msh+1)*9,2) nodes_coord((Num_msh+1)*9,3)]' ; 0 ; 0 ; 0];                            
P_10_dis=[ H(1,:)*[nodes_coord((Num_msh+1)*10,1) nodes_coord((Num_msh+1)*10,2) nodes_coord((Num_msh+1)*10,3) ]';H(2,:)*[nodes_coord((Num_msh+1)*10,1) nodes_coord((Num_msh+1)*10,2) nodes_coord((Num_msh+1)*10,3) ]'; H(3,:)*[nodes_coord((Num_msh+1)*10,1) nodes_coord((Num_msh+1)*10,2) nodes_coord((Num_msh+1)*10,3)]' ; 0 ; 0 ; 0];                            
P_11_dis=[ H(1,:)*[nodes_coord((Num_msh+1)*11,1) nodes_coord((Num_msh+1)*11,2) nodes_coord((Num_msh+1)*11,3) ]';H(2,:)*[nodes_coord((Num_msh+1)*11,1) nodes_coord((Num_msh+1)*11,2) nodes_coord((Num_msh+1)*11,3) ]'; H(3,:)*[nodes_coord((Num_msh+1)*11,1) nodes_coord((Num_msh+1)*11,2) nodes_coord((Num_msh+1)*11,3)]' ; 0 ; 0 ; 0];                            
P_12_dis=[ H(1,:)*[nodes_coord((Num_msh+1)*12,1) nodes_coord((Num_msh+1)*12,2) nodes_coord((Num_msh+1)*12,3) ]';H(2,:)*[nodes_coord((Num_msh+1)*12,1) nodes_coord((Num_msh+1)*12,2) nodes_coord((Num_msh+1)*12,3) ]'; H(3,:)*[nodes_coord((Num_msh+1)*12,1) nodes_coord((Num_msh+1)*12,2) nodes_coord((Num_msh+1)*12,3)]' ; 0 ; 0 ; 0];  

%Imposed Boundary Condition (For elasticity tensor calculation)
H_elast = sym('epo',[3 3]);
P_1_dis_elast=[ H_elast(1,:)*[nodes_coord(Num_msh+1,1) nodes_coord(Num_msh+1,2) nodes_coord(Num_msh+1,3) ]';H_elast(2,:)*[nodes_coord(Num_msh+1,1) nodes_coord(Num_msh+1,2) nodes_coord(Num_msh+1,3)]';H_elast(3,:)*[nodes_coord(Num_msh+1,1) nodes_coord(Num_msh+1,2) nodes_coord(Num_msh+1,3)]'; 0 ; 0 ; 0];
P_2_dis_elast=[ H_elast(1,:)*[nodes_coord((Num_msh+1)*2,1) nodes_coord((Num_msh+1)*2,2) nodes_coord((Num_msh+1)*2,3) ]';H_elast(2,:)*[nodes_coord((Num_msh+1)*2,1) nodes_coord((Num_msh+1)*2,2) nodes_coord((Num_msh+1)*2,3) ]'; H_elast(3,:)*[nodes_coord((Num_msh+1)*2,1) nodes_coord((Num_msh+1)*2,2) nodes_coord((Num_msh+1)*2,3)]' ; 0 ; 0 ; 0];                            
P_3_dis_elast=[ H_elast(1,:)*[nodes_coord((Num_msh+1)*3,1) nodes_coord((Num_msh+1)*3,2) nodes_coord((Num_msh+1)*3,3) ]';H_elast(2,:)*[nodes_coord((Num_msh+1)*3,1) nodes_coord((Num_msh+1)*3,2) nodes_coord((Num_msh+1)*3,3) ]'; H_elast(3,:)*[nodes_coord((Num_msh+1)*3,1) nodes_coord((Num_msh+1)*3,2) nodes_coord((Num_msh+1)*3,3)]' ; 0 ; 0 ; 0];                            
P_4_dis_elast=[ H_elast(1,:)*[nodes_coord((Num_msh+1)*4,1) nodes_coord((Num_msh+1)*4,2) nodes_coord((Num_msh+1)*4,3) ]';H_elast(2,:)*[nodes_coord((Num_msh+1)*4,1) nodes_coord((Num_msh+1)*4,2) nodes_coord((Num_msh+1)*4,3) ]'; H_elast(3,:)*[nodes_coord((Num_msh+1)*4,1) nodes_coord((Num_msh+1)*4,2) nodes_coord((Num_msh+1)*4,3)]' ; 0 ; 0 ; 0];                            
P_5_dis_elast=[ H_elast(1,:)*[nodes_coord((Num_msh+1)*5,1) nodes_coord((Num_msh+1)*5,2) nodes_coord((Num_msh+1)*5,3) ]';H_elast(2,:)*[nodes_coord((Num_msh+1)*5,1) nodes_coord((Num_msh+1)*5,2) nodes_coord((Num_msh+1)*5,3) ]'; H_elast(3,:)*[nodes_coord((Num_msh+1)*5,1) nodes_coord((Num_msh+1)*5,2) nodes_coord((Num_msh+1)*5,3)]' ; 0 ; 0 ; 0];                            
P_6_dis_elast=[ H_elast(1,:)*[nodes_coord((Num_msh+1)*6,1) nodes_coord((Num_msh+1)*6,2) nodes_coord((Num_msh+1)*6,3) ]';H_elast(2,:)*[nodes_coord((Num_msh+1)*6,1) nodes_coord((Num_msh+1)*6,2) nodes_coord((Num_msh+1)*6,3) ]'; H_elast(3,:)*[nodes_coord((Num_msh+1)*6,1) nodes_coord((Num_msh+1)*6,2) nodes_coord((Num_msh+1)*6,3)]' ; 0 ; 0 ; 0];                            
P_7_dis_elast=[ H_elast(1,:)*[nodes_coord((Num_msh+1)*7,1) nodes_coord((Num_msh+1)*7,2) nodes_coord((Num_msh+1)*7,3) ]';H_elast(2,:)*[nodes_coord((Num_msh+1)*7,1) nodes_coord((Num_msh+1)*7,2) nodes_coord((Num_msh+1)*7,3) ]'; H_elast(3,:)*[nodes_coord((Num_msh+1)*7,1) nodes_coord((Num_msh+1)*7,2) nodes_coord((Num_msh+1)*7,3)]' ; 0 ; 0 ; 0];                            
P_8_dis_elast=[ H_elast(1,:)*[nodes_coord((Num_msh+1)*8,1) nodes_coord((Num_msh+1)*8,2) nodes_coord((Num_msh+1)*8,3) ]';H_elast(2,:)*[nodes_coord((Num_msh+1)*8,1) nodes_coord((Num_msh+1)*8,2) nodes_coord((Num_msh+1)*8,3) ]'; H_elast(3,:)*[nodes_coord((Num_msh+1)*8,1) nodes_coord((Num_msh+1)*8,2) nodes_coord((Num_msh+1)*8,3)]' ; 0 ; 0 ; 0];                            
P_9_dis_elast=[ H_elast(1,:)*[nodes_coord((Num_msh+1)*9,1) nodes_coord((Num_msh+1)*9,2) nodes_coord((Num_msh+1)*9,3) ]';H_elast(2,:)*[nodes_coord((Num_msh+1)*9,1) nodes_coord((Num_msh+1)*9,2) nodes_coord((Num_msh+1)*9,3) ]'; H_elast(3,:)*[nodes_coord((Num_msh+1)*9,1) nodes_coord((Num_msh+1)*9,2) nodes_coord((Num_msh+1)*9,3)]' ; 0 ; 0 ; 0];                            
P_10_dis_elast=[ H_elast(1,:)*[nodes_coord((Num_msh+1)*10,1) nodes_coord((Num_msh+1)*10,2) nodes_coord((Num_msh+1)*10,3) ]';H_elast(2,:)*[nodes_coord((Num_msh+1)*10,1) nodes_coord((Num_msh+1)*10,2) nodes_coord((Num_msh+1)*10,3) ]'; H_elast(3,:)*[nodes_coord((Num_msh+1)*10,1) nodes_coord((Num_msh+1)*10,2) nodes_coord((Num_msh+1)*10,3)]' ; 0 ; 0 ; 0];                            
P_11_dis_elast=[ H_elast(1,:)*[nodes_coord((Num_msh+1)*11,1) nodes_coord((Num_msh+1)*11,2) nodes_coord((Num_msh+1)*11,3) ]';H_elast(2,:)*[nodes_coord((Num_msh+1)*11,1) nodes_coord((Num_msh+1)*11,2) nodes_coord((Num_msh+1)*11,3) ]'; H_elast(3,:)*[nodes_coord((Num_msh+1)*11,1) nodes_coord((Num_msh+1)*11,2) nodes_coord((Num_msh+1)*11,3)]' ; 0 ; 0 ; 0];                            
P_12_dis_elast=[ H_elast(1,:)*[nodes_coord((Num_msh+1)*12,1) nodes_coord((Num_msh+1)*12,2) nodes_coord((Num_msh+1)*12,3) ]';H_elast(2,:)*[nodes_coord((Num_msh+1)*12,1) nodes_coord((Num_msh+1)*12,2) nodes_coord((Num_msh+1)*12,3) ]'; H_elast(3,:)*[nodes_coord((Num_msh+1)*12,1) nodes_coord((Num_msh+1)*12,2) nodes_coord((Num_msh+1)*12,3)]' ; 0 ; 0 ; 0];                            



%%  Creating displacement vector to solve (u) & Entry of Boundary conditios

P_BC_imp=[P_1_dis;P_2_dis;P_3_dis;P_4_dis;P_5_dis;P_6_dis;P_7_dis;P_8_dis;P_9_dis;P_10_dis;P_11_dis;P_12_dis];

P_BC_imp_elast=[P_1_dis_elast;P_2_dis_elast;P_3_dis_elast;P_4_dis_elast;P_5_dis_elast;P_6_dis_elast;P_7_dis_elast;P_8_dis_elast;P_9_dis_elast;P_10_dis_elast;P_11_dis_elast;P_12_dis_elast];

%% plotting the geometry and the mesh
scatter3(nodes_coord(:,1),nodes_coord(:,2),nodes_coord(:,3),'DisplayName','Red')
nodes_coord;
hold on
j=0;
for Num_beams=1:12
    for i=0+j: Num_msh-1+j
    line([nodes_coord((i+1),1),nodes_coord((i+2),1)] , [nodes_coord((i+1),2),nodes_coord((i+2),2)] , [nodes_coord((i+1),3),nodes_coord((i+2),3)]) 
    end 
    j=j+Num_msh+1;
end

% Plot the imaginary border of the unit cell

line([L L-(2*L*cos(beta))],[0 (2*L*sin(beta))],[0 0],'Color','green','LineStyle','--')
line([-L L-(2*L*cos(beta))],[0 (2*L*sin(beta))],[0 0],'Color','green','LineStyle','--')
line([-L -L+(2*L*cos((2*pi-beta)))],[0 (2*L*sin((2*pi-beta)))],[0 0],'Color','green','LineStyle','--')
line([L -L+(2*L*cos((2*pi-beta)))],[0 (2*L*sin((2*pi-beta)))],[0 0],'Color','green','LineStyle','--')

line([L+(L*cos((7*pi)/6)*sqrt(3)/3) L-(2*L*cos(beta))+(L*cos((7*pi)/6)*sqrt(3)/3)],[(L*sin((7*pi)/6)*sqrt(3)/3) (2*L*sin(beta))+(L*sin((7*pi)/6)*sqrt(3)/3)],[L*sqrt(2/3) L*sqrt(2/3)],'Color','green','LineStyle','--')
line([-L+(L*cos((7*pi)/6)*sqrt(3)/3) L-(2*L*cos(beta))+(L*cos((7*pi)/6)*sqrt(3)/3)],[(L*sin((7*pi)/6)*sqrt(3)/3) (2*L*sin(beta))+(L*sin((7*pi)/6)*sqrt(3)/3)],[L*sqrt(2/3) L*sqrt(2/3)],'Color','green','LineStyle','--')
line([-L+(L*cos((7*pi)/6)*sqrt(3)/3) -L+(2*L*cos((2*pi-beta)))+(L*cos((7*pi)/6)*sqrt(3)/3)],[(L*sin((7*pi)/6)*sqrt(3)/3) (2*L*sin((2*pi-beta)))+(L*sin((7*pi)/6)*sqrt(3)/3)],[L*sqrt(2/3) L*sqrt(2/3)],'Color','green','LineStyle','--')
line([L+(L*cos((7*pi)/6)*sqrt(3)/3) -L+(2*L*cos((2*pi-beta)))+(L*cos((7*pi)/6)*sqrt(3)/3)],[(L*sin((7*pi)/6)*sqrt(3)/3) (2*L*sin((2*pi-beta)))+(L*sin((7*pi)/6)*sqrt(3)/3)],[L*sqrt(2/3) L*sqrt(2/3)],'Color','green','LineStyle','--')

line([L+(L*cos(pi/6)*sqrt(3)/3) L-(2*L*cos(beta))+(L*cos(pi/6)*sqrt(3)/3)],[(L*sin(pi/6)*sqrt(3)/3) (2*L*sin(beta))+(L*sin(pi/6)*sqrt(3)/3)],[-L*sqrt(2/3) -L*sqrt(2/3)],'Color','green','LineStyle','--')
line([-L+(L*cos(pi/6)*sqrt(3)/3) L-(2*L*cos(beta))+(L*cos(pi/6)*sqrt(3)/3)],[(L*sin(pi/6)*sqrt(3)/3) (2*L*sin(beta))+(L*sin(pi/6)*sqrt(3)/3)],[-L*sqrt(2/3) -L*sqrt(2/3)],'Color','green','LineStyle','--')
line([-L+(L*cos(pi/6)*sqrt(3)/3) -L+(2*L*cos((2*pi-beta)))+(L*cos(pi/6)*sqrt(3)/3)],[(L*sin(pi/6)*sqrt(3)/3) (2*L*sin((2*pi-beta)))+(L*sin(pi/6)*sqrt(3)/3)],[-L*sqrt(2/3) -L*sqrt(2/3)],'Color','green','LineStyle','--')
line([L+(L*cos(pi/6)*sqrt(3)/3) -L+(2*L*cos((2*pi-beta)))+(L*cos(pi/6)*sqrt(3)/3)],[(L*sin(pi/6)*sqrt(3)/3) (2*L*sin((2*pi-beta)))+(L*sin(pi/6)*sqrt(3)/3)],[-L*sqrt(2/3) -L*sqrt(2/3)],'Color','green','LineStyle','--')

line([L+(L*cos((7*pi)/6)*sqrt(3)/3) L+(L*cos(pi/6)*sqrt(3)/3)],[(L*sin((7*pi)/6)*sqrt(3)/3) (L*sin(pi/6)*sqrt(3)/3) ],[L*sqrt(2/3) -L*sqrt(2/3)],'Color','green','LineStyle','--')
line([-L+(L*cos((7*pi)/6)*sqrt(3)/3) -L+(L*cos(pi/6)*sqrt(3)/3) ],[(L*sin((7*pi)/6)*sqrt(3)/3)  (L*sin(pi/6)*sqrt(3)/3) ],[L*sqrt(2/3) -L*sqrt(2/3)],'Color','green','LineStyle','--')
line([-L+(2*L*cos((2*pi-beta)))+(L*cos((7*pi)/6)*sqrt(3)/3)  -L+(2*L*cos((2*pi-beta)))+(L*cos(pi/6)*sqrt(3)/3)],[(2*L*sin((2*pi-beta)))+(L*sin((7*pi)/6)*sqrt(3)/3)  (2*L*sin((2*pi-beta)))+(L*sin(pi/6)*sqrt(3)/3)],[L*sqrt(2/3) -L*sqrt(2/3)],'Color','green','LineStyle','--')
line([L-(2*L*cos(beta))+(L*cos((7*pi)/6)*sqrt(3)/3)  L-(2*L*cos(beta))+(L*cos(pi/6)*sqrt(3)/3)],[(2*L*sin(beta))+(L*sin((7*pi)/6)*sqrt(3)/3)  (2*L*sin(beta))+(L*sin(pi/6)*sqrt(3)/3) ],[L*sqrt(2/3) -L*sqrt(2/3)],'Color','green','LineStyle','--')

axis equal

%% Calculate the elementary Stiffness Matrix for each beam

Ar=b_I*h_I;  % The cross-section of the beam
Iy=(1/12)*b_I*(h_I^3);  % Second Moment of area w.r.t y axis
Iz=(1/12)*h_I*(b_I^3);  % Second Moment of area w.r.t z axis

a_k=b_I/2;c_k=h_I/2; %Parameters of the torsional stiffness factor K
K_TOR= a_k*c_k^3*( (16/3)-(3.36*(c_k/a_k)*(1-(c_k^4/(12*a_k^4)))) );   %The torsional stiffness factor K
J_b=(b_I*h_I*((b_I^2)+(h_I^2)))/12;  % polar moment of inertia 
G1_4=E1_4/(2*(1+v_PR));  G2_5=E2_5/(2*(1+v_PR));   G3_6=E3_6/(2*(1+v_PR)); % modulus ofrigidity

% The parameters for each beam 
E_S=[E1_4 E2_5 E3_6 E1_4 E2_5 E3_6 E1_4 E2_5 E3_6 E1_4 E2_5 E3_6];
G_S=[G1_4 G2_5 G3_6 G1_4 G2_5 G3_6 G1_4 G2_5 G3_6 G1_4 G2_5 G3_6];
Areas=[Ar Ar Ar Ar Ar Ar Ar Ar Ar Ar Ar Ar];
I_Sy=[Iy Iy Iy Iy Iy Iy Iy Iy Iy Iy Iy Iy];
I_Sz=[Iz Iz Iz Iz Iz Iz Iz Iz Iz Iz Iz Iz];
J_S=[J_b J_b J_b J_b J_b J_b J_b J_b J_b J_b J_b J_b];

si=6*(Num_msh+1); % Size of the elementry matrix
d=si-6;  % size of elementry matrix without the common node

j=1;
for OO=2:Num_msh+1:(((Num_msh+1)*12)-Num_msh+1)

% Transformation matrix
l_t=(nodes_coord(OO,1))/(L/Num_msh);
m_t=(nodes_coord(OO,2))/(L/Num_msh);
n_t=(nodes_coord(OO,3))/(L/Num_msh);
D_t=sqrt((l_t^2)+(m_t^2));

Tr=sparse([1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 6 6 6 7 7 7 8 8 8 9 9 9 10 10 10 11 11 11 12 12 12],[1 2 3 1 2 3 ...
    1 2 3 4 5 6 4 5 6 4 5 6 7 8 9 7 8 9 7 8 9 10 11 12 10 11 12 10 11 12],...
    [l_t m_t n_t -m_t/D_t l_t/D_t 0 (-l_t*n_t)/D_t (-m_t*n_t)/D_t D_t l_t m_t n_t -m_t/D_t l_t/D_t 0 (-l_t*n_t)/D_t (-m_t*n_t)/D_t D_t...
    l_t m_t n_t -m_t/D_t l_t/D_t 0 (-l_t*n_t)/D_t (-m_t*n_t)/D_t D_t l_t m_t n_t -m_t/D_t l_t/D_t 0 (-l_t*n_t)/D_t (-m_t*n_t)/D_t D_t],12,12);

% Local stiffness matrix for  element

Ke_s=sparse([1 1 2 2 2 2 3 3 3 3 4 4 5 5 5 5 6 6 6 6 7 7 8 8 8 8 9 9 9 9 10 10 11 11 11 11 12 12 12 12]...
    ,[1 7 2 6 8 12 3 5 9 11 4 10 3 5 9 11 2 6 8 12 1 7 2 6 8 12 3 5 9 11 4 10 3 5 9 11 2 6 8 12],...
    [ Areas(j)*E_S(j)/(length/Num_msh)  -Areas(j)*E_S(j)/(length/Num_msh)  12*E_S(j)*I_Sz(j)/(length/Num_msh)^3   6*E_S(j)*I_Sz(j)/(length/Num_msh)^2   -12*E_S(j)*I_Sz(j)/(length/Num_msh)^3   6*E_S(j)*I_Sz(j)/(length/Num_msh)^2  ...
    12*E_S(j)*I_Sy(j)/(length/Num_msh)^3   -6*E_S(j)*I_Sy(j)/(length/Num_msh)^2   -12*E_S(j)*I_Sy(j)/(length/Num_msh)^3   -6*E_S(j)*I_Sy(j)/(length/Num_msh)^2  ...
    G_S(j)*K_TOR/(length/Num_msh)   -G_S(j)*K_TOR/(length/Num_msh)  -6*E_S(j)*I_Sy(j)/(length/Num_msh)^2   4*E_S(j)*I_Sy(j)/(length/Num_msh)  6*E_S(j)*I_Sy(j)/(length/Num_msh)^2   2*E_S(j)*I_Sy(j)/(length/Num_msh)  ...
    6*E_S(j)*I_Sz(j)/(length/Num_msh)^2   4*E_S(j)*I_Sz(j)/(length/Num_msh)  -6*E_S(j)*I_Sz(j)/(length/Num_msh)^2   2*E_S(j)*I_Sz(j)/(length/Num_msh) ...
    -Areas(j)*E_S(j)/(length/Num_msh)  Areas(j)*E_S(j)/(length/Num_msh)  ...
    -12*E_S(j)*I_Sz(j)/(length/Num_msh)^3  -6*E_S(j)*I_Sz(j)/(length/Num_msh)^2  12*E_S(j)*I_Sz(j)/(length/Num_msh)^3  -6*E_S(j)*I_Sz(j)/(length/Num_msh)^2  ...
    -12*E_S(j)*I_Sy(j)/(length/Num_msh)^3   6*E_S(j)*I_Sy(j)/(length/Num_msh)^2  12*E_S(j)*I_Sy(j)/(length/Num_msh)^3   6*E_S(j)*I_Sy(j)/(length/Num_msh)^2  ...
    -G_S(j)*K_TOR/(length/Num_msh)   G_S(j)*K_TOR/(length/Num_msh)  ...
    -6*E_S(j)*I_Sy(j)/(length/Num_msh)^2  2*E_S(j)*I_Sy(j)/(length/Num_msh)  6*E_S(j)*I_Sy(j)/(length/Num_msh)^2  4*E_S(j)*I_Sy(j)/(length/Num_msh)  ...
     6*E_S(j)*I_Sz(j)/(length/Num_msh)^2  2*E_S(j)*I_Sz(j)/(length/Num_msh)  -6*E_S(j)*I_Sz(j)/(length/Num_msh)^2  4*E_S(j)*I_Sz(j)/(length/Num_msh)],12,12);
    
K_e=Tr' * Ke_s * Tr ;     %Local stiffness matrix after Transformation

% Create the stiffness Matrix for each beam
KK=sparse(6*(Num_msh+1),6*(Num_msh+1));
for i=0:12:(12*(Num_msh-1))
KK((i/2)+1:(i/2)+12,(i/2)+1:(i/2)+12)=KK((i/2)+1:(i/2)+12,(i/2)+1:(i/2)+12)+K_e;
end

   if j == 1
       K1 = KK;
   elseif j == 2
       K2 = KK;
   elseif j == 3
       K3 = KK;
   elseif j == 4
       K4 = KK;
   elseif j == 5
       K5 = KK;
   elseif j == 6
       K6 = KK;
   elseif j == 7
       K7 = KK;
   elseif j == 8
       K8 = KK;
   elseif j == 9
       K9 = KK;
   elseif j == 10
       K10 = KK;
   elseif j == 11
       K11 = KK;
   else
       K12 = KK;
   end 

j=j+1;
end


%% Assemble of the global stiffness Matrix for the whole unit-cell

K_glob=sparse((78*Num_msh)-(6*(Num_msh-1)),(78*Num_msh)-(6*(Num_msh-1)));
K_glob(1:si,1:si)=K1;
K_glob(1:6,1:6)=K1(1:6,1:6)+K2(1:6,1:6)+K3(1:6,1:6)+K4(1:6,1:6)+K5(1:6,1:6)+K6(1:6,1:6)+K7(1:6,1:6)+K8(1:6,1:6)+K9(1:6,1:6)+K10(1:6,1:6)+K11(1:6,1:6)+K12(1:6,1:6);

K_glob(si+1:si+d,si+1:si+d)=K2(7:end,7:end);
K_glob(1:6,si+1:si+d)=K2(1:6,7:end);
K_glob(si+1:si+d,1:6)=K2(7:end,1:6);

K_glob(si+(1*d)+1:si+(2*d),si+(1*d)+1:si+(2*d))=K3(7:end,7:end);
K_glob(1:6,si+(1*d)+1:si+(2*d))=K3(1:6,7:end);
K_glob(si+(1*d)+1:si+(2*d),1:6)=K3(7:end,1:6);

K_glob(si+(2*d)+1:si+(3*d),si+(2*d)+1:si+(3*d))=K4(7:end,7:end);
K_glob(1:6,si+(2*d)+1:si+(3*d))=K4(1:6,7:end);
K_glob(si+(2*d)+1:si+(3*d),1:6)=K4(7:end,1:6);

K_glob(si+(3*d)+1:si+(4*d),si+(3*d)+1:si+(4*d))=K5(7:end,7:end);
K_glob(1:6,si+(3*d)+1:si+(4*d))=K5(1:6,7:end);
K_glob(si+(3*d)+1:si+(4*d),1:6)=K5(7:end,1:6);

K_glob(si+(4*d)+1:si+(5*d),si+(4*d)+1:si+(5*d))=K6(7:end,7:end);
K_glob(1:6,si+(4*d)+1:si+(5*d))=K6(1:6,7:end);
K_glob(si+(4*d)+1:si+(5*d),1:6)=K6(7:end,1:6);

K_glob(si+(5*d)+1:si+(6*d),si+(5*d)+1:si+(6*d))=K7(7:end,7:end);
K_glob(1:6,si+(5*d)+1:si+(6*d))=K7(1:6,7:end);
K_glob(si+(5*d)+1:si+(6*d),1:6)=K7(7:end,1:6);

K_glob(si+(6*d)+1:si+(7*d),si+(6*d)+1:si+(7*d))=K8(7:end,7:end);
K_glob(1:6,si+(6*d)+1:si+(7*d))=K8(1:6,7:end);
K_glob(si+(6*d)+1:si+(7*d),1:6)=K8(7:end,1:6);

K_glob(si+(7*d)+1:si+(8*d),si+(7*d)+1:si+(8*d))=K9(7:end,7:end);
K_glob(1:6,si+(7*d)+1:si+(8*d))=K9(1:6,7:end);
K_glob(si+(7*d)+1:si+(8*d),1:6)=K9(7:end,1:6);

K_glob(si+(8*d)+1:si+(9*d),si+(8*d)+1:si+(9*d))=K10(7:end,7:end);
K_glob(1:6,si+(8*d)+1:si+(9*d))=K10(1:6,7:end);
K_glob(si+(8*d)+1:si+(9*d),1:6)=K10(7:end,1:6);

K_glob(si+(9*d)+1:si+(10*d),si+(9*d)+1:si+(10*d))=K11(7:end,7:end);
K_glob(1:6,si+(9*d)+1:si+(10*d))=K11(1:6,7:end);
K_glob(si+(9*d)+1:si+(10*d),1:6)=K11(7:end,1:6);

K_glob(si+(10*d)+1:si+(11*d),si+(10*d)+1:si+(11*d))=K12(7:end,7:end);
K_glob(1:6,si+(10*d)+1:si+(11*d))=K12(1:6,7:end);
K_glob(si+(10*d)+1:si+(11*d),1:6)=K12(7:end,1:6);


%% Calculating the periodic field

j=1;Ke_s=0;
for OO=2:Num_msh+1:(((Num_msh+1)*12)-Num_msh+1)

% Transformation matrix
l_t=(nodes_coord(OO,1))/(L/Num_msh);
m_t=(nodes_coord(OO,2))/(L/Num_msh);
n_t=(nodes_coord(OO,3))/(L/Num_msh);
D_t=sqrt((l_t^2)+(m_t^2));

Tr=sparse([1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 6 6 6 7 7 7 8 8 8 9 9 9 10 10 10 11 11 11 12 12 12],[1 2 3 1 2 3 ...
    1 2 3 4 5 6 4 5 6 4 5 6 7 8 9 7 8 9 7 8 9 10 11 12 10 11 12 10 11 12],...
    [l_t m_t n_t -m_t/D_t l_t/D_t 0 (-l_t*n_t)/D_t (-m_t*n_t)/D_t D_t l_t m_t n_t -m_t/D_t l_t/D_t 0 (-l_t*n_t)/D_t (-m_t*n_t)/D_t D_t...
    l_t m_t n_t -m_t/D_t l_t/D_t 0 (-l_t*n_t)/D_t (-m_t*n_t)/D_t D_t l_t m_t n_t -m_t/D_t l_t/D_t 0 (-l_t*n_t)/D_t (-m_t*n_t)/D_t D_t],12,12);

% Local stiffness matrix for  element
Ke_s=sparse([1 1 2 2 2 2 3 3 3 3 4 4 5 5 5 5 6 6 6 6 7 7 8 8 8 8 9 9 9 9 10 10 11 11 11 11 12 12 12 12]...
    ,[1 7 2 6 8 12 3 5 9 11 4 10 3 5 9 11 2 6 8 12 1 7 2 6 8 12 3 5 9 11 4 10 3 5 9 11 2 6 8 12],...
    [ Areas(j)*E_S(j)/(length)  -Areas(j)*E_S(j)/(length)  12*E_S(j)*I_Sz(j)/(length)^3   6*E_S(j)*I_Sz(j)/(length)^2   -12*E_S(j)*I_Sz(j)/(length)^3   6*E_S(j)*I_Sz(j)/(length)^2  ...
    12*E_S(j)*I_Sy(j)/(length)^3   -6*E_S(j)*I_Sy(j)/(length)^2   -12*E_S(j)*I_Sy(j)/(length)^3   -6*E_S(j)*I_Sy(j)/(length)^2  ...
    G_S(j)*K_TOR/(length)   -G_S(j)*K_TOR/(length)  -6*E_S(j)*I_Sy(j)/(length)^2   4*E_S(j)*I_Sy(j)/(length)  6*E_S(j)*I_Sy(j)/(length)^2   2*E_S(j)*I_Sy(j)/(length)  ...
    6*E_S(j)*I_Sz(j)/(length)^2   4*E_S(j)*I_Sz(j)/(length)  -6*E_S(j)*I_Sz(j)/(length)^2   2*E_S(j)*I_Sz(j)/(length) ...
    -Areas(j)*E_S(j)/(length)  Areas(j)*E_S(j)/(length)  ...
    -12*E_S(j)*I_Sz(j)/(length)^3  -6*E_S(j)*I_Sz(j)/(length)^2  12*E_S(j)*I_Sz(j)/(length)^3  -6*E_S(j)*I_Sz(j)/(length)^2  ...
    -12*E_S(j)*I_Sy(j)/(length)^3   6*E_S(j)*I_Sy(j)/(length)^2  12*E_S(j)*I_Sy(j)/(length)^3   6*E_S(j)*I_Sy(j)/(length)^2  ...
    -G_S(j)*K_TOR/(length)   G_S(j)*K_TOR/(length)  ...
    -6*E_S(j)*I_Sy(j)/(length)^2  2*E_S(j)*I_Sy(j)/(length)  6*E_S(j)*I_Sy(j)/(length)^2  4*E_S(j)*I_Sy(j)/(length)  ...
     6*E_S(j)*I_Sz(j)/(length)^2  2*E_S(j)*I_Sz(j)/(length)  -6*E_S(j)*I_Sz(j)/(length)^2  4*E_S(j)*I_Sz(j)/(length)],12,12);

K_e=Tr' * Ke_s * Tr ;     %Local stiffness matrix after Transformation

% Create the stiffness Matrix for each beam
   if j == 1
       K1_per = K_e;
   elseif j == 2
       K2_per = K_e;
   elseif j == 3
       K3_per = K_e;
   elseif j == 4
       K4_per = K_e;
   elseif j == 5
       K5_per = K_e;
   elseif j == 6
       K6_per = K_e;
   elseif j == 7
       K7_per = K_e;
   elseif j == 8
       K8_per = K_e;
   elseif j == 9
       K9_per = K_e;
   elseif j == 10
       K10_per = K_e;
   elseif j == 11
       K11_per = K_e;
   else
       K12_per = K_e;
   end 

j=j+1;
end

si_per=12; % Size of the elementry matrix
d_per=6;  % size of elementry matrix without the common node
%asemble the stiffness for calculating the periodic displacement field 
K_per=sparse(78,78);
K_per(1:si_per,1:si_per)=K1_per;
K_per(1:6,1:6)=K1_per(1:6,1:6)+K2_per(1:6,1:6)+K3_per(1:6,1:6)+K4_per(1:6,1:6)+K5_per(1:6,1:6)+K6_per(1:6,1:6)+K7_per(1:6,1:6)+K8_per(1:6,1:6)+K9_per(1:6,1:6)+K10_per(1:6,1:6)+K11_per(1:6,1:6)+K12_per(1:6,1:6);

K_per(si_per+1:si_per+d_per,si_per+1:si_per+d_per)=K2_per(7:end,7:end);
K_per(1:6,si_per+1:si_per+d_per)=K2_per(1:6,7:end);
K_per(si_per+1:si_per+d_per,1:6)=K2_per(7:end,1:6);

K_per(si_per+(1*d_per)+1:si_per+(2*d_per),si_per+(1*d_per)+1:si_per+(2*d_per))=K3_per(7:end,7:end);
K_per(1:6,si_per+(1*d_per)+1:si_per+(2*d_per))=K3_per(1:6,7:end);
K_per(si_per+(1*d_per)+1:si_per+(2*d_per),1:6)=K3_per(7:end,1:6);


K_per(si_per+(2*d_per)+1:si_per+(3*d_per),si_per+(2*d_per)+1:si_per+(3*d_per))=K4_per(7:end,7:end);
K_per(1:6,si_per+(2*d_per)+1:si_per+(3*d_per))=K4_per(1:6,7:end);
K_per(si_per+(2*d_per)+1:si_per+(3*d_per),1:6)=K4_per(7:end,1:6);

K_per(si_per+(3*d_per)+1:si_per+(4*d_per),si_per+(3*d_per)+1:si_per+(4*d_per))=K5_per(7:end,7:end);
K_per(1:6,si_per+(3*d_per)+1:si_per+(4*d_per))=K5_per(1:6,7:end);
K_per(si_per+(3*d_per)+1:si_per+(4*d_per),1:6)=K5_per(7:end,1:6);

K_per(si_per+(4*d_per)+1:si_per+(5*d_per),si_per+(4*d_per)+1:si_per+(5*d_per))=K6_per(7:end,7:end);
K_per(1:6,si_per+(4*d_per)+1:si_per+(5*d_per))=K6_per(1:6,7:end);
K_per(si_per+(4*d_per)+1:si_per+(5*d_per),1:6)=K6_per(7:end,1:6);

K_per(si_per+(5*d_per)+1:si_per+(6*d_per),si_per+(5*d_per)+1:si_per+(6*d_per))=K7_per(7:end,7:end);
K_per(1:6,si_per+(5*d_per)+1:si_per+(6*d_per))=K7_per(1:6,7:end);
K_per(si_per+(5*d_per)+1:si_per+(6*d_per),1:6)=K7_per(7:end,1:6);

K_per(si_per+(6*d_per)+1:si_per+(7*d_per),si_per+(6*d_per)+1:si_per+(7*d_per))=K8_per(7:end,7:end);
K_per(1:6,si_per+(6*d_per)+1:si_per+(7*d_per))=K8_per(1:6,7:end);
K_per(si_per+(6*d_per)+1:si_per+(7*d_per),1:6)=K8_per(7:end,1:6);

K_per(si_per+(7*d_per)+1:si_per+(8*d_per),si_per+(7*d_per)+1:si_per+(8*d_per))=K9_per(7:end,7:end);
K_per(1:6,si_per+(7*d_per)+1:si_per+(8*d_per))=K9_per(1:6,7:end);
K_per(si_per+(7*d_per)+1:si_per+(8*d_per),1:6)=K9_per(7:end,1:6);

K_per(si_per+(8*d_per)+1:si_per+(9*d_per),si_per+(8*d_per)+1:si_per+(9*d_per))=K10_per(7:end,7:end);
K_per(1:6,si_per+(8*d_per)+1:si_per+(9*d_per))=K10_per(1:6,7:end);
K_per(si_per+(8*d_per)+1:si_per+(9*d_per),1:6)=K10_per(7:end,1:6);

K_per(si_per+(9*d_per)+1:si_per+(10*d_per),si_per+(9*d_per)+1:si_per+(10*d_per))=K11_per(7:end,7:end);
K_per(1:6,si_per+(9*d_per)+1:si_per+(10*d_per))=K11_per(1:6,7:end);
K_per(si_per+(9*d_per)+1:si_per+(10*d_per),1:6)=K11_per(7:end,1:6);

K_per(si_per+(10*d_per)+1:si_per+(11*d_per),si_per+(10*d_per)+1:si_per+(11*d_per))=K12_per(7:end,7:end);
K_per(1:6,si_per+(10*d_per)+1:si_per+(11*d_per))=K12_per(1:6,7:end);
K_per(si_per+(10*d_per)+1:si_per+(11*d_per),1:6)=K12_per(7:end,1:6);


IU=eye(24);
IM=eye(36);
III=eye(6);
ZZZ=0*eye(6);
IL=[ZZZ III ZZZ ; ZZZ ZZZ III ; III ZZZ ZZZ];

Z0=zeros(78,42);
Z0(1:24,1:24)=IU;
Z0(25:60,7:42)=IM;
Z0(61:78,25:42)=IL;


%Addind the aditional coditions of the nulity of the average displacment.
I_sol=[1 0 0 0 0 0 ; 0 1 0 0 0 0 ; 0 0 1 0 0 0];
c=[I_sol I_sol I_sol I_sol I_sol I_sol I_sol];
hh=Z0' * K_per *Z0;


qq=[0;0;0;0;0;0;P_BC_imp];
qq_elast=[0;0;0;0;0;0;P_BC_imp_elast];

ff= -Z0' * K_per * qq;
ff_elast= -Z0' * K_per * qq_elast;
% using lagrange multiplyer
zero=zeros(3,3);
tt=[hh c';c zero];
fff=[ff;0;0;0];
fff_elast=[ff_elast;0;0;0];

XX = tt\fff;
XX_elast = tt\fff_elast;

XX=XX(1:42);
XX_elast=XX_elast(1:42);

P_BC_TOTAL=Z0*XX+qq;
P_BC_TOTAL_elast=Z0*XX_elast+qq_elast;

P_BC=P_BC_TOTAL(7:end);

%% Create Reduced system to solve because there are prescriped values.
  k=K_glob;

if Num_msh == 1
    u_red=P_BC_TOTAL;
   f=zeros(78,1);

   for i=1:72

   f=(u_red(end,1)*k(:,end))+f;
   k(end,:)=[];
   k(:,end)=[];
   u_red(end)=[];
   f(end)=[];
   end
   X_red = k\-f;
      X=X_red;
 
else 
    
 % Create the unkown victor (X) -->  A X = B
 % Entry the prescribed value (the boundary conditions)
 % For the unknown values entry (123.6) very especial value to acurately determination 
 % of that elemnt in the next step which is elemnate the known values
    hh=1;
    u_red=ones(size(K_glob,2),1)*123.6;   
    for i=4+d-3:d:size(K_glob,2)
    u_red(i)=P_BC(hh);
    u_red(i+1)=P_BC(hh+1);
    u_red(i+2)=P_BC(hh+2);
    u_red(i+3)=P_BC(hh+3);
    u_red(i+4)=P_BC(hh+4);
    u_red(i+5)=P_BC(hh+5);
    hh=hh+6;
    end
    u_red;
    
  
    k=K_glob;
    f=zeros(size(K_glob,2),1);
    for i=1:size(K_glob,2)
    if u_red(i)~=123.6
    fo=K_glob(:,i)*u_red(i);
    f=f+fo;
    k(:,i)=zeros(size(K_glob,2),1);
    k(i,:)=zeros(size(K_glob,2),1);
    k(i,i)=1;
    end
        
    end


X_red = k \ -f;

% eleminate the prescribed values
X_red(si-5:si)=[];
X_red(((d+1)+(d-6)):((d+6)+(d-6)))=[];
X_red(((d+1)+(2*(d-6))):(d+6)+(2*(d-6)))=[];
X_red(((d+1)+(3*(d-6))):(d+6)+(3*(d-6)))=[];
X_red(((d+1)+(4*(d-6))):(d+6)+(4*(d-6)))=[];
X_red(((d+1)+(5*(d-6))):(d+6)+(5*(d-6)))=[];

X_red(((d+1)+(6*(d-6))):(d+6)+(6*(d-6)))=[];
X_red(((d+1)+(7*(d-6))):(d+6)+(7*(d-6)))=[];
X_red(((d+1)+(8*(d-6))):(d+6)+(8*(d-6)))=[];
X_red(((d+1)+(9*(d-6))):(d+6)+(9*(d-6)))=[];
X_red(((d+1)+(10*(d-6))):(d+6)+(10*(d-6)))=[];
X_red(((d+1)+(11*(d-6))):(d+6)+(11*(d-6)))=[];
end
X=X_red;
t=1; h=7; dd=2;
% nodes_coord_after --> each row contain (three direction + three slope)
% for each node
for z=1:12    
for i=2:Num_msh
nodes_coord_after(dd,1)=nodes_coord(dd,1)+X(h);
nodes_coord_after(dd,2)=nodes_coord(dd,2)+X(h+1);
nodes_coord_after(dd,3)=nodes_coord(dd,3)+X(h+2);

nodes_coord_after(dd,4)=X(h+3);
nodes_coord_after(dd,5)=X(h+4);
nodes_coord_after(dd,6)=X(h+5);

h=h+6;dd=dd+1;
end
nodes_coord_after(dd,1)=nodes_coord(dd,1)+P_BC(t);
nodes_coord_after(dd,2)=nodes_coord(dd,2)+P_BC(t+1);
nodes_coord_after(dd,3)=nodes_coord(dd,3)+P_BC(t+2);
nodes_coord_after(dd,4)=P_BC(t+3);
nodes_coord_after(dd,5)=P_BC(t+4);
nodes_coord_after(dd,6)=P_BC(t+5);
dd=dd+2;t=t+6;
end
hold on
s=1;
for i=1:12   % delete dublicate values
    nodes_coord_after(s,:)=[];
    s=s+Num_msh;
end
nodes_coord_after=[X(1)+XX(1) X(2)+XX(2) X(3)+XX(3) X(4) X(5) X(6);nodes_coord_after];

hold on
% plot the solved node
scatter3(nodes_coord_after(:,1),nodes_coord_after(:,2),nodes_coord_after(:,3),'red');

% reform the arrays of node to use it in postprocessin
node_after=zeros(size(nodes_coord,1),6);

co=1; ci=2;
for i=1:12 
    node_after(co,:)= nodes_coord_after(1,:);
    co=co+1;
    for j=1: Num_msh
      node_after(co,:)=nodes_coord_after(ci,:);  
      co=co+1;ci=ci+1;  
    end 
end

%% Plot Interpolation line
G=2; j=1;
for OO=2:Num_msh+1:(((Num_msh+1)*12)-Num_msh+1)
    % Transformation matrix
l_t=(nodes_coord(OO,1))/(L/Num_msh);
m_t=(nodes_coord(OO,2))/(L/Num_msh);
n_t=(nodes_coord(OO,3))/(L/Num_msh);
D_t=sqrt((l_t^2)+(m_t^2));

Tr=sparse([1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 6 6 6 7 7 7 8 8 8 9 9 9 10 10 10 11 11 11 12 12 12],[1 2 3 1 2 3 ...
    1 2 3 4 5 6 4 5 6 4 5 6 7 8 9 7 8 9 7 8 9 10 11 12 10 11 12 10 11 12],...
    [l_t m_t n_t -m_t/D_t l_t/D_t 0 (-l_t*n_t)/D_t (-m_t*n_t)/D_t D_t l_t m_t n_t -m_t/D_t l_t/D_t 0 (-l_t*n_t)/D_t (-m_t*n_t)/D_t D_t...
    l_t m_t n_t -m_t/D_t l_t/D_t 0 (-l_t*n_t)/D_t (-m_t*n_t)/D_t D_t l_t m_t n_t -m_t/D_t l_t/D_t 0 (-l_t*n_t)/D_t (-m_t*n_t)/D_t D_t],12,12);
Trr=Tr(1:6,1:6);
 
u=100;
v=100;
w=100;
u_general=0;
v_general=0;
w_general=0;
l=length/Num_msh;
for j=2:Num_msh+1;
 
p_inter=3;    
for x=linspace(0,l,p_inter)
    ze= x/l;
q1_glo=[node_after(G-1,1);node_after(G-1,2);node_after(G-1,3);node_after(G-1,4);node_after(G-1,5);node_after(G-1,6)];
q2_glo=[node_after(G,1);node_after(G,2);node_after(G,3);node_after(G,4);node_after(G,5);node_after(G,6)];


q1_loc=Trr * q1_glo ;
q2_loc=Trr * q2_glo ;

u1_m=q1_loc(1);v1_m=q1_loc(2);w1_m=q1_loc(3);m1_m1=q1_loc(4);m1_m2=q1_loc(5);m1_m3=q1_loc(6);
u2_m=q2_loc(1);v2_m=q2_loc(2);w2_m=q2_loc(3);m2_m1=q2_loc(4);m2_m2=q2_loc(5);m2_m3=q2_loc(6);


u_general11=(1-ze)*u1_m + (ze*u2_m);
v_general11= (1-(3*(ze^2))+(2*(ze^3)))*v1_m +(ze*((1-ze)^2)*l*m1_m2) + (((3*(ze^2))-(2*(ze^3)))*v2_m) +((ze^2)*(ze-1)*l*m2_m2)   ;
w_general11= (1-(3*(ze^2))+(2*(ze^3)))*w1_m +(ze*((1-ze)^2)*l*m1_m2) + (((3*(ze^2))-(2*(ze^3)))*w2_m) +((ze^2)*(ze-1)*l*m2_m2)   ;

new_q=Trr'*[u_general11;v_general11;w_general11;0;0;0];

u_general=new_q(1); v_general=new_q(2); w_general=new_q(3);

u=[u,u_general];
v=[v,v_general];
w=[w,w_general];

end
G=G+1;
end
G=G+1;
u(1)=[];
v(1)=[];
w(1)=[];
intpolation=plot3(u,v,w,'c','LineWidth',2);
hold on
end
hold on
title('After solving') 
hhhh = zeros(3, 1);
hhhh(1) = plot(NaN,NaN,'or');
hhhh(2) = plot(NaN,NaN,'c');
hhhh(3) = plot(NaN,NaN,'ob');
hhhh(4) = plot(NaN,NaN,'--g');

xlabel('X_{Global} ') 
ylabel('Y_{Global} ') 
zlabel('Z_{Global} ') 
legend(hhhh, 'Solved Nodes','Interpolation line ','Before solving', 'Cell Border-Before' );
hold off

%% Plot Multi-Cells

if Plot_Multi_Cell== true
figure
axis equal
title('Multi-Cell')
xlabel('x Axis') 
ylabel('y Axis') 
zlabel('z Axis') 
sh=0;shh=0;
sv=0;svv=0;
sz=0;szz=0;

for up=1:Num_up
sh=0;
sv=0;
sz=0;

for vv=1:Num_ver
for h=0:Num_hor-1
    hold on   
    scatter3(nodes_coord_after(:,1)+(((2*L)+(P_BC(1)-P_BC(19)))*h)+sh+shh  ,  nodes_coord_after(:,2)+((P_BC(2)-P_BC(20))*h)+sv+svv, nodes_coord_after(:,3)+((P_BC(3)-P_BC(21))*h)+sz+szz,'red'); 
G=2;
for length=[L L L L L L L L L L L L]
u=100;
v=100;
w=100;
u_general=0;
v_general=0;
w_general=0;
l=length/Num_msh;
for j=2:Num_msh+1;
    i=1;
for x=0:l/2:l
    ze= x/l;
u_general(i)=(1-ze)*node_after(G-1,1) + (ze*node_after(G,1)) + (((2*L)+(P_BC(1)-P_BC(19)))*h)+sh+shh ;
v_general(i)=((1-(3*ze^2)+(2*ze^3))*(node_after(G-1,2)))+( l * ze * ((1-ze)^2) * node_after(G-1,4) ) + (( (3*ze^2)-(2*ze^3) )*node_after(G,2))+ ( l*ze*ze*(ze-1)*node_after(G,4) )+   ((P_BC(2)-P_BC(20))*h)+sv+svv;
w_general(i)=((1-(3*ze^2)+(2*ze^3))*(node_after(G-1,3)))+( l * ze * ((1-ze)^2) * node_after(G-1,5) ) + (( (3*ze^2)-(2*ze^3) )*node_after(G,3))+ ( l*ze*ze*(ze-1)*node_after(G,5) )+((P_BC(3)-P_BC(21))*h)+sz+szz;

i=i+1;
end
u=[u,u_general];
v=[v,v_general];
w=[w,w_general];

G=G+1;

end
G=G+1;
u(1)=[];
v(1)=[];
w(1)=[];
intpolation=plot3(u,v,w,'c','LineWidth',2);
hold on
end
    
end
sh=((2*L*cos(beta))+((P_BC(7))-(P_BC(25)))) * vv;
sv=((2*L*sin(beta))+((P_BC(8))-(P_BC(26)))) * vv;
sz=((P_BC(9))-(P_BC(27)))* vv ;
end

shh=((2*L*cos(pi/6)*sqrt(3)/3)+(2*P_BC(55))) * up ;
svv=((2*L*sin(pi/6)*sqrt(3)/3)+(2*P_BC(56))) * up ;
szz=((-2*L*sqrt(2/3))+(2*(P_BC(57)))) * up ;

end
axis equal
end

%% Post-Processing
if PostProcessing== true
   
% Assemble total displacement node vector 
dis_nod=1234.6;
cen_nod=X(1:6);
cv=7;
for j=1:12
    solve=X(cv:cv+(d-7));
    boundary=P_BC(1+((j-1)*6):6+((j-1)*6));
    NEW=[cen_nod ; solve ; boundary];
    dis_nod=[dis_nod ; NEW];
    cv=(cv+(d-7))+1;
end
dis_nod(1)=[];

%calculate and plot the normal force
figure;
p.a = 'r' ;p.b = 'c'; p.c = 'k';p.d = 'or';p.e = 'oc';p.f = 'ok'; %structre for color plotting
p.g = 'sm';p.h = 'sg';p.i = 'sb';p.j = '^m';p.k = '^g';p.l = '^b';
y = struct2cell(p);
ys=[10 10 10 5 7 11 9 11 12 8 10 12];
G=1; o=1; rol=1; ang=1; j=1;
OO_tr=2:Num_msh+1:(((Num_msh+1)*12)-Num_msh+1);
for length=[L L L L L L L L L L L L]
i=1;
l=length/Num_msh;
for ty=1:Num_msh;
% Transformation matrix
l_t=(nodes_coord(OO_tr(j),1))/(L/Num_msh);
m_t=(nodes_coord(OO_tr(j),2))/(L/Num_msh);
n_t=(nodes_coord(OO_tr(j),3))/(L/Num_msh);
D_t=sqrt((l_t^2)+(m_t^2));

Lmbda_3_3=[l_t m_t n_t ; -m_t/D_t l_t/D_t 0 ; (-l_t*n_t)/D_t (-m_t*n_t)/D_t D_t ];
I_t=zeros(3,3);
Tr=[Lmbda_3_3 I_t I_t I_t ; I_t Lmbda_3_3 I_t I_t ; I_t I_t Lmbda_3_3 I_t ; I_t I_t I_t Lmbda_3_3 ] ;

d_dis_shear=[dis_nod(G); dis_nod(G+1) ; dis_nod(G+2);dis_nod(G+3) ;dis_nod(G+4) ;dis_nod(G+5) ...
    ;dis_nod(G+6) ;dis_nod(G+7) ;dis_nod(G+8) ;dis_nod(G+9) ;dis_nod(G+10) ;dis_nod(G+11)] ;   

hof=Tr*d_dis_shear;

ss=hof(1);
ss2=hof(7);
    
N_general(i)=((((-E_S(rol)*Areas(rol))/l)*ss) + (((E_S(rol)*Areas(rol))/l)*ss2));
G=G+6;i=i+1;
cos(angles(ang));
end
G=G+6; j=j+1;
xxx = linspace(0,length,size(N_general,2));

plot(xxx,N_general,y{o},'MarkerSize',ys(o))
o=o+1; rol=rol+1;
title('Normal force')
xlabel('Length of beam (m)')
ylabel('Normal force (N)')
hold on
end
legend('Beam [1]','Beam [2]','Beam [3]','Beam [4]','Beam [5]','Beam [6]','Beam [7]','Beam [8]','Beam [9]','Beam [10]','Beam [11]','Beam [12]')
hold off

%calculate and plot the transverse force in y direction
figure;
o=1; G=1; rol=1; ang=1; j=1;
OO_tr=2:Num_msh+1:(((Num_msh+1)*12)-Num_msh+1);
for length=[L L L L L L L L L L L L]
i=1;
l=length/Num_msh;
for ty=1:Num_msh;
% Transformation matrix
l_t=(nodes_coord(OO_tr(j),1))/(L/Num_msh);
m_t=(nodes_coord(OO_tr(j),2))/(L/Num_msh);
n_t=(nodes_coord(OO_tr(j),3))/(L/Num_msh);
D_t=sqrt((l_t^2)+(m_t^2));

Lmbda_3_3=[l_t m_t n_t ; -m_t/D_t l_t/D_t 0 ; (-l_t*n_t)/D_t (-m_t*n_t)/D_t D_t ];
I_t=zeros(3,3);
Tr=[Lmbda_3_3 I_t I_t I_t ; I_t Lmbda_3_3 I_t I_t ; I_t I_t Lmbda_3_3 I_t ; I_t I_t I_t Lmbda_3_3 ] ;

d_dis_shear=[dis_nod(G); dis_nod(G+1) ; dis_nod(G+2);dis_nod(G+3) ;dis_nod(G+4) ;dis_nod(G+5) ...
    ;dis_nod(G+6) ;dis_nod(G+7) ;dis_nod(G+8) ;dis_nod(G+9) ;dis_nod(G+10) ;dis_nod(G+11)] ;   

hof=Tr*d_dis_shear;
vs=hof(2);
vs2=hof(8);
ms=hof(6);
ms2=hof(12);
 
Ty_general(i)= ((-12*E_S(rol)*I_Sz(rol)/(l^3))*vs) - ((6*E_S(rol)*I_Sz(rol)/(l^2))*ms) + ((12*E_S(rol)*I_Sz(rol)/(l^3))*vs2)-((6*E_S(rol)*I_Sz(rol)/(l^2))*ms2);

G=G+6;i=i+1;
cos(angles(ang));
end
G=G+6; j=j+1;
xxx = linspace(0,length,size(Ty_general,2));

plot(xxx,Ty_general,y{o},'MarkerSize',ys(o))
o=o+1; rol=rol+1;
title('Shear force in Y direction')
xlabel('Length of beam (m)')
ylabel('shear force (N)')
hold on
end
legend('Beam [1]','Beam [2]','Beam [3]','Beam [4]','Beam [5]','Beam [6]','Beam [7]','Beam [8]','Beam [9]','Beam [10]','Beam [11]','Beam [12]')
hold off

%calculate and plot the transverse force in z direction
figure;
o=1; G=1; rol=1; ang=1; j=1;
OO_tr=2:Num_msh+1:(((Num_msh+1)*12)-Num_msh+1);
for length=[L L L L L L L L L L L L]
i=1;
l=length/Num_msh;
for ty=1:Num_msh;
% Transformation matrix
l_t=(nodes_coord(OO_tr(j),1))/(L/Num_msh);
m_t=(nodes_coord(OO_tr(j),2))/(L/Num_msh);
n_t=(nodes_coord(OO_tr(j),3))/(L/Num_msh);
D_t=sqrt((l_t^2)+(m_t^2));

Lmbda_3_3=[l_t m_t n_t ; -m_t/D_t l_t/D_t 0 ; (-l_t*n_t)/D_t (-m_t*n_t)/D_t D_t ];
I_t=zeros(3,3);
Tr=[Lmbda_3_3 I_t I_t I_t ; I_t Lmbda_3_3 I_t I_t ; I_t I_t Lmbda_3_3 I_t ; I_t I_t I_t Lmbda_3_3 ] ;

d_dis_shear=[dis_nod(G); dis_nod(G+1) ; dis_nod(G+2);dis_nod(G+3) ;dis_nod(G+4) ;dis_nod(G+5) ...
    ;dis_nod(G+6) ;dis_nod(G+7) ;dis_nod(G+8) ;dis_nod(G+9) ;dis_nod(G+10) ;dis_nod(G+11)] ;   

hofz=Tr*d_dis_shear;
vs=hofz(3);
vs2=hofz(9);
ms=hofz(5);
ms2=hofz(11);
 
Tz_general(i)= ((-12*E_S(rol)*I_Sz(rol)/(l^3))*vs) + ((6*E_S(rol)*I_Sz(rol)/(l^2))*ms) + ((12*E_S(rol)*I_Sz(rol)/(l^3))*vs2)+((6*E_S(rol)*I_Sz(rol)/(l^2))*ms2);

G=G+6;i=i+1;
cos(angles(ang));
end
G=G+6; j=j+1;
xxx = linspace(0,length,size(Tz_general,2));

plot(xxx,Tz_general,y{o},'MarkerSize',ys(o))
o=o+1; rol=rol+1;
title('Shear force in Z direction')
xlabel('Length of beam (m)')
ylabel('shear force (N)')
hold on
end
legend('Beam [1]','Beam [2]','Beam [3]','Beam [4]','Beam [5]','Beam [6]','Beam [7]','Beam [8]','Beam [9]','Beam [10]','Beam [11]','Beam [12]')
hold off


%calculate and plot bending moment in y direction
figure;
o=1; G=1; rol=1; ang=1; j=1;
OO_tr=2:Num_msh+1:(((Num_msh+1)*12)-Num_msh+1);
for length=[L L L L L L L L L L L L]
i=1;
l=length/Num_msh;
for ty=1:Num_msh;
% Transformation matrix
l_t=(nodes_coord(OO_tr(j),1))/(L/Num_msh);
m_t=(nodes_coord(OO_tr(j),2))/(L/Num_msh);
n_t=(nodes_coord(OO_tr(j),3))/(L/Num_msh);
D_t=sqrt((l_t^2)+(m_t^2));

Lmbda_3_3=[l_t m_t n_t ; -m_t/D_t l_t/D_t 0 ; (-l_t*n_t)/D_t (-m_t*n_t)/D_t D_t ];
I_t=zeros(3,3);
Tr=[Lmbda_3_3 I_t I_t I_t ; I_t Lmbda_3_3 I_t I_t ; I_t I_t Lmbda_3_3 I_t ; I_t I_t I_t Lmbda_3_3 ] ;

d_dis_shear=[dis_nod(G); dis_nod(G+1) ; dis_nod(G+2);dis_nod(G+3) ;dis_nod(G+4) ;dis_nod(G+5) ...
    ;dis_nod(G+6) ;dis_nod(G+7) ;dis_nod(G+8) ;dis_nod(G+9) ;dis_nod(G+10) ;dis_nod(G+11)] ;   

hof=Tr*d_dis_shear;
vs=hof(3);
vs2=hof(9);
ms=hof(5);
ms2=hof(11);
    
My_general(i)= ((-6*E_S(rol)*I_Sy(rol)/(l^2))*vs) + ((2*E_S(rol)*I_Sy(rol)/(l))*ms) + ((6*E_S(rol)*I_Sy(rol)/(l^2))*vs2)+((4*E_S(rol)*I_Sy(rol)/(l))*ms2);

G=G+6;i=i+1;
cos(angles(ang));
end
G=G+6; j=j+1;
xxx = linspace(0,length,size(My_general,2));

plot(xxx,My_general,y{o},'MarkerSize',ys(o))
o=o+1; rol=rol+1;
title('Bending moment in Y direction n')
xlabel('Length of beam (m)')
ylabel('Bending moment (N.m)')
hold on
end
legend('Beam [1]','Beam [2]','Beam [3]','Beam [4]','Beam [5]','Beam [6]','Beam [7]','Beam [8]','Beam [9]','Beam [10]','Beam [11]','Beam [12]')
hold off

%calculate and plot bending moment in z direction
figure;
o=1; G=1; rol=1; ang=1; j=1;
OO_tr=2:Num_msh+1:(((Num_msh+1)*12)-Num_msh+1);
for length=[L L L L L L L L L L L L]
i=1;
l=length/Num_msh;
for ty=1:Num_msh;
% Transformation matrix
l_t=(nodes_coord(OO_tr(j),1))/(L/Num_msh);
m_t=(nodes_coord(OO_tr(j),2))/(L/Num_msh);
n_t=(nodes_coord(OO_tr(j),3))/(L/Num_msh);
D_t=sqrt((l_t^2)+(m_t^2));

Lmbda_3_3=[l_t m_t n_t ; -m_t/D_t l_t/D_t 0 ; (-l_t*n_t)/D_t (-m_t*n_t)/D_t D_t ];
I_t=zeros(3,3);
Tr=[Lmbda_3_3 I_t I_t I_t ; I_t Lmbda_3_3 I_t I_t ; I_t I_t Lmbda_3_3 I_t ; I_t I_t I_t Lmbda_3_3 ] ;

d_dis_shear=[dis_nod(G); dis_nod(G+1) ; dis_nod(G+2);dis_nod(G+3) ;dis_nod(G+4) ;dis_nod(G+5) ...
    ;dis_nod(G+6) ;dis_nod(G+7) ;dis_nod(G+8) ;dis_nod(G+9) ;dis_nod(G+10) ;dis_nod(G+11)] ;   

hof=Tr*d_dis_shear;
vs=hof(2);
vs2=hof(8);
ms=hof(6);
ms2=hof(12);

Mz_general(i)= ((6*E_S(rol)*I_Sz(rol)/(l^2))*vs) + ((2*E_S(rol)*I_Sz(rol)/(l))*ms) + ((-6*E_S(rol)*I_Sz(rol)/(l^2))*vs2)+((4*E_S(rol)*I_Sz(rol)/(l))*ms2);

G=G+6;i=i+1;
cos(angles(ang));
end
G=G+6;
j=j+1;
xxx = linspace(0,length,size(Mz_general,2));

plot(xxx,Mz_general,y{o},'MarkerSize',ys(o))
o=o+1; rol=rol+1;
title('Bending moment in Z direction n')
xlabel('Length of beam (m)')
ylabel('Bending moment (N.m)')
hold on
end
legend('Beam [1]','Beam [2]','Beam [3]','Beam [4]','Beam [5]','Beam [6]','Beam [7]','Beam [8]','Beam [9]','Beam [10]','Beam [11]','Beam [12]')
hold off


%calculate and plot torsional moment
figure;
o=1; G=1; rol=1; ang=1; j=1;
OO_tr=2:Num_msh+1:(((Num_msh+1)*12)-Num_msh+1);
for length=[L L L L L L L L L L L L]
i=1;
l=length/Num_msh;
for ty=1:Num_msh;
% Transformation matrix
l_t=(nodes_coord(OO_tr(j),1))/(L/Num_msh);
m_t=(nodes_coord(OO_tr(j),2))/(L/Num_msh);
n_t=(nodes_coord(OO_tr(j),3))/(L/Num_msh);
D_t=sqrt((l_t^2)+(m_t^2));

Lmbda_3_3=[l_t m_t n_t ; -m_t/D_t l_t/D_t 0 ; (-l_t*n_t)/D_t (-m_t*n_t)/D_t D_t ];
I_t=zeros(3,3);
Tr=[Lmbda_3_3 I_t I_t I_t ; I_t Lmbda_3_3 I_t I_t ; I_t I_t Lmbda_3_3 I_t ; I_t I_t I_t Lmbda_3_3 ] ;

d_dis_shear=[dis_nod(G); dis_nod(G+1) ; dis_nod(G+2);dis_nod(G+3) ;dis_nod(G+4) ;dis_nod(G+5) ...
    ;dis_nod(G+6) ;dis_nod(G+7) ;dis_nod(G+8) ;dis_nod(G+9) ;dis_nod(G+10) ;dis_nod(G+11)] ;

hof=Tr*d_dis_shear;
vts=hof(4);
vts2=hof(10);
 
tor_general(i)= ( ((-G1_4 * J_b)/l) *vts)  + ( ((G1_4 * J_b)/l) *vts2);

G=G+6;i=i+1;
cos(angles(ang));
end
G=G+6; j=j+1;
xxx = linspace(0,length,size(tor_general,2));

plot(xxx,tor_general,y{o},'MarkerSize',ys(o))
o=o+1; rol=rol+1;
title('Torsion in x direction n')
xlabel('Length of beam (m)')
ylabel('Torsion moment (N.m)')

hold on
end
legend('Beam [1]','Beam [2]','Beam [3]','Beam [4]','Beam [5]','Beam [6]','Beam [7]','Beam [8]','Beam [9]','Beam [10]','Beam [11]','Beam [12]')
hold off
end



if Elasticity_Tensor==true

Volume_uc= 4* sqrt(2)* L^3 ; 
elast=(P_BC_TOTAL_elast' * K_per * P_BC_TOTAL_elast);

C_1111=double((1/(2*Volume_uc))*diff(elast,H_elast(1,1),H_elast(1,1))); 
C_1122=double((1/(2*Volume_uc))*diff(elast,H_elast(1,1),H_elast(2,2))); 
C_1133=double((1/(2*Volume_uc))*diff(elast,H_elast(1,1),H_elast(3,3)));
C_1123=double((1/(2*Volume_uc))*diff(elast,H_elast(1,1),H_elast(2,3)));
C_1113=double((1/(2*Volume_uc))*diff(elast,H_elast(1,1),H_elast(1,3)));
C_1112=double((1/(2*Volume_uc))*diff(elast,H_elast(1,1),H_elast(1,2)));

C_2222=double((1/(2*Volume_uc))*diff(elast,H_elast(2,2),H_elast(2,2)));
C_2233=double((1/(2*Volume_uc))*diff(elast,H_elast(2,2),H_elast(3,3)));
C_2223=double((1/(2*Volume_uc))*diff(elast,H_elast(2,2),H_elast(2,3)));
C_2213=double((1/(2*Volume_uc))*diff(elast,H_elast(2,2),H_elast(1,3)));
C_2212=double((1/(2*Volume_uc))*diff(elast,H_elast(2,2),H_elast(1,2)));

C_3333=double((1/(2*Volume_uc))*diff(elast,H_elast(3,3),H_elast(3,3)));
C_3323=double((1/(2*Volume_uc))*diff(elast,H_elast(3,3),H_elast(2,3)));
C_3313=double((1/(2*Volume_uc))*diff(elast,H_elast(3,3),H_elast(1,3)));
C_3312=double((1/(2*Volume_uc))*diff(elast,H_elast(3,3),H_elast(1,2)));

C_2323=double((1/(2*Volume_uc))*diff(elast,H_elast(2,3),H_elast(2,3)));
C_2313=double((1/(2*Volume_uc))*diff(elast,H_elast(2,3),H_elast(1,3)));
C_2312=double((1/(2*Volume_uc))*diff(elast,H_elast(2,3),H_elast(1,2)));

C_1313=double((1/(2*Volume_uc))*diff(elast,H_elast(1,3),H_elast(1,3)));
C_1312=double((1/(2*Volume_uc))*diff(elast,H_elast(1,3),H_elast(1,2)));

C_1212=double((1/(2*Volume_uc))*diff(elast,H_elast(1,2),H_elast(1,2)));

C_Elasticity_Tensor=[C_1111 C_1122 C_1133 C_1123*sqrt(2) C_1113*sqrt(2) C_1112*sqrt(2) ; ...
     C_1122 C_2222 C_2233 C_2223*sqrt(2) C_2213*sqrt(2) C_2212*sqrt(2); ...
     C_1133 C_2233 C_3333 C_3323*sqrt(2) C_3313*sqrt(2) C_3312*sqrt(2); ...
     C_1123*sqrt(2) C_2223*sqrt(2) C_3323*sqrt(2) 2*C_2323 2*C_2313 2*C_2312; ...
     C_1113*sqrt(2) C_2213*sqrt(2) C_3313*sqrt(2) 2*C_2313 2*C_1313 2*C_1312; ...
     C_1112*sqrt(2) C_2212*sqrt(2) C_3312*sqrt(2) 2*C_2312 2*C_1312 2*C_1212]
S_et=inv(C_Elasticity_Tensor);

for theta_et= 0:0.2:2*pi
for phi_et= 0:0.2:2*pi

n=[ (sin(phi_et)*cos(theta_et))^2  (sin(phi_et)*sin(theta_et))^2  (cos(phi_et))^2   sqrt(2)*sin(phi_et)*sin(theta_et)*cos(phi_et)     sqrt(2)*sin(phi_et)*cos(theta_et)*cos(phi_et)   sqrt(2)*sin(phi_et)*cos(theta_et)*sin(phi_et)*sin(theta_et) ];

E_ym=1/(n*S_et*n');

x_et = E_ym*sin(phi_et)*cos(theta_et);
y_et = E_ym*sin(phi_et)*sin(theta_et);
z_et = E_ym*cos(phi_et);

plot3(x_et,y_et,z_et,'ko','MarkerSize',15,'MarkerFaceColor','r')
title('Youngs modulus (Pa) plotting in spherical coordinates')
hold on
grid on
axis equal
end
end
end



