clear; clc; close all; 

%% Parameters Entry
% Numbering without brackets for Nodes / with brackets [] for Beams

%                          3               2
%                            \           / 
%                             \         /
%                            [3]      [2] 
%                               \     /
%                                \   /            \
%                                 \ /beta)    (gama\
%                4 ------[4]------ 0 ------[1]------ 1
%                                 / \       L
%                                /   \
%                               /     \
%                             [5]     [6]
%                             /         \
%                            /           \
%                          5               6

% The Constructing parameters                                
beta=pi/3;  % Angle (1 0 2) in (rad)                           
gama=pi/3; % Angle (0 1 2) in (rad)                             
coord_SNude=[0,0]; % Coordinate of the starting node    
Num_msh=1; % Number of elements per beam                      
L=1; %the length of the beam in (meter)                                   


% The three imposed tensors 
%         [ eps_11  0 ]         [ 0       0 ]         [ 0   eps_12 ]
%    H_1= |           |  , H_2= |           | ,  H_3= |            |
%         [ 0       0 ]         [ 0   eps_22]         [ eps_12   0 ]
%
%  H_1 -->  Tension or Comprision (+/-) in the horizontal direction
%  H_2 -->  Tension or Comprision (+/-) in the vertical direction
%  H_3 -->  Displacement(Counter clock or clock wise +/-) in the shear direction
H=[1 0 ; 0  0 ];    % Filling as choices from above


% Beam Parameters
%A1=0.785; A2=0.785; A3=0.785; A4=0.785; A5=0.785; A6=0.785; % The cross-sectional area for each beam 
E1_4=2*10^11;   E2_5=2*10^11;   E3_6=2*10^11;   % Modulus of Elasticity for each beam in (Pascale)
% Second moment of area for each beam
b_I=0.2; %Width of the beam in (meter)
h_I=0.2; %Hight of the beam in (meter)

% Multi_Cell plot
Num_hor = 3;
Num_ver = 3;


PostProcessing = true;   % Enable (true) or Disable (false)
Plot_Multi_Cell = true;   % Enable (true) or Disable (false)
Elasticity_Tensor = true;  % Enable (true) or Disable (false)

%% Creating geometry
% Calculation of Geometry Parameters
alpha=pi-(beta+gama);
a=(sin(gama)*L)/(sin(alpha));
b=(sin(beta)*L)/(sin(alpha));

angles=[0 beta pi-gama pi pi+beta (2*pi)-gama ];

j=1;
nodes_coord=coord_SNude; % matrix will contain nudes in the rows and their coordinates in the columns
for length=[L a b L a b] 
    for i=0: Num_msh
      node_coord=[coord_SNude(1)+((length/Num_msh)*i*cos(angles(j))) ,coord_SNude(2)+((length/Num_msh)*i*sin(angles(j)))];  
      nodes_coord=[nodes_coord;node_coord];  
    end 
    j=j+1;
end
nodes_coord(1,:) = []; % Remove the first duplicate row 

%% Imposed Boundary Condition

P_1_dis=[ H(1,:)*[nodes_coord(Num_msh+1,1) nodes_coord(Num_msh+1,2) ]';H(2,:)*[nodes_coord(Num_msh+1,1) nodes_coord(Num_msh+1,2) ]'; 0 ];
P_2_dis=[ H(1,:)*[nodes_coord((Num_msh+1)*2,1) nodes_coord((Num_msh+1)*2,2) ]';H(2,:)*[nodes_coord((Num_msh+1)*2,1) nodes_coord((Num_msh+1)*2,2) ]'; 0 ];                            
P_3_dis=[ H(1,:)*[nodes_coord((Num_msh+1)*3,1) nodes_coord((Num_msh+1)*3,2) ]';H(2,:)*[nodes_coord((Num_msh+1)*3,1) nodes_coord((Num_msh+1)*3,2) ]'; 0 ];                            
P_4_dis=[ H(1,:)*[nodes_coord((Num_msh+1)*4,1) nodes_coord((Num_msh+1)*4,2) ]';H(2,:)*[nodes_coord((Num_msh+1)*4,1) nodes_coord((Num_msh+1)*4,2) ]'; 0 ];                            
P_5_dis=[ H(1,:)*[nodes_coord((Num_msh+1)*5,1) nodes_coord((Num_msh+1)*5,2) ]';H(2,:)*[nodes_coord((Num_msh+1)*5,1) nodes_coord((Num_msh+1)*5,2) ]'; 0 ];                            
P_6_dis=[ H(1,:)*[nodes_coord((Num_msh+1)*6,1) nodes_coord((Num_msh+1)*6,2) ]';H(2,:)*[nodes_coord((Num_msh+1)*6,1) nodes_coord((Num_msh+1)*6,2) ]'; 0 ];                            


%Imposed Boundary Condition (For elasticity tensor calculation)
H_elast = sym('epo',[2 2]);
P_1_dis_elast=[ H_elast(1,:)*[nodes_coord(Num_msh+1,1) nodes_coord(Num_msh+1,2) ]';H_elast(2,:)*[nodes_coord(Num_msh+1,1) nodes_coord(Num_msh+1,2) ]'; 0 ];
P_2_dis_elast=[ H_elast(1,:)*[nodes_coord((Num_msh+1)*2,1) nodes_coord((Num_msh+1)*2,2) ]';H_elast(2,:)*[nodes_coord((Num_msh+1)*2,1) nodes_coord((Num_msh+1)*2,2) ]'; 0 ];                            
P_3_dis_elast=[ H_elast(1,:)*[nodes_coord((Num_msh+1)*3,1) nodes_coord((Num_msh+1)*3,2) ]';H_elast(2,:)*[nodes_coord((Num_msh+1)*3,1) nodes_coord((Num_msh+1)*3,2) ]'; 0 ];                            
P_4_dis_elast=[ H_elast(1,:)*[nodes_coord((Num_msh+1)*4,1) nodes_coord((Num_msh+1)*4,2) ]';H_elast(2,:)*[nodes_coord((Num_msh+1)*4,1) nodes_coord((Num_msh+1)*4,2) ]'; 0 ];                            
P_5_dis_elast=[ H_elast(1,:)*[nodes_coord((Num_msh+1)*5,1) nodes_coord((Num_msh+1)*5,2) ]';H_elast(2,:)*[nodes_coord((Num_msh+1)*5,1) nodes_coord((Num_msh+1)*5,2) ]'; 0 ];                            
P_6_dis_elast=[ H_elast(1,:)*[nodes_coord((Num_msh+1)*6,1) nodes_coord((Num_msh+1)*6,2) ]';H_elast(2,:)*[nodes_coord((Num_msh+1)*6,1) nodes_coord((Num_msh+1)*6,2) ]'; 0 ];
%P_1_dis=[-5;0;0]; 
%%  Creating displacement vector to solve (u) & Entry of Boundary conditios

P_BC_imp=[P_1_dis;P_2_dis;P_3_dis;P_4_dis;P_5_dis;P_6_dis];
P_BC_imp_elast=[P_1_dis_elast;P_2_dis_elast;P_3_dis_elast;P_4_dis_elast;P_5_dis_elast;P_6_dis_elast];

%% plotting the geometry and the mesh
%scatter(nodes_coord(:,1),nodes_coord(:,2),'DisplayName','Red')
scatter(nodes_coord(:,1),nodes_coord(:,2),'bo','LineWidth',2,'MarkerFaceColor','w')
hold on
j=0;
for Num_beams=1:6
    for i=0+j: Num_msh-1+j
    line([nodes_coord((i+1),1),nodes_coord((i+2),1)] , [nodes_coord((i+1),2),nodes_coord((i+2),2)],'LineWidth',2) 
    end 
    j=j+Num_msh+1;
end

% draw the imaginary border of the unit cell
line([L L-(2*b*cos(gama))],[0 (2*b*sin(gama))],'Color','green','LineStyle','--')
line([-L L-(2*b*cos(gama))],[0 (2*b*sin(gama))],'Color','green','LineStyle','--')
line([-L -L+(2*b*cos((2*pi-gama)))],[0 (2*b*sin((2*pi-gama)))],'Color','green','LineStyle','--')
line([L -L+(2*b*cos((2*pi-gama)))],[0 (2*b*sin((2*pi-gama)))],'Color','green','LineStyle','--')
axis equal

%% Calculate the elementary Stiffness Matrix for each beam

Ar=b_I*h_I;
I=(1/12)*b_I*(h_I^3);
Areas=[Ar Ar Ar Ar Ar Ar];
E_S=[E1_4 E2_5 E3_6 E1_4 E2_5 E3_6];
I_S=[I I I I I I];

si=3*(Num_msh+1); % Size of the elementry matrix
d=si-3;  % size of elementry matrix without the common node

j=1;
for length=[L a b L a b]
% Transformation matrix
Tr=sparse([1 1 2 2 3 4 4 5 5 6],[1 2 1 2 3 4 5 4 5 6],[cos(angles(j))  sin(angles(j)) -sin(angles(j))...
 cos(angles(j)) 1  cos(angles(j)) sin(angles(j)) -sin(angles(j)) cos(angles(j)) 1],6,6);

% Local stiffness matrix for  element
Ke_s=sparse([1 1 2 2 2 2 3 3 3 3 4 4 5 5 5 5 6 6 6 6],[1 4 2 3 5 6 2 3 5 6 1 4 2 3 5 6 2 3 5 6],[Areas(j)*E_S(j)/(length/Num_msh) -Areas(j)*E_S(j)/(length/Num_msh) ...
12*E_S(j)*I_S(j)/(length/Num_msh)^3  6*E_S(j)*I_S(j)/(length/Num_msh)^2  -12*E_S(j)*I_S(j)/(length/Num_msh)^3  6*E_S(j)*I_S(j)/(length/Num_msh)^2   6*E_S(j)*I_S(j)/(length/Num_msh)^2  4*E_S(j)*I_S(j)/(length/Num_msh)...
-6*E_S(j)*I_S(j)/(length/Num_msh)^2  2*E_S(j)*I_S(j)/(length/Num_msh)    -Areas(j)*E_S(j)/(length/Num_msh)  Areas(j)*E_S(j)/(length/Num_msh)  -12*E_S(j)*I_S(j)/(length/Num_msh)^3  -6*E_S(j)*I_S(j)/(length/Num_msh)^2   12*E_S(j)*I_S(j)/(length/Num_msh)^3 -6*E_S(j)*I_S(j)/(length/Num_msh)^2 ...
6*E_S(j)*I_S(j)/(length/Num_msh)^2   2*E_S(j)*I_S(j)/(length/Num_msh)   -6*E_S(j)*I_S(j)/(length/Num_msh)^2  4*E_S(j)*I_S(j)/(length/Num_msh)],6,6);

K_e=Tr' * Ke_s * Tr ;     %Local stiffness matrix after Transformation

% Create the stiffness Matrix for one beam
KK=sparse(3*(Num_msh+1),3*(Num_msh+1));
for i=0:6:(6*(Num_msh-1))
KK((i/2)+1:(i/2)+6,(i/2)+1:(i/2)+6)=KK((i/2)+1:(i/2)+6,(i/2)+1:(i/2)+6)+K_e;
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
   else
       K6 = KK;
   end 

j=j+1;
end

%% Assemble of the global stiffness Matrix for the whole unit-cell with descritization nodes

K_glob=sparse(3*(Num_msh+1),3*(Num_msh+1));
K_glob(1:si,1:si)=K1;
K_glob(1:3,1:3)=K1(1:3,1:3)+K2(1:3,1:3)+K3(1:3,1:3)+K4(1:3,1:3)+K5(1:3,1:3)+K6(1:3,1:3);

K_glob(si+1:si+d,si+1:si+d)=K2(4:end,4:end);
K_glob(1:3,si+1:si+d)=K2(1:3,4:end);
K_glob(si+1:si+d,1:3)=K2(4:end,1:3);

K_glob(si+(1*d)+1:si+(2*d),si+(1*d)+1:si+(2*d))=K3(4:end,4:end);
K_glob(1:3,si+(1*d)+1:si+(2*d))=K3(1:3,4:end);
K_glob(si+(1*d)+1:si+(2*d),1:3)=K3(4:end,1:3);


K_glob(si+(2*d)+1:si+(3*d),si+(2*d)+1:si+(3*d))=K4(4:end,4:end);
K_glob(1:3,si+(2*d)+1:si+(3*d))=K4(1:3,4:end);
K_glob(si+(2*d)+1:si+(3*d),1:3)=K4(4:end,1:3);

K_glob(si+(3*d)+1:si+(4*d),si+(3*d)+1:si+(4*d))=K5(4:end,4:end);
K_glob(1:3,si+(3*d)+1:si+(4*d))=K5(1:3,4:end);
K_glob(si+(3*d)+1:si+(4*d),1:3)=K5(4:end,1:3);

K_glob(si+(4*d)+1:si+(5*d),si+(4*d)+1:si+(5*d))=K6(4:end,4:end);
K_glob(1:3,si+(4*d)+1:si+(5*d))=K6(1:3,4:end);
K_glob(si+(4*d)+1:si+(5*d),1:3)=K6(4:end,1:3);

Ke_s=0;
%% periodic

j=1;
for length=[L a b L a b]

% Transformation matrix
Tr=sparse([1 1 2 2 3 4 4 5 5 6],[1 2 1 2 3 4 5 4 5 6],[cos(angles(j))  sin(angles(j)) -sin(angles(j))...
 cos(angles(j)) 1  cos(angles(j)) sin(angles(j)) -sin(angles(j)) cos(angles(j)) 1],6,6);

% Local stiffness matrix for  element
Ke_s=sparse([1 1 2 2 2 2 3 3 3 3 4 4 5 5 5 5 6 6 6 6],[1 4 2 3 5 6 2 3 5 6 1 4 2 3 5 6 2 3 5 6],[Areas(j)*E_S(j)/(length) -Areas(j)*E_S(j)/(length) ...
12*E_S(j)*I_S(j)/(length)^3  6*E_S(j)*I_S(j)/(length)^2  -12*E_S(j)*I_S(j)/(length)^3  6*E_S(j)*I_S(j)/(length)^2   6*E_S(j)*I_S(j)/(length)^2  4*E_S(j)*I_S(j)/(length)...
-6*E_S(j)*I_S(j)/(length)^2  2*E_S(j)*I_S(j)/(length)    -Areas(j)*E_S(j)/(length)  Areas(j)*E_S(j)/(length)  -12*E_S(j)*I_S(j)/(length)^3  -6*E_S(j)*I_S(j)/(length)^2   12*E_S(j)*I_S(j)/(length)^3 -6*E_S(j)*I_S(j)/(length)^2 ...
6*E_S(j)*I_S(j)/(length)^2   2*E_S(j)*I_S(j)/(length)   -6*E_S(j)*I_S(j)/(length)^2  4*E_S(j)*I_S(j)/(length)],6,6);
K_e=Tr' * Ke_s * Tr ;     %Local stiffness matrix after Transformation
% Create the stiffness Matrix for one beam
KK=K_e;
   if j == 1
       K11 = KK;
   elseif j == 2
       K22 = KK;
   elseif j == 3
       K33 = KK;
   elseif j == 4
       K44 = KK;
   elseif j == 5
       K55 = KK;
   else
       K66 = KK;
   end 

j=j+1;
end

% assemble the stiffness matrix but each beam as one element without any discretization point.
K_per=sparse(3*(1+1),3*(1+1));
K_per(1:6,1:6)=K11;
K_per(1:3,1:3)=K11(1:3,1:3)+K22(1:3,1:3)+K33(1:3,1:3)+K44(1:3,1:3)+K55(1:3,1:3)+K66(1:3,1:3);
K_per(6+1:6+3,6+1:6+3)=K22(4:end,4:end);
K_per(1:3,6+1:6+3)=K22(1:3,4:end);
K_per(6+1:6+3,1:3)=K22(4:end,1:3);
K_per(6+(1*3)+1:6+(2*3),6+(1*3)+1:6+(2*3))=K33(4:end,4:end);
K_per(1:3,6+(1*3)+1:6+(2*3))=K33(1:3,4:end);
K_per(6+(1*3)+1:6+(2*3),1:3)=K33(4:end,1:3);
K_per(6+(2*3)+1:6+(3*3),6+(2*3)+1:6+(3*3))=K44(4:end,4:end);
K_per(1:3,6+(2*3)+1:6+(3*3))=K44(1:3,4:end);
K_per(6+(2*3)+1:6+(3*3),1:3)=K44(4:end,1:3);
K_per(6+(3*3)+1:6+(4*3),6+(3*3)+1:6+(4*3))=K55(4:end,4:end);
K_per(1:3,6+(3*3)+1:6+(4*3))=K55(1:3,4:end);
K_per(6+(3*3)+1:6+(4*3),1:3)=K55(4:end,1:3);
K_per(6+(4*3)+1:6+(5*3),6+(4*3)+1:6+(5*3))=K66(4:end,4:end);
K_per(1:3,6+(4*3)+1:6+(5*3))=K66(1:3,4:end);
K_per(6+(4*3)+1:6+(5*3),1:3)=K66(4:end,1:3);

Id=eye(12);
II=eye(9);
Z0=zeros(21,12);
Z0(1:12,1:12)=Id;
Z0(13:end,4:end)=II;

c=[1 0 0 1 0 0 1 0 0 1 0 0 ; 0 1 0 0 1 0 0 1 0 0 1 0];
hh=Z0' * K_per *Z0;
Resttt = rref(hh);
[R,p] = rref(hh);


qq=[0;0;0;P_BC_imp];
qq_elast=[0;0;0;P_BC_imp_elast];

ff= -Z0' * K_per * qq;
ff_elast= -Z0' * K_per * qq_elast;
 
zero=zeros(2,2);
tt=[hh c';c zero];
fff=[ff;0;0];
fff_elast=[ff_elast;0;0];


XX = tt\fff;
XX_elast=tt\fff_elast;

XX=XX(1:12);
XX_elast=XX_elast(1:12);

P_BC_TOTAL=Z0*XX+qq;
P_BC_TOTAL_elast=Z0*XX_elast+qq_elast;

P_BC=P_BC_TOTAL(4:end);

%% Create Reduced system to solve
  k=K_glob;

if Num_msh == 1
   u_red=P_BC_TOTAL;
   f=zeros(21,1);

   for i=1:18

   f=(u_red(end,1)*k(:,end))+f;
   k(end,:)=[];
   k(:,end)=[];
   u_red(end)=[];
   f(end)=[];
   end
   X_red = k\-f;
      X=X_red;
% Define the prescribed boundary nodes 
else 
    hh=1;
    u_red=ones(size(K_glob,2),1)*123.6;
    for i=4+d-3:d:size(K_glob,2)
    u_red(i)=P_BC(hh);
    u_red(i+1)=P_BC(hh+1);
    u_red(i+2)=P_BC(hh+2);
    hh=hh+3;
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
    
   X_red = k\-f;

   X_red(si-2:si)=[];
   X_red(((d+1)+(d-3)):((d+3)+(d-3)))=[];
   X_red(((d+1)+(2*(d-3))):(d+3)+(2*(d-3)))=[];
   X_red(((d+1)+(3*(d-3))):(d+3)+(3*(d-3)))=[];
   X_red(((d+1)+(4*(d-3))):(d+3)+(4*(d-3)))=[];
   X_red(((d+1)+(5*(d-3))):(d+3)+(5*(d-3)))=[];

   X=X_red;
end 

%  Assemble of the total displacement vector
dis_nod=1234.6;
cen_nod=X(1:3);
cv=4;
for j=1:6
    solve=X(cv:cv+(d-4));
    boundary=P_BC(1+((j-1)*3):3+((j-1)*3));
    NEW=[cen_nod ; solve ; boundary];
    dis_nod=[dis_nod ; NEW];
    cv=(cv+(d-4))+1;
end
dis_nod(1)=[];


t=1;
h=4;
dd=2;
for z=1:6
    
for i=2:Num_msh
nodes_coord_after(dd,1)=nodes_coord(dd,1)+X(h);
nodes_coord_after(dd,2)=nodes_coord(dd,2)+X(h+1);
nodes_coord_after(dd,3)=X(h+2);

h=h+3;dd=dd+1;
end
nodes_coord_after(dd,1)=nodes_coord(dd,1)+P_BC(t);
nodes_coord_after(dd,2)=nodes_coord(dd,2)+P_BC(t+1);
nodes_coord_after(dd,3)=P_BC(t+2);
dd=dd+2;t=t+3;
end
hold on
s=1;
for i=1:6   % delete dublicate values
    nodes_coord_after(s,:)=[];
    s=s+Num_msh;
end
nodes_coord_after=[X(1)+XX(1) X(2)+XX(2) X(3);nodes_coord_after];

hold on
sol_plt=scatter(nodes_coord_after(:,1),nodes_coord_after(:,2),'red','LineWidth',2,'MarkerFaceColor','w');

node_after=zeros(size(nodes_coord,1),3);

co=1;
ci=2;
for i=1:6 
    node_after(co,:)= nodes_coord_after(1,:);
    co=co+1;
    for j=1: Num_msh
      node_after(co,:)=nodes_coord_after(ci,:);  
      co=co+1;ci=ci+1;  
    end 
end

%% Plot Interpolation line

countrt=1;
o=1;oo=1; G=1; ang=1; j=1;
for length=[L a b L a b ]
i=1;
l=length/Num_msh;
u=100;x_coords=100;
v=100;y_coords=100;
u_general=0;
v_general=0;
for ty=1:Num_msh;
d_dis_shear=[dis_nod(G); dis_nod(G+1) ; dis_nod(G+2);dis_nod(G+3) ;dis_nod(G+4) ;dis_nod(G+5)  ] ;   

us=d_dis_shear(1);
us2=d_dis_shear(4);
vs=d_dis_shear(2);
vs2=d_dis_shear(5);
ms=d_dis_shear(3);
ms2=d_dis_shear(6);

p_inter=10;
for x=linspace(0,l,p_inter)

mat_pass=[cos(angles(j)) sin(angles(j)) 0 ; -sin(angles(j)) cos(angles(j)) 0 ; 0 0 1];

q1_glo=[us;vs;ms];
q2_glo=[us2;vs2;ms2];
q1_loc=mat_pass*q1_glo;
q2_loc=mat_pass*q2_glo;

u_trac=[q1_loc(1);q2_loc(1)];
vt_flex=[q1_loc(2);q1_loc(3);q2_loc(2);q2_loc(3)];
u_general_local = [1-x/l x/l]*u_trac;
v_general_local = [1-(3*x^2)/l^2+(2*x^3)/l^3 x-(2*x^2)/l+x^3/l^2 (3*x^2)/l^2-(2*x^3)/l^3 -(x^2)/l+(x^3)/l^2]*vt_flex;

new_g=mat_pass'*[u_general_local;v_general_local;0];
u_general=new_g(1);
v_general=new_g(2);

i=i+1;
countrt=countrt+1;
u=[u,u_general];
v=[v,v_general];
end

x_coord=linspace(nodes_coord(o,1),nodes_coord((o+1),1),p_inter);%x_coord(end)=[];
y_coord=linspace(nodes_coord(o,2),nodes_coord((o+1),2),p_inter);%y_coord(end)=[];

x_coords=[x_coords,x_coord];
y_coords=[y_coords,y_coord];

o=o+1;
k_so1=u;

G=G+3;
end
o=o+1; G=G+3; j=j+1;

u(1)=[]; x_coords(1)=[];
v(1)=[]; y_coords(1)=[];

interol_x=u+x_coords;
interol_y=v+y_coords;

intpolation=plot(interol_x,interol_y,'c','LineWidth',3);
hold on
end

title('Unit cell') 
hhhh = zeros(3, 1);
hhhh(1) = plot(NaN,NaN,'or','LineWidth',2,'MarkerFaceColor','w');
hhhh(2) = plot(NaN,NaN,'c','LineWidth',3);
hhhh(3) = plot(NaN,NaN,'ob','LineWidth',2,'MarkerFaceColor','w');
hhhh(4) = plot(NaN,NaN,'--g');

xlabel('X_{Global} ') 
ylabel('Y_{Global} ') 
legend(hhhh, 'Solved Nodes','Interpolation line ','Before solving', 'Cell Border-Before' );
hold off
hold off

x_soso=linspace(nodes_coord(1,1),nodes_coord((2),1),10);

%% Plot Multi-Cells

if Plot_Multi_Cell== true
figure
axis equal
title('Multi-Cell')
sh=0;
sv=0;
for vv=1:Num_ver
for h=0:Num_hor-1
    hold on   
    scatter(nodes_coord_after(:,1)+(((2*L)+(P_BC(1)-P_BC(10)))*h)+sh  ,  nodes_coord_after(:,2)+((P_BC(2)-P_BC(11))*h)+sv,'or','LineWidth',2,'MarkerFaceColor','w');
    G=2;j=1;
for length=[L a b L a b]
u=0;
v=0;
u_general=0;
v_general=0;
l=length/Num_msh;
for ty=2:Num_msh+1;
    i=1;
for  x=linspace(0,l,p_inter)
    
    mat_pass=[cos(angles(j)) sin(angles(j)) 0 ; -sin(angles(j)) cos(angles(j)) 0 ; 0 0 1];

q1_glo=[node_after(G-1,1);node_after(G-1,2);node_after(G-1,3)];
q2_glo=[node_after(G,1);node_after(G,2);node_after(G,3)];
q1_loc=mat_pass*q1_glo;
q2_loc=mat_pass*q2_glo;

u1_m=q1_loc(1);v1_m=q1_loc(2);m1_m=q1_loc(3);
u2_m=q2_loc(1);v2_m=q2_loc(2);m2_m=q2_loc(3);
    ze= x/l;
u_generalll=(1-ze)*u1_m + (ze*u2_m) ;
v_generalll= (1-(3*(ze^2))+(2*(ze^3)))*v1_m +(ze*((1-ze)^2)*l*m1_m) + (((3*(ze^2))-(2*(ze^3)))*v2_m) +((ze^2)*(ze-1)*l*m2_m)   ;

new_q=mat_pass'*[u_generalll;v_generalll;0];

u_general(i)=((new_q(1)+ (((2*L)+(P_BC(1)-P_BC(10))))*h))+sh;
v_general(i)=((new_q(2)+((P_BC(2)-P_BC(11)))*h))+sv;
i=i+1;
end
u=[u,u_general];
v=[v,v_general];
G=G+1;

end
G=G+1;j=j+1;
u(1)=[];
v(1)=[];
intpolation=plot(u,v,'c','LineWidth',3);
hold on
end
      
end
sh=((2*a*cos(beta))+((P_BC(4))-(P_BC(13)))) * vv;
sv=((2*a*sin(beta))+((P_BC(5))-(P_BC(14)))) * vv;
end
xlabel('X_{Global} ') 
ylabel('Y_{Global} ')
axis equal
end

%% Post-Processing
if PostProcessing== true
   
% Calculate shear force
figure;
o=1; y=['+r','*m','oc','b','g','k']; G=1; rol=1; ang=1; j=1;
for length=[L a b L a b ]
i=1;
l=length/Num_msh;
for ty=1:Num_msh;
Tr=[cos(angles(j)) sin(angles(j)) 0 0 0 0 ; -sin(angles(j)) cos(angles(j)) 0 0 0 0 ; 0 0 1 0 0 0 ; ...
    0 0 0 cos(angles(j)) sin(angles(j)) 0; 0 0 0 -sin(angles(j)) cos(angles(j)) 0; 0 0 0 0 0 1 ];    
d_dis_shear=[dis_nod(G); dis_nod(G+1) ; dis_nod(G+2);dis_nod(G+3) ;dis_nod(G+4) ;dis_nod(G+5)  ] ;   

hof=Tr*d_dis_shear;
vs=hof(2);
vs2=hof(5);
ms=hof(3);
ms2=hof(6);
  
T_general(i)= ((-12*E_S(rol)*I_S(rol)/(l^3))*vs) - ((6*E_S(rol)*I_S(rol)/(l^2))*ms) + ((12*E_S(rol)*I_S(rol)/(l^3))*vs2)-((6*E_S(rol)*I_S(rol)/(l^2))*ms2);

G=G+3;i=i+1;
end
G=G+3;
j=j+1;
xxx = linspace(0,length,size(T_general,2));

plot(xxx,T_general,y(o))
o=o+1;
rol=rol+1;
title('Shear force')
xlabel('Length of beam (m)')
ylabel('shear force (N)')
hold on
end
legend('Beam [1]','Beam [2]','Beam [3]','Beam [4]','Beam [5]','Beam [6]')
hold off


% Calculate Normal force
figure;
o=1; y=['+r','*m','oc','b','g','k']; G=1; rol=1; ang=1; j=1;
for length=[L a b L a b]
i=1;
v=0;
l=length/Num_msh;
for ty=2:Num_msh+1;
Tr=[cos(angles(j)) sin(angles(j)) 0 0 0 0 ; -sin(angles(j)) cos(angles(j)) 0 0 0 0 ; 0 0 1 0 0 0 ; ...
    0 0 0 cos(angles(j)) sin(angles(j)) 0; 0 0 0 -sin(angles(j)) cos(angles(j)) 0; 0 0 0 0 0 1 ];    
d_dis=[dis_nod(G); dis_nod(G+1) ; dis_nod(G+2);dis_nod(G+3) ;dis_nod(G+4) ;dis_nod(G+5)  ] ;   

hof=Tr*d_dis;
ss=hof(1);
ss2=hof(4);
    
N_general(i)=((((-E_S(rol)*Areas(rol))/l)*ss) + (((E_S(rol)*Areas(rol))/l)*ss2));

G=G+3;i=i+1;
cos(angles(ang));
end
G=G+3;

xxx = linspace(0,length,size(N_general,2));

plot(xxx,N_general,y(o))
o=o+1;
j=j+1;

rol=rol+1;
title('Normal force')
xlabel('Length of beam (m)')
ylabel('Normal force (N)')
hold on
end
legend('Beam [1]','Beam [2]','Beam [3]','Beam [4]','Beam [5]','Beam [6]')
hold off

% Calculate Bending moment
figure;
o=1; y=['+r','*m','oc','b','g','k']; G=1; rol=1; ang=1; j=1;
for length=[L a b L a b]
i=1; v=0;
l_s=length/Num_msh;
for ty=1:Num_msh;
Tr=[cos(angles(j)) sin(angles(j)) 0 0 0 0 ; -sin(angles(j)) cos(angles(j)) 0 0 0 0 ; 0 0 1 0 0 0 ; ...
    0 0 0 cos(angles(j)) sin(angles(j)) 0; 0 0 0 -sin(angles(j)) cos(angles(j)) 0; 0 0 0 0 0 1 ];    
d_dis=[dis_nod(G); dis_nod(G+1) ; dis_nod(G+2);dis_nod(G+3) ;dis_nod(G+4) ;dis_nod(G+5)  ] ;   

hof=Tr*d_dis;
vs=hof(2);
vs2=hof(5);
ms=hof(3);
ms2=hof(6);

for x_s=0:0.1:l_s;   
M_general(i)= (E_S(rol)*I_S(rol)*((-6*(l_s^-2))+(12*x_s*(l_s^-3)))*vs) + (E_S(rol)*I_S(rol)*((-4*(l_s^-1))+(6*x_s*(l_s^-2)))*ms) + (E_S(rol)*I_S(rol)*((6*(l_s^-2))-(12*x_s*(l_s^-3)))*vs2) + (E_S(rol)*I_S(rol)*((-2*(l_s^-1))+(6*x_s*(l_s^-2)))*ms2);
i=i+1;
end
G=G+3;
cos(angles(ang));
end
G=G+3;

Mx = linspace(1,size(M_general,2),Num_msh);
fg=ceil(Mx);
MM= M_general(1);
for oo=2:Num_msh
MM=[MM M_general(fg(oo))];
end

xxx = linspace(0,length,size(MM,2));

plot(xxx,MM,y(o))
o=o+1;
j=j+1;

rol=rol+1;
title('Bending Moment')
xlabel('Length of beam (m)')
ylabel('Bending Moment (N.m)')
hold on
end
legend('Beam [1]','Beam [2]','Beam [3]','Beam [4]','Beam [5]','Beam [6]')
hold off
end


%% Elasticity Tensor
if Elasticity_Tensor==true

Area_surface= 4*L*a*sin(beta);
elast=(P_BC_TOTAL_elast' * K_per * P_BC_TOTAL_elast);

C_1111=(1/(2*Area_surface))*diff(elast,H_elast(1,1),H_elast(1,1));
C_2222=(1/(2*Area_surface))*diff(elast,H_elast(2,2),H_elast(2,2));
C_1122=(1/(2*Area_surface))*diff(elast,H_elast(1,1),H_elast(2,2));
C_1112=(1/(2*Area_surface))*diff(elast,H_elast(1,1),H_elast(1,2));
C_2212=(1/(2*Area_surface))*diff(elast,H_elast(2,2),H_elast(1,2));
C_1212=(1/(2*Area_surface))*diff(elast,H_elast(1,2),H_elast(1,2));

C_1=double(C_1111);
C_2=double(C_2222);
C_3=double(C_1122);
C_4=double(C_1112);
C_5=double(C_2212);
C_6=double(C_1212);

C_Elasticity_Tensor=[C_1 C_3 C_4*sqrt(2) ; C_3 C_2 C_5*sqrt(2) ;C_4*sqrt(2) C_5*sqrt(2) C_6*2]
S_et=inv(C_Elasticity_Tensor);

i_et=1;
for theta_et= 0:0.01:2*pi
n=[cos(theta_et)^2 sin(theta_et)^2 cos(theta_et)*sin(theta_et)*sqrt(2)];

rr=n*S_et*n';
E_ym(i_et)=1/rr;

i_et=i_et+1;
end

figure;
theta_plt= 0:0.01:2*pi;
polarplot(theta_plt,E_ym,'linewidth',2)
title('Young modulus (Pa) ')
% xlabel('Angle')
% ylabel('Young modulus (pa)')
pax = gca;
pax.ThetaAxisUnits = 'radians';

end



