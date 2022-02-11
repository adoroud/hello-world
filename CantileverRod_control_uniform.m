%% Inverse Dynamic Control of Cantilever Rod using Cosserat-Rod model (Actuation Forces + Gravity Force + External Forces)
%Azadeh version created 2/9/22
%Bending for grasping-Silicone

%visialization: finally a uniform and tapering 3D arm with 3 physical
%segments

%original Cosserat Code from Paper John Till, Caleb Rucker 2018
% Real-Time Dynamics of Soft and Continuum Robots based on Cosserat-Rod Models

function output= CantileverRod_control_uniform
%% Parameters - Section 2.1
global fseg lseg rhoAg C q t result ud vd state_out segment N kapa_out lseg
% N L STEPS 
%% Video recording ()
% myVideo = VideoWriter('grasp'); %open video file
% myVideo.FrameRate = 1;  %can adjust this, 5 - 10 works well for me
% open(myVideo)
%% Time descretization
del_t = 1; %0.025; % Time step= 0.01 is good ds=5ms
% del_t = 0.0336; %0.0250
STEPS = 100;
t= zeros(STEPS,1);
%% Space discretezition
scale = 1;
physical_segment = 1;
decretized_sections = 10;
N = physical_segment*decretized_sections;      %40 ; Number of discretization for rod = 4 * number of nodes per segment
L = physical_segment*0.185*scale; % Length of rod %L>>r cross-section diameter or sqrt(A)
ds = L/N;    % Discretization of curve parameter and step size in RK4
% ds = del_t;
%% Main input and outputs
result = zeros(STEPS,19,N+1);
state_out = zeros(STEPS,6,N+1);
% seg_num = zeros(physical_segment,1);
circle = zeros(301,3);  %visualization cylindrical view
p_c = zeros(301,3,N+1);  %surface of the cylinder
% pressure = 0.01*20*6894.76;  %psi to N/m2
% pressure_2 = 0.04*20*6894.76;  %psi to N/m2
%% vectors directions in 3D
e1 = [ 1 ; 0 ; 0 ];
e2 = [ 0 ; 1 ; 0 ];
e3 = [ 0 ; 0 ; 1 ];
%% desired motions-initial
ud = zeros(STEPS,3,N+1);
vd = zeros(STEPS,3,N+1);

kapa = zeros(1,N+1);
% u_xd = zeros(1,N+1);  %bending element of strain
% ut_xd = zeros(1,N+1);
% utt_xd = zeros(1,N+1);
% 
% u_yd = zeros(1,N+1);  %bending element of strain
% ut_yd = zeros(1,N+1);
% utt_yd = zeros(1,N+1);
% 
% u_zd = zeros(1,N+1);  %torsion element of strain
% ut_zd = zeros(1,N+1);
% utt_zd = zeros(1,N+1);

u_d = zeros(3,N+1);
ut_d = zeros(3,N+1);
utt_d = zeros(3,N+1); 

v_d = zeros(3,N+1); 
vt_d = zeros(3,N+1); 
vtt_d = zeros(3,N+1); 
%% External forces and torques to tip (free-end)- global frame
F_tip = zeros(STEPS,3); %0.001*9.81*[0 ; 0 ; 0 ]; % Force on tip based on gram
M_tip = zeros(STEPS,3); %[ 0 ; 0 ; 0 ]; % Momentum on tip
g = 1*[0 ; 0 ; -9.81] ; % Gravitational force
%% Boundary & Initial Conditions - Section 2.4   global frame
p0 = [ 0 ; 0 ; 0 ]; % Global position in cartesian coordinates
h0 = [ 1 ; 0 ; 0 ; 0 ]; % Quaternions for material orientation
q0 = [ 0 ; 0 ; 0 ]; % Velocity (linear) in the local frame (w.r.t time)
w0 = [ 0 ; 0 ; 0 ]; % Angular velocity in the local frame (w.r.t time)
n0 = [ 0 ; 0 ; 0 ]; % initial global internal force
m0 = [ 0 ; 0 ; 0 ]; % initial global internal moment
vstar = [ 0 ; 0 ; 1 ];    % initial position velocity w.r.t arc length (local frame)
ustar = [ 0 ; 0 ; 0 ];   % initial curvature for a straight rod (rotation derivative w.r.t arc length) (local frame)
%% Dependent Parameter Calculations
rho = 800;   %507kg/m3 % 2330; % 1.3468e3;    % Hydrogel density=1.1 - Beam density=M(3e-3)/V(2.2275e-6)=1.3468e3
% rhow = 1000;

% C = 0.0041*eye(3); % Square law drag air resistance- in our case damping coefficient in water
% C=0.0262*eye(3);  %old value
C=zeros(3,3);

% r0 = 0.003*scale;         % Hydrogel voxel center of mass distance to the backbone
% r = 0.003*e2;    % Hydrogel voxel center of mass distance to the backbone
% r = 0.003*e3;    % Hydrogel voxel center of mass distance to the backbone
r1 = 0.034*[cos(-pi/4);sin(-pi/4);0];   %pressure tube number 1 distance to the backbone
r2 = 0.034*[cos(pi/4);sin(pi/4);0];   %pressure tube number 2 distance to the backbone
%% uniform arm
r0 = 0.075*sqrt(2)/2;
% beam_height = 1e-3*scale; %backbone or beam 11mm total %x direction
% beam_width = 4.5e-3*scale; %y direction
% beam_length = L;  %z direction
%% Tapering arm
% rp=r0;  %proximal segment radius
% rd=0.01;  %distal segment radius
%% Cross-section area
A0 = pi*r0^2;  %circular  %rectangular: beam_height*beam_width; % Cross-sectional area
% radius = sqrt(A); %L>>r cross-section diameter approximation or sqrt(A)
% delta_L = L/radius;  %beam_length/beam_width
A =A0; %in the silicone arm this won't change significantly and it is passive
%% Young modulus of Hydrogel beam, was E=(L=0.045/A=11*4.5e-6)*(2.2061)=2.005e3
% Hydrogel from Ximin data 5kPa at rest, 30kPa at actuation
% New data from Roozbeh, diff hydrogels from 20kPa to 80kPa at rest
% E for Rubber=10-100 MPa, Nylon=2-4 GPa, 
%E = (L/A)*(0.0542)*1;   %E=271 when L is considered for stiffness (0.0542 in mm scale)
Em = 0.28*10^6; %188*10^9;   %0.3*1000*scale^2*(0.0542*scale)/A0;  %E=6000 when L is not considered for rigidity
Gm = Em/(2*(1 + 0.4));   %Shear modulus

Ixx = (pi/4)*r0^4;  %circular   %rectangular: (beam_height*beam_width^3)/12;
Iyy = (pi/4)*r0^4;  %circular   %rectangular: (beam_width*beam_height^3)/12;
Izz=Ixx+Iyy;
J = diag([Ixx  Iyy  Izz]); % Second mass moment of inertia tensor

Kse = diag([Gm*A0, Gm*A0, 10*Em*A0]); % Stiffness matrix for shear(d1,d2) rigidity and extension(d3) or axial rigidity
Kbt = diag([Em*J(1,1), Em*J(2,2), Gm*J(3,3)]); % Stiffness matrix for bending(d1,d2)rigidity and twisting(d3) rigidity

% Bse = zeros(3); % viscous damping matrix for shear & extension
% Bbt = zeros(3); % viscous damping matrix for bending & torsion
tau=0; %2*1.5; %
Bse = Kse* diag([tau, tau, tau]);
Bbt = Kbt* diag([tau, tau, tau]);
%% BDF (backward diffrentiation formula) Coefficients 
%%(an implicit method of time derivative approximation)
%%Discretization in time for PDE (Section 2.2 Eq.(5)-b)
alpha=0;    %BDF2
% alpha=-0.5; %Trapezoidal
% alpha=-0.2; %found to be accurate and stable in Till 2019 paper
%d1 = alpha/(1+alpha); %alpha >> d1
c0 = (1.5+alpha)/(del_t*(1+alpha)); %c0 = 1.5/del_t; 
c1 = -2/del_t;
c2 = (0.5+alpha)/(del_t*(1+alpha)); %c2 = 0.5/del_t; 
%% Expressions extracted from simulation loop
Kse_plus_c0_Bse_inv = (Kse+c0*Bse)^-1; % Eq.(7)-b
Kbt_plus_c0_Bbt_inv = (Kbt+c0*Bbt)^-1; % Eq.(7)-b
Kse_vstar = Kse*vstar; % Eq.(7)-b
Kbt_ustar = Kbt*ustar; % Eq.(7)-b
rhoA = rho*A0; % Eq.(7)-a
rhoJ = rho*J; % % Eq.(1)
rhoAg = rho*A0*g; % Eq.(7)-a
% rhowAg = (rho - rhow)*A0*g;  %with water drag effect
%% Initialize to straight configuration
% y and z are general variables as in Eq.(12) >> Generalized PDE
% y : = [ p ; h ; n ; m ; q ; w ] and z : = [ v ; u ]
y = [zeros(2,N+1); linspace(0,L,N+1); zeros(16,N+1)] ;
z = [zeros(2,N+1); ones(1,N+1); zeros(3,N+1)];
y_prev = y;
z_prev = z;
%% Main Simulation Loop - Section 2.4
% visualize(); % Plot initial rod
% plotRobot();
% Force_sensor = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
% Force_sensor = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
G = zeros(6,1) ; % Shooting method initial guess
%%
    for i = 1:STEPS
    
    t(i)= (i-1) * del_t;
    %% 4 motions
    %     [u_d,ut_d,utt_d,v_d,vt_d,vtt_d]= grasp(t(i));
    %     [u_d,ut_d,utt_d,v_d,vt_d,vtt_d]= twist(t(i));        
    %     [v_d,vt_d,vtt_d]= shear(t(i));
    %     [u_d,ut_d,utt_d,v_d,vt_d,vtt_d]= extension(t(i));       
    %% main original code solving the coesserat rod forward dynamics
        % Set history terms - Eq.(5)-b
        % y_prev ~ _()^(i-2)y(t)
        yh = c1*y + c2*y_prev;
        zh = c1*z + c2*z_prev;
        
        %store from previous step for the next step
        y_prev = y;
        z_prev = z;
        
        % Midpoints are linearly interpolated for RK4
        yh_int = 0.5*(yh(:, 1:end-1) + yh(:, 2:end));
        zh_int = 0.5*(zh(:, 1:end-1) + zh(:, 2:end));
        
        %Standard Simple Shooting Method (SSM) using fsolve
          
        G = fsolve(@getResidual,G); 
        
        result(i,:,:) = y;    %saving results
        state_out(i,:,:) = z;  %saving u,v output during time and space
%         kapa_out(i,:) = kapa;
        plotRobot();
    end
%% getResidual Function
    function E = getResidual(G)   %takes initial guess and calculate E, give it to next step as a new guess
        %Reaction force and moment are guessed values
        n0 = G(1:3); m0 = G(4:6);
        initial_guess(i,:) = G;
        %Cantilever boundary condition-initial point
        y(:,1) = [p0; h0; n0; m0; q0; w0] ; 
        %Fourth-Order Runge-Kutta Integration 
        %(implicit integrating scheme) semi-discretized model: space not time
        for j = 1:N
            yj = y(:,j) ; %yhj_int = yh_int(:,j); %ythj_int = yth_int(:,j);
            [k1, z(:,j+1), ~] = ODE(yj, yh(:,j), zh(:,j),j);
            [k2, ~, ~] = ODE(yj + k1*ds/2, yh_int(:,j), zh_int(:, j), j);
            [k3, ~, ~] = ODE(yj + k2*ds/2, yh_int(:,j), zh_int(:, j), j);
            [k4, ~, R] = ODE(yj + k3*ds , yh(:,j+1), zh(:,j+1), j);
            %where y is calculated using RK4 spatial integration
            y(:, j+1) = yj + ds*(k1 + 2*(k2+k3) + k4)/6 ;
            %y(:, j+1) = yj + ds*k1; % Euler's Method (explicit time integrating method)
            
            t_c = 0:pi/150:2*pi;

        for counter = 1:length(t_c)
            x_c = r0*cos(t_c(counter));
            y_c = r0*sin(t_c(counter));
            z_c = 0;
            circle(counter,:) = [x_c;y_c;z_c]; 
            p_c(counter,:,j+1) = (y(1:3,j+1)+R*circle(counter,:)');
        end
%inverse kinematic funtion being called here to give theta value of segment number j    
%             if mod(j,N/physical_segment)==0
%                segment = j/(N/physical_segment);
%                 inv_kin(y(1:3,j+1),R,segment);
%             end
         end
       %disp('iteration') 
       %At each time and space step, inside the PDEs/ODEs more solve inverse dynamics
       %Here by solving the RK4 for PDes/ODEs at each time step solve forward dynamics 
       %from initial guess of initial space point to the tip and by constrainting the 
       %free-end boundary condition the output provides data for fsolve (BVP) to find 
       %the new guess for the next timing step
       %Cantilever boundary conditions @ end point -tip
       %residual error of the tip for states n and m as force and moment at 
       %free-end that should be kept at zero except any external values that may apply. 
       %The fsolve will be solved for E=0, then both BVP will be
       %satisfied.(TPBVP:two point boundary value problem)
        nL = y(8:10, N+1); mL = y(11:13 ,N+1);
        F_tip(i,:) = [0;0;0]'+ rhoAg'- (R*C*q.*abs(q))'+ fseg(i,:,end);  %0.001*9.81*[0;0;0];  %rhoAg' +
        M_tip(i,:) = [0;0;0]'+ lseg(i,:,end);  %0.001*[0.01;0;0];
        E = [F_tip(i,:)'-nL; M_tip(i,:)'-mL]; 
        free_end(i,:) = E;
    end

%% ODE Function Eq.(12)-b - solving PDE/ODE 
    function [ys, z, R] = ODE(y, yh, zh, node)
        %% discritization and segment relationship based on 4 real segments
        if mod(node,N/physical_segment)==0
            segment = node/(N/physical_segment);
            seg_num(segment,1) = node;
        else
            segment = 1+(node-mod(node,N/physical_segment))/(N/physical_segment);
        end
        %% defining initial input data to the segment for inverse loop coming from previous forward loop
        %at time step i for seg j
        p = y(1:3);
        h = y(4:7); 
        n = y(8:10); 
        m = y(11:13);
        q = y(14:16); 
        w = y(17:19);
        vh = zh(1:3); 
        uh = zh(4:6);
        %Quaternion to Rotation - Eq.(10)- Kinematics - Task space variables
        h1 = h(1); h2 = h(2); h3 = h(3); h4 = h(4);
        R = eye(3)+ 2/( h'*h )* ...
            [-h3^2- h4^2 , h2*h3 - h4*h1 , h2*h4+h3*h1;
            h2*h3+h4*h1 , -h2^2 - h4^2 , h3*h4 - h2*h1;
            h2*h4 - h3*h1 , h3*h4 + h2*h1 , -h2^2 - h3^2];

        %% this is actually the main forward dynamic calculation (force/moment to motion)
        %Solved - Eq.(6)- Constitutive Law - Configuration space variables-local frame
        v = Kse_plus_c0_Bse_inv*(R'*n + Kse_vstar - Bse*vh);  
        u = Kbt_plus_c0_Bbt_inv*(R'*m + Kbt_ustar - Bbt*uh);  
%         if y(1)==0
%             kapa=0;
%         else
            p_local = R'*p;
            kapa(node+1) = -2*p_local(1)/(p_local(1)^2+p_local(3)^2);   %curvature_about_x_axis
%             % forward solution by the robot
%             u = [0; -kapa(node+1); 0];
%             v = [0; 0; (dL)/ds];
%         end
%         error_of_curvature_about_x_axis = u(2)-4*kapa;
        %% the rest...
        z = [v; u];
        %Time Derivatives - Eq.(5)-- BDF-alpha- an implicit method
        %using backward diffrentiation method 
        yt = c0*y + yh; % derivative wrt time
        zt = c0*z + zh; % derivative wrt time
        vt = zt(1:3); % derivative wrt time
        ut = zt(4:6); % derivative wrt time
        qt = yt(14:16); % derivative wrt time
        wt = yt(17:19); % derivative wrt time
        %% The final input to the robotic arm 
        %%these data are the arc length spatial derivatives in order to be used in RK4 integrator
        %Rod State Derivatives - Eq.(7) - Dynamics
        %diff. kinematics calculations for final Cosserat PDE equations
        ps = R*v;
        Rs = R*[0, -u(3), u(2);
                u(3), 0, -u(1);
                -u(2), u(1), 0];  %cross(R,u) or R*uhat;
            
        %% Weight and Square-Law-Drag - Eq.(3)
%          f = rhoAg; 
         fenv = rhoAg - R*C*q.*abs(q); %without water drag effect  
%          fenv = rhowAg - R*C*q.*abs(q);   %with water drag effect
         lenv = cross([0;0;0],fenv); %zero distance as long as it is applied to the backbone directly
        %% when external loads being applied (like from the object being detected by FSR)
        fext = [0;0;0];   %external force to the segment (constant and local)
        lext = cross(r0*[1;0;0],R*fext);   %external moment to the segment (constant and global) in desired direction
       %% open loop control 20psi step input
%         fseg(i,:,segment) = [0;0;0]; %-(pressure_1*A*Rs*e3 + pressure_2*A*Rs*e3);
%         lseg(i,:,segment) =  R*(cross(r1+r2,pressure*A*e3)); %this is not per unit length

%         fseg(i,:,segment) = -(pressure_1*A*Rs*e3 + pressure_2*A*Rs*e3);   
%         lseg(i,:,segment) =  -(pressure_1*A*R*(cross((v+cross(u,r1)),e3)+cross(r1,cross(u,e3))) + pressure_2*A*R*(cross((v+cross(u,r2)),e3)+cross(r2,cross(u,e3))));
       %% PD control-shear -part of solving inverse dynamics- global frame- where we close the loop
%         [u_d,ut_d,utt_d,v_d,vt_d,vtt_d]=grasp(segment);
        grasp();
        ud(i,:,node+1) = u_d;
        vd(i,:,node+1) = v_d;
        fseg(i,:,segment) = -R*( Kse*(v_d-v) + Bse*(vt_d-vt) + vtt_d) - fenv -fext;  
        lseg(i,:,segment) = -30*R*( Kbt*(u_d-u) + Bbt*(ut_d-ut) + utt_d) -lenv -lext; %
%         fseg(i,:,segment) = R*[(4/3)*-70 0 0;0 (4/3)*-5*0.001 0;0 0 1]*( Kse*(v_d(:,segment)-v) + Bse*(vt_d(:,segment)-vt) + vtt_d(:,segment)) - fenv -fext;  
%         lseg(i,:,segment) = R*[1 0 0;0 1 0;0 0 1]*( Kbt*(u_d(:,segment)-u) + Bbt*(ut_d(:,segment)-ut) - utt_d(:,segment)) -lext;
        %         [0.001/2 0 0;0 0.01 0;0 0 1]  %for sin case
%         [120 0 0;0 50 0;0 0 1]  %for const case
% [(4/3)*-70 0 0;0 (4/3)*-50 0;0 0 12] %for extension
%I have to add a function here that translate the controller input to the low-level controller to control actuators
% the input is f & l, the output should be sent to the robot
        %% this is actually the main inverse dynamic calculation (motion to force/moment)
        % n=R*(Kse*(v-vstar)+Bse*vt) %global frame
        % m=R*(Kbt*(u-ustar)+Bbt*ut) %global frame
%% calculating the input to the segment for the next loop considering all dynamics
       
        ns = R*rhoA*(cross(w,q)+qt) - fseg(i,:,segment)' -fenv - fext;
        ms = R*( cross(w, rhoJ*w) + rhoJ*wt) - cross(ps, n) - lseg(i,:,segment)' - lext -lenv;
   
        %compatibility
        qs = vt - cross(u , q) + cross(w, v) ;
        ws = ut - cross(u ,w) ;
        %Quaternion Derivative - Eq.(9)
        hs = [0, -u(1), -u(2), -u(3);
              u(1), 0, u(3), -u(2);
              u(2), -u(3), 0 , u(1);
              u(3), u(2), -u(1), 0]*h/2;
        ys = [ps; hs; ns; ms; qs; ws];   %the segments' tip data at this time step
    end
%% Grasp related input definition, switching conditions based on sensor feedback 
    function grasp()
      
    %             lambda = L;%/2;0.020;%(N/2)*ds;          
    %             k_l=2*pi/lambda;
%         alpha= 8/physical_segment;
        T_period = 100;
        frequency = 1/T_period; 

%         if Force_sensor(i)==0
%            contact= 0;
%         elseif Force_sensor(i)==1
%            contact= 1;
%         end

%     for k=1:N+1  %useful for assigne different amounts of input to each section or segment

             if segment==1

                 alpha= 1; %6/physical_segment;

                 u_xd =  0; %alpha*(sin(2*pi*frequency*t)); %+pi/2 %angle of changes
                 ut_xd =  0; %alpha*(2*pi*frequency)*cos(2*pi*frequency*t);
                 utt_xd = 0; %-alpha*(2*pi*frequency)^2*sin(2*pi*frequency*t);
        %          
                 u_yd =  alpha*(sin(2*pi*frequency*t(i))); %exp(alpha*t(i)); %
                 ut_yd = alpha*(2*pi*frequency)*cos(2*pi*frequency*t(i)); %alpha*exp(alpha*t(i));
                 utt_yd = -alpha*(2*pi*frequency)^2*sin(2*pi*frequency*t(i));%alpha*alpha*exp(alpha*t(i));

                 u_zd =  0; %alpha*sin((2*pi)*(frequency*t(i)))*cos((2*pi)*(frequency*t(i) - k*ds/(lambda)));
                 ut_zd =  0; %alpha*(2*pi*frequency)*cos((2*pi)*(frequency*t(i)))*cos((2*pi)*(frequency*t(i) - k*ds/(lambda)))+sin((2*pi)*(frequency*t(i)))*(-2*pi*frequency)*sin((2*pi)*(frequency*t(i)-k*ds/(lambda)));
                 utt_zd = 0; %-alpha*(2*pi*frequency)^2*sin((2*pi)*(frequency*t(i)));
                 
                 v_zd = 1; %+0.2*t(i)/STEPS; 
                 vt_zd = 0; %.2/STEPS;
                 vtt_zd = 0; 

             elseif segment==2

                 alpha= 5; %6/physical_segment;

                 u_xd =  0; %alpha*(sin(2*pi*frequency*t)); %+pi/2 %angle of changes
                 ut_xd =  0; %alpha*(2*pi*frequency)*cos(2*pi*frequency*t);
                 utt_xd = 0; %-alpha*(2*pi*frequency)^2*sin(2*pi*frequency*t);
        %          
                 u_yd =  alpha*(sin(2*pi*frequency*t(i))); %exp(alpha*t(i)); %
                 ut_yd = alpha*(2*pi*frequency)*cos(2*pi*frequency*t(i)); %alpha*exp(alpha*t(i));
                 utt_yd = -alpha*(2*pi*frequency)^2*sin(2*pi*frequency*t(i));%alpha*alpha*exp(alpha*t(i));

                 u_zd =  0; %alpha*sin((2*pi)*(frequency*t(i)))*cos((2*pi)*(frequency*t(i) - k*ds/(lambda)));
                 ut_zd =  0; %alpha*(2*pi*frequency)*cos((2*pi)*(frequency*t(i)))*cos((2*pi)*(frequency*t(i) - k*ds/(lambda)))+sin((2*pi)*(frequency*t(i)))*(-2*pi*frequency)*sin((2*pi)*(frequency*t(i)-k*ds/(lambda)));
                 utt_zd = 0; %-alpha*(2*pi*frequency)^2*sin((2*pi)*(frequency*t(i)));

                 v_zd = 1; %1+0.2*t(i)/STEPS;
                 vt_zd = 0; %0.2/STEPS;
                 vtt_zd = 0; %0;

             elseif segment==3

                 alpha= 0.6;%6/physical_segment;

                 u_xd=  0; %alpha*(sin(2*pi*frequency*t)); %+pi/2 %angle of changes
                 ut_xd =  0; %alpha*(2*pi*frequency)*cos(2*pi*frequency*t);
                 utt_xd = 0; %-alpha*(2*pi*frequency)^2*sin(2*pi*frequency*t);
        %          
                 u_yd =  alpha*(sin(2*pi*frequency*t(i))); %exp(alpha*t(i)); %
                 ut_yd = alpha*(2*pi*frequency)*cos(2*pi*frequency*t(i)); %alpha*exp(alpha*t(i));
                 utt_yd = -alpha*(2*pi*frequency)^2*sin(2*pi*frequency*t(i));%alpha*alpha*exp(alpha*t(i));

                 u_zd =  0; %alpha*sin((2*pi)*(frequency*t(i)))*cos((2*pi)*(frequency*t(i) - k*ds/(lambda)));
                 ut_zd=  0; %alpha*(2*pi*frequency)*cos((2*pi)*(frequency*t(i)))*cos((2*pi)*(frequency*t(i) - k*ds/(lambda)))+sin((2*pi)*(frequency*t(i)))*(-2*pi*frequency)*sin((2*pi)*(frequency*t(i)-k*ds/(lambda)));
                 utt_zd = 0; %-alpha*(2*pi*frequency)^2*sin((2*pi)*(frequency*t(i)));

                 v_zd = 1; %1+0.2*t(i)/STEPS;
                 vt_zd = 0; %0.2/STEPS;
                 vtt_zd =0; %0;
             end
%     end

    %       for k=1:N  
    %          %%exactly the same as Roozbeh's travel wave experiment (continuous signal)
    %          %%the change is from 1 to 0.3 because it is dz not voltage and we
    %          %%have shrinking limit of one-third of the max
% %              v_zd = 1+0.2*t(i)/STEPS; %0.65 + 0.35*cos((2*pi) * ( frequency*t(i) - k*ds/(lambda)));
% %              vt_zd =  0.2/STEPS; %-0.35*(2*pi*frequency)*sin((2*pi) * ( frequency*t(i) - k*ds/(lambda)));
% %              vtt_zd = 0; %-0.35*(2*pi*frequency)^2*cos((2*pi) * ( frequency*t(i) - k*ds/(lambda)));
    %       end
    %% input command for bending (u_x,u_y) and torsion (u_z)-one at an experiment
        u_d = [u_xd;u_yd;u_zd];
        ut_d = [ut_xd;ut_yd;ut_zd];
        utt_d = [utt_xd;utt_yd;utt_zd];
        
        v_d = [0;0;v_zd];  %1.1*t(i)/STEPS
        vt_d = [0;0;vt_zd];
        vtt_d = [0;0;vtt_zd];
    end
%% Visualization of the robotic arm been drawn in real-time 
    function plotRobot()
        p_c(:,:,1) = circle(:,:);   %first segment's first point
    %     p_c(:,:,N+1) = p_c(:,:,N);    %last segment's last point
        for j=1:N+1
            p_c1(j,:) = p_c(:,1,j)';
            p_c2(j,:) = p_c(:,2,j)';
            p_c3(j,:) = p_c(:,3,j)';
%         end
%         for j=1:N+1
    
            if j<= N/physical_segment
                
                 plot3(100*p_c1(j,:), 100*p_c2(j,:), 100*p_c3(j,:),'Color',[1 0 0],'LineWidth',3);
    
            elseif j<= 2*N/physical_segment
                
                 plot3(100*p_c1(j,:), 100*p_c2(j,:), 100*p_c3(j,:),'Color',[0 0 1],'LineWidth',3);
                 
            elseif j<= 3*N/physical_segment
    
                 plot3(100*p_c1(j,:), 100*p_c2(j,:), 100*p_c3(j,:),'Color',[1 0 0],'LineWidth',3);
                 
            else
    
                 plot3(100*p_c1(j,:), 100*p_c2(j,:), 100*p_c3(j,:),'Color',[0 0 1],'LineWidth',3);
                 
            end          
            hold on;
            if mod(j-1,N/physical_segment)==0
                    plot3(100*p_c1(j,:), 100*p_c2(j,:), 100*p_c3(j,:),'Marker','o','Markersize',4,'MarkerEdgeColor','k'); %'Color',[1 0 0],'LineWidth',1,
            end
        end    
    %             hold on;
    %             hold off;
                grid on;          
                ax = gca;
                ax.XDir = 'reverse';
                ax.YDir = 'reverse';
                ax.ZDir = 'normal';
                axis([-50 50 -50 50 0  100]) ;
                xlabel('X(cm)');ylabel('Y(cm)');zlabel('Z(cm)');
    %           daspect ([ 2 1 1 ]); 
                view(0,0);
    %             view(15,10);
    %             view(-150,20);
                drawnow;
    %             if i==1
    %                 savefig('shearg1.fig')
    %             elseif i==3
    %                 savefig('shearg2.fig')
    %             elseif i==8
    %                 savefig('shearg3.fig')
    %             elseif i==10
    %                 savefig('shearg4.fig')
    %             end
                hold off;
        %get frame
%         frame = getframe(gcf);
%         writeVideo(myVideo, frame);
    end
%   close(myVideo)
  output= result;
end