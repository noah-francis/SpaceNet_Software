%%
clc;clear;close all;
%% MATLAB example PF from website
%11/8/2020
%Keith Poletti implementing
%% Initialize elements

load CSIM_Position&Veloctiy_ECIJ2000.mat
load ECI_True.mat
load TDoA_ECI_out.mat
load TDoA_DOP.mat

%Load in the true values for the entire orbit
State_Orbit = cell2mat(finalCell);
Dates = State_Orbit(:,1:6);
%Creates the state Vector for the entire orbit
xTrueOrbit = [State_Orbit(:,7:9),State_Orbit(:,11:13)];


% Load in the True values at each measurement
TofMeas = ECI_True(:,1:6); %Time of Measurement
R_TrueMeas = ECI_True(:,7:9)*10^(-3);%True Position when measurement happens
vTrueMeas = ECI_True(:,11:13)*10^(-3);%True Velocity when measurement happens
xTrueMeas = [R_TrueMeas, vTrueMeas];
indMeas = ECI_True(:,end); % Indices in the entire orbit state where a measurement happens


%Load in the measurements
yMeas=ECI_Tdoa(:,7:9)*10^(-3);
sig_XYZ=(yMeas-xTrueMeas(:,1:3));
k=find(yMeas(:,1)>10^4);
% yMeas(k,:)=(yMeas(k-1,:)+yMeas(k+1,:))/2;
sig_r=ECI_Tdoa(:,14);%*10^3;


% For debugging
% yMeas=R_TrueMeas;


rng(1); % Fix the random number generator for reproducible results
%Add process Noise to True state?
% xTrue = xTrue; %+ 0.1^2 * randn(length(timeVector),2);

%Add Measurement Noise
% yTrue = xTrue(:,1);
% yMeas = yTrue + sqrt(R) * randn(size(yTrue)); % sqrt(R): Standard deviation of noise


%Vector to initialize the particle Filter with in 3 sigma bounds

R_Init= [(yMeas(1,:) - 3*sig_r(1))', ...
         (yMeas(1,:) + 3*sig_r(1))'];
     
     
 

%CHEATERS way of admissble region DO NOT USE IN FINAL VERSION

%%
[R_Gibbs ,V_ECI] = ODGibbsMethod([yMeas(7,:)',yMeas(17,:)',yMeas(32,:)']);
 xTrueMeas(16,4:6)
 Vstd=std(diff(yMeas(1:32,:))/10)';
[~ ,V_ECIerr] = ODGibbsMethod([yMeas(1,:)',yMeas(2,:)',yMeas(3,:)']);
%%
V_CHEAT=[V_ECI-3*Vstd,V_ECI+3*Vstd];
% [val,ind]=min(abs(vTrueMeas),[],1);
% V_CHEAT=[sign(vTrueMeas(ind(1),1)),sign(vTrueMeas(ind(2),2)),sign(vTrueMeas(ind(3),3))].*val;
% [val,ind]=max(abs(vTrueMeas),[],1);
% 
% V_CHEAT=[V_CHEAT',([sign(vTrueMeas(ind(1),1)),sign(vTrueMeas(ind(2),2)),sign(vTrueMeas(ind(3),3))].*val)'];

% V_CHEAT=mean(diff(yMeas(1:32,:))/10);     
% Vstd=std(diff(yMeas(1:32,:))/10).^2;


% Rstd=mean(sig_XYZ).^2;
% 
% COVInit=diag([sig_XYZ(1,:).^2,Vstd]);
Initial=[yMeas(1,:)';vTrueMeas(1,:)'];
pf = particleFilter(@OrbitalParticleFilterStateFcn,@OrbitalPFMeasurementLikelihoodFcn);
initialize(pf, 5000,[R_Init;V_CHEAT]);

%% Estimate
vidfile = VideoWriter('OrbitalparticlesGoZoom.mp4','MPEG-4');
vidfile.FrameRate=3;
open(vidfile)


xCorrectedPF = zeros(size(xTrueMeas));
COV=zeros(size(xTrueMeas));

%Make a sphere for Earth
[N,E,W]=sphere;

subplot(1,2,1)
plot3(pf.Particles(1,:)*10^(3),pf.Particles(2,:)*10^(3),pf.Particles(3,:)*10^(3),'or','MarkerIndices',floor(linspace(1,length(pf.Particles),10)),'linewidth',0.5)
hold on;
plot3(pf.State(1)*10^(3),pf.State(2)*10^(3),pf.State(3)*10^(3),'*c','linewidth',3)
plot3(yMeas(1,1)*10^(3),yMeas(1,2)*10^(3),yMeas(1,3)*10^(3),'*g','linewidth',3)
plot3(xTrueMeas(1,1)*10^(3),xTrueMeas(1,2)*10^(3),xTrueMeas(1,3)*10^(3),'*m','linewidth',3)
plot3(xTrueOrbit(indMeas(1):indMeas(1)+580,1)*10^3,xTrueOrbit(indMeas(2):indMeas(2)+580,2)*10^3,xTrueOrbit(indMeas(3):indMeas(3)+580,3)*10^3,'c','linewidth',2)
surf(6378*10^3*N,6378*10^3*E,6378*10^3*W)
xlabel('Position [x]')
ylabel('Position [y]')
zlabel('Position [z]')
titlestr=strcat('Iteration=',num2str(0));
title(titlestr)
view(30,15)
hold off;
legend('Particles','Expected state','Y_m_e_a_s','X_true','location','best')
subplot(1,2,2)
plot3(pf.Particles(1,:)*10^(3)-xTrueMeas(1,1)*10^(3),pf.Particles(2,:)*10^(3)-xTrueMeas(1,2)*10^(3),pf.Particles(3,:)*10^(3)-xTrueMeas(1,3)*10^(3),'.r','linewidth',0.5)
hold on;
plot3(pf.State(1)*10^(3)-xTrueMeas(1,1)*10^(3),pf.State(2)*10^(3)-xTrueMeas(1,2)*10^(3),pf.State(3)*10^(3)-xTrueMeas(1,3)*10^(3),'*c','linewidth',3)
plot3(yMeas(1,1)*10^(3)-xTrueMeas(1,1)*10^(3),yMeas(1,2)*10^(3)-xTrueMeas(1,2)*10^(3),yMeas(1,3)*10^(3)-xTrueMeas(1,3)*10^(3),'*g','linewidth',3)
plot3(0,0,0,'*m','linewidth',3);
xlabel('Position [x]')
ylabel('Position [y]')
zlabel('Position [z]')
title('Zoomed in with True Position at the Origin')
% xlim([-2 2])
% ylim([-2 2])

drawnow
hold off;
im=getframe(gcf);
writeVideo(vidfile, im)

timeVector=linspace(0,330,31);
for k=1:31
    % Use measurement y[k] to correct the particles for time k
    [~,tmp] = getStateEstimate(pf);
    tmp=sqrt(diag(tmp));
    xCorrectedPF(k,:) = correct(pf,yMeas(k,:),k); % Filter updates and stores Particles[k|k], Weights[k|k]
    
%      xCorrectedPF(k,:) = correct(pf,X_TrueMeas(k,1:3)); % Filter updates and stores Particles[k|k], Weights[k|k]
    
    % The result is x[k|k]: Estimate of states at time k, utilizing
    % measurements up to time k. This estimate is the mean of all particles
    % because StateEstimationMethod was 'mean'.
    %
    % Now, predict particles at next time step. These are utilized in the
    % next correct command
    [~,tmp] = getStateEstimate(pf);
    COV(k,:)=diag(tmp);
    figure(1);
    subplot(1,2,1)
    plot3(pf.Particles(1,:)*10^(3),pf.Particles(2,:)*10^(3),pf.Particles(3,:)*10^(3),'.r')
    hold on;
    plot3(pf.State(1)*10^(3),pf.State(2)*10^(3),pf.State(3)*10^(3),'*c','linewidth',2)
    %plot Measurement
    plot3(yMeas(k,1)*10^(3),yMeas(k,2)*10^(3),yMeas(k,3)*10^(3),'*g','linewidth',2)
    %plot k true state
    plot3(xTrueMeas(k,1)*10^(3),xTrueMeas(k,2)*10^(3),xTrueMeas(k,3)*10^(3),'*m','linewidth',2)
    %plot orbit
    plot3(xTrueOrbit(indMeas(1):indMeas(1)+580,1)*10^3,xTrueOrbit(indMeas(2):indMeas(2)+580,2)*10^3,xTrueOrbit(indMeas(3):indMeas(3)+580,3)*10^3,'c','linewidth',2)
    
    surf(6378*10^3*N,6378*10^3*E,6378*10^3*W)
%      camtarget(pf.State(1:3)*10^3)\
    legend('Particles','Expected state','Y_m_e_a_s','X_true','location','best')
    view(35,15)
    xlabel('Position X[m]')
    ylabel('Position Y[m]')
    zlabel('Position Z[m]')
    titlestr=strcat('Iteration=',num2str(k));
    title(titlestr)
    hold off
    subplot(1,2,2)
    plot3(pf.Particles(1,:)*10^(3)-xTrueMeas(k,1)*10^(3),pf.Particles(2,:)*10^(3)-xTrueMeas(k,2)*10^(3),pf.Particles(3,:)*10^(3)-xTrueMeas(k,3)*10^(3),'.r','Markerindices',floor(linspace(1,length(pf.Particles),500)),'linewidth',0.5)
    hold on;
    plot3(pf.State(1)*10^(3)-xTrueMeas(k,1)*10^(3),pf.State(2)*10^(3)-xTrueMeas(k,2)*10^(3),pf.State(3)*10^(3)-xTrueMeas(k,3)*10^(3),'*c','linewidth',3)
    plot3(yMeas(k,1)*10^(3)-xTrueMeas(k,1)*10^(3),yMeas(k,2)*10^(3)-xTrueMeas(k,2)*10^(3),yMeas(k,3)*10^(3)-xTrueMeas(k,3)*10^(3),'*g','linewidth',3)
    plot3(0,0,0,'*m','linewidth',3);
    xlabel('Position [x]')
    ylabel('Position [y]')
    zlabel('Position [z]')
    title('Zoomed in with True Position at the Origin')
    % xlim([-2 2])
    % ylim([-2 2])

    drawnow
    hold off;
        im=getframe(gcf);
        writeVideo(vidfile, im);
    predict(pf,DOP(k)); % Filter updates and stores Particles[k+1|k]
end
close(vidfile)
%% Plot
figure();
sgtitle('Position Measurement and predicitions')
subplot(3,1,1);

plot(timeVector,xTrueMeas(1:k,1),timeVector,xCorrectedPF(1:k,1),timeVector,yMeas(1:k,1));
hold on
grid on
plot(timeVector,xCorrectedPF(1:k,1)-3*sqrt(COV(1:k,1)),'-m',timeVector,xCorrectedPF(1:k,1)+3*sqrt(COV(1:k,1)),'-m')
legend('True','Particlte filter estimate','Measured','3 \sigma')
% ylim([-1 1]);
xlim([0 max(timeVector)])
ylabel('X [km]');


subplot(3,1,2);
plot(timeVector,xTrueMeas(1:k,2),timeVector,xCorrectedPF(1:k,2),timeVector,yMeas(1:k,2));
hold on
grid on
plot(timeVector,xCorrectedPF(1:k,2)-3*sqrt(COV(1:k,2)),'-m',timeVector,xCorrectedPF(1:k,2)+3*sqrt(COV(1:k,2)),'-m')
% ylim([-1.5 1.5]);
xlim([0 max(timeVector)])
xlabel('Time [s]');
ylabel('Y [km]');


subplot(3,1,3);
plot(timeVector,xTrueMeas(1:k,3),timeVector,xCorrectedPF(1:k,3),timeVector,yMeas(1:k,3));
hold on
grid on
plot(timeVector,xCorrectedPF(1:k,3)-3*sqrt(COV(1:k,3)),'-m',timeVector,xCorrectedPF(1:k,3)+3*sqrt(COV(1:k,3)),'-m')
% ylim([-1.5 1.5]);
xlim([0 max(timeVector)])
xlabel('Time [s]');
ylabel('Z[km]');

figure();
sgtitle('Velocity Measurement and predicitions')
subplot(3,1,1);
plot(timeVector,xTrueMeas(1:k,4),timeVector,xCorrectedPF(1:k,4));
hold on
grid on
plot(timeVector,xCorrectedPF(1:k,4)-3*sqrt(COV(1:k,4)),'-m',timeVector,xCorrectedPF(1:k,4)+3*sqrt(COV(1:k,4)),'-m')
legend('True','Particlte filter estimate','3 \sigma')
% ylim([-1 1]);
xlim([0 max(timeVector)])
ylabel('$\dot{X}[km/s]$', 'Interpreter','latex');

subplot(3,1,2);
plot(timeVector,xTrueMeas(1:k,5),timeVector,xCorrectedPF(1:k,5));
hold on
grid on
plot(timeVector,xCorrectedPF(1:k,5)-3*sqrt(COV(1:k,5)),'-m',timeVector,xCorrectedPF(1:k,5)+3*sqrt(COV(1:k,5)),'-m')
% ylim([-1.5 1.5]);
xlim([0 max(timeVector)])
xlabel('Time [s]');
ylabel('$\dot{Y}[km/s]$', 'Interpreter','latex');


subplot(3,1,3);
plot(timeVector,xTrueMeas(1:k,6),timeVector,xCorrectedPF(1:k,6));
hold on
grid on
plot(timeVector,xCorrectedPF(1:k,6)-3*sqrt(COV(1:k,6)),'-m',timeVector,xCorrectedPF(1:k,6)+3*sqrt(COV(1:k,6)),'-m')
% ylim([-1.5 1.5]);
xlim([0 max(timeVector)])
xlabel('Time [s]');
ylabel('$\dot{Z}[km/s]$', 'Interpreter','latex');

%%
figure
sgtitle('Position Residuals vs time')
subplot(3,1,1)
plot(timeVector,xTrueMeas(1:k,1)-xCorrectedPF(1:k,1));

ylabel('R_x [km]')
subplot(3,1,2)
plot(timeVector,xTrueMeas(1:k,2)-xCorrectedPF(1:k,2));

ylabel('R_y [km]')
subplot(3,1,3)
plot(timeVector,xTrueMeas(1:k,3)-xCorrectedPF(1:k,3));
xlabel('Time [s]')
ylabel('R_z [km]')



figure
sgtitle('Velocity Residuals vs time')
subplot(3,1,1)
plot(timeVector,xTrueMeas(1:k,4)-xCorrectedPF(1:k,4));

ylabel('$\dot{R_x}$ [km/s]', 'Interpreter','latex')

subplot(3,1,2)
plot(timeVector,xTrueMeas(1:k,5)-xCorrectedPF(1:k,5));

ylabel('$\dot{ R_y }$ [km/s]', 'Interpreter','latex')
subplot(3,1,3)
plot(timeVector,xTrueMeas(1:k,6)-xCorrectedPF(1:k,6));
xlabel('Time [s]')
ylabel('$\dot{R_x}$ [km/s]', 'Interpreter','latex')
%% Orbital propagtion to compare outputs
% [~,xTrue]=ode45(@vdp1,timeVector,[2;0]);
mu=3.986004418*10^14;
Period=2*pi * sqrt((575e3+6378e3)^3/mu);
timeVector=0:10:Period;
input=xCorrectedPF(17,:)'*10^(3);
[~,X_PF]=ode45(@OrbitalODE,timeVector,input);
[~,X_Gibbs]=ode45(@OrbitalODE,timeVector',[R_Gibbs ;V_ECI]*10^(3));
xTrue=xTrueOrbit(17:length(X_PF)+17,:)*10^3;

figure

plot3(X_PF(:,1),X_PF(:,2),X_PF(:,3),'linewidth',2)
hold on
plot3(X_Gibbs(:,1),X_Gibbs(:,2),X_Gibbs(:,3),'linewidth',2)
plot3(xTrue(:,1),xTrue(:,2),xTrue(:,3),'linewidth',2)
surf(6378*10^3*N,6378*10^3*E,6378*10^3*W)
    xlabel('X Position [km]')
    ylabel('Y Position [km]')
    zlabel('Z Position [km]')
legend('Particle Filter','Gibb''s Method','True Orbit')

%% ODE
function sols = OrbitalODE(t,X)
    mu=3.986004418*10^14;
    R=[X(1); X(2); X(3)];
    V=[X(4); X(5); X(6)];
    % Acceleration based on 
    A=-mu*R/(norm(R))^3;
    
    sols=[V;A];

end
function[t,out] = RK4(t,X,h,ODE)
k1=ODE(t,X)';
k2=ODE(t+h/2,X+h*k1/2)';
k3=ODE(t+h/2,X+h*k2/2)';
k4=ODE(t+h,X+h*k3)';
out=X+h*(k1+2*k2+2*k3+k4)/6;
t=t+h;
end
%{
Kyle D. Kemble
University of Colorado Boulder
Dept. Aerospace Engineering Sciences

ASEN-5050: Spaceflight Dynamics

Orbit Determination by Gibbs Method

ODGibbsMethod.m

%}
function [R_ECI V_ECI] = ODGibbsMethod(r_eci)
%%
%{
INPUTS:
r_eci - [km] 3 X 3 matrix where each columnn represents a spacecraft ECI
        Position vector in chronological order.

OUTPUTS:
R_ECI - The middle input position vector
V_ECI - The Estimated Velocity vector at the middle position vector
%}

%% Universal Variables
r_earth = 6378.14; %[km]
u_earth = 398600.45; %[km^3/s^2]
omega_earth = 7.2921158553e-5*180/pi; %[deg/sec]

%% Gibbs Method
%{
r_eci = [0         0        0;
         0        -4464.696 5740.323;
         6378.137 -5102.509 3198.068];
%}         

tol = 1e-2;
%Magnitude of the Radius Vectors
r1 = norm(r_eci(:,1));
r2 = norm(r_eci(:,2));
r3 = norm(r_eci(:,3));

%Cross Product Constants
c12 = cross(r_eci(:,1),r_eci(:,2));
c23 = cross(r_eci(:,2),r_eci(:,3));
c31 = cross(r_eci(:,3),r_eci(:,1));

%Check of Observations are coplaner
ierr = 0;
alpha_coplaner = 90 - acosd(dot(c23,r_eci(:,1))/(norm(c23)*norm(r_eci(:,1))));
if abs(alpha_coplaner) > tol
    ierr = 1;
    disp(strcat('Estimated Position Vectors are outside the ',num2str(tol),' [deg] Angular Tolerance'))
end

%N-Vector
N = r1*c23 + r2*c31 + r3*c12;

%D-Vector
D = c12 + c23 + c31;

%S-Vector
S = r_eci(:,1)*(r2 - r3) + r_eci(:,2)*(r3 - r1) + r_eci(:,3)*(r1 - r2);

%B-Vector
B = cross(D,r_eci(:,2));

L_g = sqrt(u_earth/(norm(N)*norm(D)));

%Velocity Vector at the 2nd Position Indicie
V_ECI = L_g/r2*B + L_g*S;
R_ECI = r_eci(:,2);
end
function [xx,yy,zz] = earth_sphere(varargin)
%EARTH_SPHERE Generate an earth-sized sphere.
%   [X,Y,Z] = EARTH_SPHERE(N) generates three (N+1)-by-(N+1)
%   matrices so that SURFACE(X,Y,Z) produces a sphere equal to 
%   the radius of the earth in kilometers. The continents will be
%   displayed.
%
%   [X,Y,Z] = EARTH_SPHERE uses N = 50.
%
%   EARTH_SPHERE(N) and just EARTH_SPHERE graph the earth as a 
%   SURFACE and do not return anything.
%
%   EARTH_SPHERE(N,'mile') graphs the earth with miles as the unit rather
%   than kilometers. Other valid inputs are 'ft' 'm' 'nm' 'miles' and 'AU'
%   for feet, meters, nautical miles, miles, and astronomical units
%   respectively.
%
%   EARTH_SPHERE(AX,...) plots into AX instead of GCA.
% 
%  Examples: 
%    earth_sphere('nm') produces an earth-sized sphere in nautical miles
%
%    earth_sphere(10,'AU') produces 10 point mesh of the Earth in
%    astronomical units
%
%    h1 = gca;
%    earth_sphere(h1,'mile')
%    hold on
%    plot3(x,y,z)
%      produces the Earth in miles on axis h1 and plots a trajectory from
%      variables x, y, and z
%   Clay M. Thompson 4-24-1991, CBM 8-21-92.
%   Will Campbell, 3-30-2010
%   Copyright 1984-2010 The MathWorks, Inc. 
%% Input Handling
[cax,args,nargs] = axescheck(varargin{:}); % Parse possible Axes input
error(nargchk(0,2,nargs)); % Ensure there are a valid number of inputs
% Handle remaining inputs.
% Should have 0 or 1 string input, 0 or 1 numeric input
j = 0;
k = 0;
n = 50; % default value
units = 'km'; % default value
for i = 1:nargs
    if ischar(args{i})
        units = args{i};
        j = j+1;
    elseif isnumeric(args{i})
        n = args{i};
        k = k+1;
    end
end
if j > 1 || k > 1
    error('Invalid input types')
end
%% Calculations
% Scale factors
Scale = {'km' 'm'  'mile'            'miles'           'nm'              'au'                 'ft';
         1    1000 0.621371192237334 0.621371192237334 0.539956803455724 6.6845871226706e-009 3280.839895};
% Identify which scale to use
try
    myscale = 6378.1363*Scale{2,strcmpi(Scale(1,:),units)};
catch %#ok<*CTCH>
    error('Invalid units requested. Please use m, km, ft, mile, miles, nm, or AU')
end
     
% -pi <= theta <= pi is a row vector.
% -pi/2 <= phi <= pi/2 is a column vector.
theta = (-n:2:n)/n*pi;
phi = (-n:2:n)'/n*pi/2;
cosphi = cos(phi); cosphi(1) = 0; cosphi(n+1) = 0;
sintheta = sin(theta); sintheta(1) = 0; sintheta(n+1) = 0;
x = myscale*cosphi*cos(theta);
y = myscale*cosphi*sintheta;
z = myscale*sin(phi)*ones(1,n+1);
%% Plotting
if nargout == 0
    cax = newplot(cax);
    % Load and define topographic data
    load('topo.mat','topo','topomap1');
    % Rotate data to be consistent with the Earth-Centered-Earth-Fixed
    % coordinate conventions. X axis goes through the prime meridian.
    % http://en.wikipedia.org/wiki/Geodetic_system#Earth_Centred_Earth_Fixed_.28ECEF_or_ECF.29_coordinates
    %
    % Note that if you plot orbit trajectories in the Earth-Centered-
    % Inertial, the orientation of the contintents will be misleading.
    topo2 = [topo(:,181:360) topo(:,1:180)]; %#ok<NODEF>
    
    % Define surface settings
    props.FaceColor= 'texture';
    props.EdgeColor = 'none';
    props.FaceLighting = 'phong';
    props.Cdata = topo2;
    % Create the sphere with Earth topography and adjust colormap
    surface(x,y,z,props,'parent',cax)
    colormap(topomap1)
% Replace the calls to surface and colormap with these lines if you do 
% not want the Earth's topography displayed.
%     surf(x,y,z,'parent',cax)
%     shading flat
%     colormap gray
    
    % Refine figure
    axis equal
    xlabel(['X [' units ']'])
    ylabel(['Y [' units ']'])
    zlabel(['Z [' units ']'])
    view(127.5,30)
else
    xx = x; yy = y; zz = z;
end
end
