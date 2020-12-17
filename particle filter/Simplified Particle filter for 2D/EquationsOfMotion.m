%% Housekeeping
clc;clear;close all;
%% Constants
R_E = 6378; %[km]
R_sat = R_E + 575; %[km]
mu = 398600; %[km^3/s^2]
c = 299792458 * 10^(-3); %(km/s)
%Line of sight max distance for a unit to see the satellite
LOS = 1500; %[km]
%period of satellite
T_sat = 2 * pi * R_sat ^( 3/2)/ sqrt(mu); % [s]

%mean motion
n = 2 * pi / T_sat; %[rad/s]

% Magnitude of velocity
V_sat = n * R_sat;

% place sensor unit at latitude of 40 deg for Boulder

latS1 = 40 * pi/180; % [rad]

% place sensor unit at 35 deg
latS2 = 35 * pi/180; % [rad]

latS3 = 45 * pi/180; % [rad]

% define the vector locations of sensor units
rS1 = R_E * [ cos(latS1) , sin(latS1)];
rS2 = R_E * [ cos(latS2) , sin(latS2)];
rS3 = R_E * [ cos(latS3) , sin(latS3)];

%% Earth circle
X_earth = R_E*cos(linspace(0,2*pi,1000));
Y_earth = R_E*sin(linspace(0,2*pi,1000));
%% Propagate to check orbit
phi_0 = 0; 

% define intital satellite vector
rSatVec = R_sat * [ cos(phi_0) , sin(phi_0)];
vSatVec = V_sat * [ -sin(phi_0) , cos(phi_0)];

X0 = [ rSatVec , vSatVec];
% tspan = [0 , T_sat];
% [ t , y] = ode45(@OrbitalODE2D, tspan, X0);

tvec=0:10:T_sat;
l = length(tvec);
X = zeros(l,4);
X(1,:) = X0;

for i = 1:l-1
    X(i+1,:) = RK4(tvec(i),X(i,:),10,@OrbitalODE2D);  
end

%% Find where the Sat is in Sight of both sensors
% distance to each unit from satellite
Dist2S1 = vecnorm(X(:,1:2)- rS1,2,2);
Dist2S2 = vecnorm(X(:,1:2)- rS2,2,2);
Dist2S3 = vecnorm(X(:,1:2)- rS3,2,2);
% find the indices that are below LOS in distance to both units

ind = find(Dist2S1 <= LOS & Dist2S2 <=LOS & Dist2S3 <=LOS);
Lock = X(ind,1:2);
%% Plots

figure; 
plot(X(:,1),X(:,2))
hold on
grid on
plot( X_earth, Y_earth,'b','linewidth',2)
plot( Lock(:,1) , Lock(:,2),'om')
plot(rS1(1),rS1(2),'*r')
plot(rS2(1),rS2(2),'*r')
plot(rS3(1),rS3(2),'*r')
xlim([-7500, 7500])
ylim([-7500, 7500])
xlabel( 'X [km]')
ylabel('Y [km]')

%%
tdoa12 = ( vecnorm( rS2 - Lock ,2,2) - vecnorm( rS1 - Lock, 2,2) )/c;
tdoa13 = ( vecnorm( rS3 - Lock, 2,2) - vecnorm( rS1 - Lock, 2,2) )/c;
% Define Tdoa with sensor 1 as the refernce
tdoa1 = [ zeros(length(ind),1), tdoa12, tdoa13];

tdoa21 = ( vecnorm( rS1 - Lock ,2,2) - vecnorm( rS2 - Lock, 2,2) )/c;
tdoa23 = ( vecnorm( rS3 - Lock, 2,2) - vecnorm( rS2 - Lock, 2,2) )/c;
% Define Tdoa with sensor 2 as the refernce
tdoa2 = [ tdoa21, zeros(length(ind),1), tdoa23];

tdoa31 = ( vecnorm( rS1 - Lock ,2,2) - vecnorm( rS3 - Lock, 2,2) )/c;
tdoa32 = ( vecnorm( rS2 - Lock, 2,2) - vecnorm( rS3 - Lock, 2,2) )/c;
% Define Tdoa with sensor 3 as the refernce
tdoa3 = [ tdoa31, tdoa32, zeros(length(ind),1)];



%% Functions
function sols = OrbitalODE2D(t,X)
    mu = 398600;
    R=[X(1) ; X(2)];
    V=[X(3) ; X(4)];
    % Accelertion based on 
    A=(-mu/(norm(R))^3) * R;
    
    sols=[V ; A ];

end

function[out] = RK4(t,X,h,ODE)
%Simple Runge-Kutta Method integrator
    k1 = ODE(t,X)';
    k2 = ODE(t+h/2,X+h*k1/2)';
    k3 = ODE(t+h/2,X+h*k2/2)';
    k4 = ODE(t+h,X+h*k3)';
    out = X+h/6*(k1+2*k2+2*k3+k4);

end
