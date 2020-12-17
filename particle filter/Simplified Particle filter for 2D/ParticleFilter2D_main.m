%% Housekeeping 
clc;clear; close all;
%% 2D Particle filter for orbit determination
%{
By Keith Poletti 12/16/2020
Simplified example of a particle filter on time delay of arrival
measurements 
%}
%% Constants
R_E = 6378; %[km]
R_sat = R_E + 575; %[km]
mu = 398600; %[km^3/s^2]
c = 299792458 * 10^(-3); %(km/s)
sig_t = 100 *10^(-9); % s
%Line of sight max distance for a unit to see the satellite
LOS = 2000; %[km]
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
%% Create truth data
phi_0 = 0; 

% define intital satellite vector
rSatVec = R_sat * [ cos(phi_0) , sin(phi_0)];
vSatVec = V_sat * [ -sin(phi_0) , cos(phi_0)];

X0 = [ rSatVec , vSatVec];
dt = 10;
tvec=0:dt:T_sat;
l = length(tvec);
X = zeros(l,4);
X(1,:) = X0;

for i = 1:l-1
    X(i+1,:) = RK4(tvec(i),X(i,:),dt,@OrbitalODE2D);  
end

%% Find where the Sat is in Sight of both sensors
% distance to each unit from satellite
Dist2S1 = vecnorm(X(:,1:2)- rS1,2,2);
Dist2S2 = vecnorm(X(:,1:2)- rS2,2,2);
Dist2S3 = vecnorm(X(:,1:2)- rS3,2,2);
% find the indices that are below LOS in distance to both units

ind = find(Dist2S1 <= LOS & Dist2S2 <=LOS & Dist2S3 <=LOS);
Lock = X(ind,1:2);

rng(1);
t1 = Dist2S1(ind) / c + sig_t * randn(length(ind),1);
t2 = Dist2S2(ind) / c + sig_t * randn(length(ind),1);
t3 = Dist2S3(ind) / c + sig_t * randn(length(ind),1);

tdoa1 = [ t1 - t1, t2 - t1, t3 - t1];
tdoa2 = [ t1 - t2, t2 - t2, t3 - t2];
tdoa3 = [ t1 - t3, t2 - t3, t3 - t3];

yMeas = tdoa1;

%% initial Guess
% calculate the initial position with TDoA
P=[rS1',rS2',rS3'];
XYZ = TDOA_calc(Lock(1,:),P,c,5,tdoa1(1,:));
phi_0 = atan2(XYZ(2) , XYZ(1));

%%
xTrue=X(:,1:2);
xCorrectedPF = zeros(length(yMeas),4);
COV=zeros(length(yMeas),4);

pf = particleFilter(@Orbit2DStateFcn,@Orbit2DPFMeasurementLikelihoodFcn);

%  Initialized based on general area of particles and velocity
% initialize( pf, 5000, [3100 6500; 2100 6000; -6 -3; 3 7.0])

%Initialized based on very specific areas and velocities

% initialize( pf, 5000, [5326.8 6026.8; 3693.8 3893.8; -5.1312 -3.1312; 5.345099 7.345099])

% Initialized based on TDoA guess for position and velocity based on
% physics 
initialize( pf, 10000, [XYZ ;-V_sat * sin(phi_0) ; V_sat * cos(phi_0)], diag([100^2*ones(1,2),0.5^2*ones(1,2)]))


% pf.ResamplingPolicy.MinEffectiveParticleRatio=0.25;
% 

plot(pf.Particles(1,:),pf.Particles(2,:),'.r')
hold on;
plot(pf.State(1),pf.State(2),'*c','linewidth',2)
% plot(yMeas(1,1),xTrue(1,2),'*g','linewidth',2)
plot(xTrue(ind(1),1),xTrue(ind(1),2),'*m','linewidth',2)

plot( X_earth, Y_earth,'b','linewidth',2)
plot(rS1(1),rS1(2),'*r','linewidth',2)
plot(rS2(1),rS2(2),'*r','linewidth',2)
plot(rS3(1),rS3(2),'*r','linewidth',2)
xlim([3100 6600])
ylim([2100 6200])
xlabel('Position [x]')
ylabel('Velocity [dx/dt]')
titlestr=strcat('Iteration=',num2str(0));
title(titlestr)

legend('Particles','Expected state','X_t_r_u_e')
drawnow
hold off;


for k = 1:length(yMeas)
   % Correct PF
   xCorrectedPF(k,:) = correct(pf,yMeas(k,:)); 
   [~,tmp] = getStateEstimate(pf);
   COV(k,:)=diag(tmp);
    figure(1);

        plot(pf.Particles(1,:),pf.Particles(2,:),'.r')
        hold on;
        plot(pf.State(1),pf.State(2),'*c','linewidth',2)
        plot(xTrue(ind(k),1),xTrue(ind(k),2),'*m','linewidth',2)
        plot( X_earth, Y_earth,'b','linewidth',2)
        plot(rS1(1),rS1(2),'*r','linewidth',2)
        plot(rS2(1),rS2(2),'*r','linewidth',2)
        plot(rS3(1),rS3(2),'*r','linewidth',2)
        xlabel('X Position [km]')
        ylabel('Y Position [km]')
        titlestr=strcat('Iteration=',num2str(k));
        title(titlestr)
        xlim([3100 6600])
        ylim([2100 6200])
        legend('Particles','Expected state','X_t_r_u_e')
        drawnow
        hold off;
%         im=getframe(gcf);
%         writeVideo(vidfile, im);
        predict(pf); % Filter updates and stores Particles[k+1|k]
end

%% Plots
timeVector = tvec(ind);

figure();
sgtitle('Position Measurement and predicitions')
subplot(2,1,1);

plot(timeVector,xTrue(ind,1),timeVector,xCorrectedPF(1:k,1));
hold on
grid on
plot(timeVector,xCorrectedPF(1:k,1)-3*sqrt(COV(1:k,1)),'-m',timeVector,xCorrectedPF(1:k,1)+3*sqrt(COV(1:k,1)),'-m')
legend('True','Particlte filter estimate','3 \sigma')
% ylim([-1 1]);
xlim([timeVector(1) max(timeVector)])
ylabel('X [km]');


subplot(2,1,2);
plot(timeVector,xTrue(ind,2),timeVector,xCorrectedPF(1:k,2));
hold on
grid on
plot(timeVector,xCorrectedPF(1:k,2)-3*sqrt(COV(1:k,2)),'-m',timeVector,xCorrectedPF(1:k,2)+3*sqrt(COV(1:k,2)),'-m')
% ylim([-1.5 1.5]);
xlim([timeVector(1) max(timeVector)])
xlabel('Time [s]');
ylabel('Y [km]');
%%
figure
sgtitle('Position Residuals vs time')
subplot(2,1,1)
plot(timeVector,xTrue(ind,1)-xCorrectedPF(1:k,1));
xlim([timeVector(1) max(timeVector)])
grid on

ylabel('R_x [km]')
subplot(2,1,2)
plot(timeVector,xTrue(ind,2)-xCorrectedPF(1:k,2));
grid on
xlim([timeVector(1) max(timeVector)])
xlabel('Time [s]');
ylabel('R_y [km]')
%%
figure();
sgtitle('Velocity Measurement and predicitions')
subplot(2,1,1);
plot(timeVector,X(ind,3),timeVector,xCorrectedPF(1:k,3));
hold on
grid on
plot(timeVector,xCorrectedPF(1:k,3)-3*sqrt(COV(1:k,3)),'-m',timeVector,xCorrectedPF(1:k,3)+3*sqrt(COV(1:k,3)),'-m')
% ylim([-1.5 1.5]);
xlim([timeVector(1) max(timeVector)])
xlabel('Time [s]');
ylabel('V_X [km]');
subplot(2,1,2);
plot(timeVector,X(ind,4),timeVector,xCorrectedPF(1:k,4));
hold on
grid on
plot(timeVector,xCorrectedPF(1:k,4)-3*sqrt(COV(1:k,4)),'-m',timeVector,xCorrectedPF(1:k,4)+3*sqrt(COV(1:k,4)),'-m')
% ylim([-1.5 1.5]);
xlim([timeVector(1) max(timeVector)])
xlabel('Time [s]');
ylabel('V_X [km]');

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
function[XYZ,err]=TDOA_calc(p_T,P,c,err_r,tdoa)
% This is based on TDoA file from MATLAB file exchange

method=1;


trials=100;
M=length(P);
in_est_error=0;
% Induce Noise into Sensor Position SEED For repeatablility
% rng(40)



% makes an initial guess for this case 575 km above LASP 
guess=P(1:2,1)+ P(1:2,1)/norm(P(1:2,1))*575;
for k=1:trials
%     err_r=sig_r*randn(size(P));
%     err_t=sig_t*randn(M-1,1);

%sensor positon vectors
% P = zeros(3,M); %includes all sensor positon vectors
% for ii=1:M
%     P(:,ii)=range_s*2*(rand(3,1)-0.5);
% end
%  
% p_T = range_T*2*(rand(3,1)-0.5);   %target positon vector
% %finding TOAs 

Perr=P+err_r;
% tdoa = tdoa +err_t;


% err_tMeters=[0;c*err_t];


    
%%% Taylor Series Expansion Solution
p_T_0 = guess + in_est_error*randn(size(guess));    %initial estimate with some error (penalty term)
d = c*tdoa';
f = zeros(M,1);
del_f = zeros(M,2);

if method==1
    


for ii=1:M
   f(ii)=norm(p_T_0-Perr(:,ii))-norm(p_T_0-Perr(:,1)); 
   del_f(ii,1) = (p_T_0(1)-Perr(1,ii))*norm(p_T_0-Perr(:,ii))^-1 - (p_T_0(1)-Perr(1,1))*norm(p_T_0-Perr(:,1))^-1;
   del_f(ii,2) = (p_T_0(2)-Perr(2,ii))*norm(p_T_0-Perr(:,ii))^-1 - (p_T_0(2)-Perr(2,1))*norm(p_T_0-Perr(:,1))^-1;
%    del_f(ii,3) = (p_T_0(3)-Perr(3,ii))*norm(p_T_0-Perr(:,ii))^-1 - (p_T_0(3)-Perr(3,1))*norm(p_T_0-Perr(:,1))^-1;    

end
%    del_f(1,:)=(p_T_0 - P(:,1))'*norm(p_T_0 - P(:,1))^(-1) - p_T_0'/norm(p_T_0);
   x_nonlin = pinv(del_f)*(d-f)+p_T_0;
   guess=x_nonlin;
   X=x_nonlin;
else
    A= -2 * [P(:,2)'- P(:,1)', d(2);
             P(:,3)'- P(:,1)', d(3); 
             P(:,4)'- P(:,1)', d(4)];
    b=[ d(2)^2-norm(P(:,2))^2+norm(P(:,1))^2;
        d(3)^2-norm(P(:,3))^2+norm(P(:,1))^2;
        d(4)^2-norm(P(:,4))^2+norm(P(:,1))^2];
%     Ainv=A'*A;
%     B=A'*b;
    X=pinv(A)*b;
%     X(1:3)=X(4)*X(1:3)/norm(X(1:3))+ X(1:3);
    
end
% P=P-err_r;


% 
% rmse(k) = norm(p_T-X)^2;
end



XYZ=X;
% err=p_T-X(1:3);
% if norm(err)>1e6
%   XYZ=NaN;
%   err=NaN;
%   DOP=NaN;
% end
end