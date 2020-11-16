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
v_TrueMeas = ECI_True(:,11:13)*10^(-3);%True Velocity when measurement happens
X_TrueMeas = [R_TrueMeas, v_TrueMeas];
indMeas = ECI_True(:,end); % Indices in the entire orbit state where a measurement happens


%Load in the measurements
yMeas=ECI_Tdoa(:,7:9)*10^(-3);
k=find(yMeas(:,1)>10^4);
% yMeas(k,:)=(yMeas(k-1,:)+yMeas(k+1,:))/2;
sig_r=ECI_Tdoa(:,1);%*10^3;


rng(1); % Fix the random number generator for reproducible results
%Add process Noise to True state?
% xTrue = xTrue; %+ 0.1^2 * randn(length(timeVector),2);

%Add Measurement Noise
% yTrue = xTrue(:,1);
% yMeas = yTrue + sqrt(R) * randn(size(yTrue)); % sqrt(R): Standard deviation of noise


%Vector to initialize the particle Filter with in 3 sigma bounds

R_Init= [(yMeas(1,:) - 420)', ...
         (yMeas(1,:) + 420)'];
     
     
%CHEATERS way of admissble region DO NOT USE IN FINAL VERSION

[~ ,V_ECI] = ODGibbsMethod([yMeas(1,:)',yMeas(16,:)',yMeas(32,:)'])
[~ ,V_ECIerr] = ODGibbsMethod([yMeas(1,:)',yMeas(5,:)',yMeas(10,:)'])
V_CHEAT=[V_ECI-2*V_ECIerr,V_ECI+2*V_ECIerr];
[val,ind]=min(abs(v_TrueMeas),[],1);
V_CHEAT=[sign(v_TrueMeas(ind(1),1)),sign(v_TrueMeas(ind(2),2)),sign(v_TrueMeas(ind(3),3))].*val;
[val,ind]=max(abs(v_TrueMeas),[],1);

V_CHEAT=[V_CHEAT',([sign(v_TrueMeas(ind(1),1)),sign(v_TrueMeas(ind(2),2)),sign(v_TrueMeas(ind(3),3))].*val)'];



Initial=[R_Init;V_CHEAT];
pf = particleFilter(@OrbitalParticleFilterStateFcn,@OrbitalPFMeasurementLikelihoodFcn);
initialize(pf, 5000,Initial);

%% Estimate
vidfile = VideoWriter('OrbitalparticlesGoZoom.mp4','MPEG-4');
vidfile.FrameRate=3;
xCorrectedPF = zeros(size(X_TrueMeas));
COV=zeros(size(X_TrueMeas));
open(vidfile)


%Make a sphere for Eart
[N,E,W]=sphere;

subplot(1,2,1)
plot3(pf.Particles(1,:)*10^(3),pf.Particles(2,:)*10^(3),pf.Particles(3,:)*10^(3),'or','MarkerIndices',floor(linspace(1,length(pf.Particles),10)),'linewidth',0.5)
hold on;
plot3(pf.State(1)*10^(3),pf.State(2)*10^(3),pf.State(3)*10^(3),'*c','linewidth',3)
plot3(yMeas(1,1)*10^(3),yMeas(1,2)*10^(3),yMeas(1,3)*10^(3),'*g','linewidth',3)
plot3(X_TrueMeas(1,1)*10^(3),X_TrueMeas(1,2)*10^(3),X_TrueMeas(1,3)*10^(3),'*m','linewidth',3)
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
plot3(pf.Particles(1,:)*10^(3)-X_TrueMeas(1,1)*10^(3),pf.Particles(2,:)*10^(3)-X_TrueMeas(1,2)*10^(3),pf.Particles(3,:)*10^(3)-X_TrueMeas(1,3)*10^(3),'.r','linewidth',0.5)
hold on;
plot3(pf.State(1)*10^(3)-X_TrueMeas(1,1)*10^(3),pf.State(2)*10^(3)-X_TrueMeas(1,2)*10^(3),pf.State(3)*10^(3)-X_TrueMeas(1,3)*10^(3),'*c','linewidth',3)
plot3(yMeas(1,1)*10^(3)-X_TrueMeas(1,1)*10^(3),yMeas(1,2)*10^(3)-X_TrueMeas(1,2)*10^(3),yMeas(1,3)*10^(3)-X_TrueMeas(1,3)*10^(3),'*g','linewidth',3)
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

timeVector=linspace(0,300,30);
for k=1:30
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
    plot3(X_TrueMeas(k,1)*10^(3),X_TrueMeas(k,2)*10^(3),X_TrueMeas(k,3)*10^(3),'*m','linewidth',2)
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
    subplot(1,2,2)
    plot3(pf.Particles(1,:)*10^(3)-X_TrueMeas(k,1)*10^(3),pf.Particles(2,:)*10^(3)-X_TrueMeas(k,2)*10^(3),pf.Particles(3,:)*10^(3)-X_TrueMeas(k,3)*10^(3),'.r','Markerindices',floor(linspace(1,length(pf.Particles),500)),'linewidth',0.5)
    hold on;
    plot3(pf.State(1)*10^(3)-X_TrueMeas(k,1)*10^(3),pf.State(2)*10^(3)-X_TrueMeas(k,2)*10^(3),pf.State(3)*10^(3)-X_TrueMeas(k,3)*10^(3),'*c','linewidth',3)
    plot3(yMeas(k,1)*10^(3)-X_TrueMeas(k,1)*10^(3),yMeas(k,2)*10^(3)-X_TrueMeas(k,2)*10^(3),yMeas(k,3)*10^(3)-X_TrueMeas(k,3)*10^(3),'*g','linewidth',3)
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
subplot(3,1,1);
plot(timeVector,X_TrueMeas(1:k,1),timeVector,xCorrectedPF(1:k,1));
hold on
grid on
plot(timeVector,xCorrectedPF(1:k,1)-3*sqrt(COV(1:k,1)),'-m',timeVector,xCorrectedPF(1:k,1)+3*sqrt(COV(1:k,1)),'-m')
legend('True','Particlte filter estimate','Measured','3 \sigma')
% ylim([-1 1]);
xlim([170 max(timeVector)])
ylabel('X [km]');


subplot(3,1,2);
plot(timeVector,X_TrueMeas(1:k,2),timeVector,xCorrectedPF(1:k,2),timeVector,yMeas(1:k,2));
hold on
grid on
plot(timeVector,xCorrectedPF(1:k,2)-3*sqrt(COV(1:k,2)),'-m',timeVector,xCorrectedPF(1:k,2)+3*sqrt(COV(1:k,2)),'-m')
% ylim([-1.5 1.5]);
xlim([170 max(timeVector)])
xlabel('Time [s]');
ylabel('Y [km]');


subplot(3,1,3);
plot(timeVector,X_TrueMeas(1:k,3),timeVector,xCorrectedPF(1:k,3),timeVector,yMeas(1:k,3));
hold on
grid on
plot(timeVector,xCorrectedPF(1:k,3)-3*sqrt(COV(1:k,3)),'-m',timeVector,xCorrectedPF(1:k,3)+3*sqrt(COV(1:k,3)),'-m')
% ylim([-1.5 1.5]);
xlim([0 max(timeVector)])
xlabel('Time [s]');
ylabel('Z[km]');


%%
figure
plot(timeVector,X_TrueMeas(1:k,1)-xCorrectedPF(1:k,1));


%% ODE
function sols = OrbitalODE(t,X)
    mu=398600;
    R=[X(1) X(2) X(3)];
    V=[X(4) X(5) X(6)];
    % Acceleration based on 
    A=-mu*R/(norm(R))^3;
    
    sols=[V,A];

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
