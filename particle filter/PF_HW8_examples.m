%%
clc;clear;close all;
%% MATLAB example PF from website
%10/30/2020
%Keith Poletti implementing
%% Initialize elements
R = 0.1^2;
T = 2*pi/50; % [s] Filter sample time
timeVector = linspace(0,2*pi,50);
[~,xTrue]=ode45(@DuffODE,timeVector,[0;1]);

rng(1); % Fix the random number generator for reproducible results
%Add process Noise to True state?
xTrue = xTrue; %+ 0.1^2 * randn(length(timeVector),2);

%Add Measurement Noise
yTrue = xTrue(:,1);
yMeas = yTrue + sqrt(R) * randn(size(yTrue)); % sqrt(R): Standard deviation of noise

pf = particleFilter(@duffParticleFilterStateFcn,@DuffPFMeasurementLikelihoodFcn);
initialize(pf, 5000,[-2 2;-2 2]);

%% Estimate
% vidfile = VideoWriter('particlesGoZoom.mp4','MPEG-4');
% vidfile.FrameRate=5;
xCorrectedPF = zeros(size(xTrue));
COV=zeros(size(xTrue));
% open(vidfile)


plot(pf.Particles(1,:),pf.Particles(2,:),'.r')
hold on;
plot(pf.State(1),pf.State(2),'*c','linewidth',2)
plot(yMeas(1,1),xTrue(1,2),'*g','linewidth',2)
plot(xTrue(1,1),xTrue(1,2),'*m','linewidth',2)
xlabel('Position [x]')
ylabel('Velocity [dx/dt]')
titlestr=strcat('Iteration=',num2str(0));
title(titlestr)
xlim([-2 2])
ylim([-2 2])
legend('Particles','Expected state','Y_m_e_a_s','X_true')
drawnow
hold off;
% im=getframe(gcf);
% writeVideo(vidfile, im)


for k=1:size(xTrue,1)
    % Use measurement y[k] to correct the particles for time k
    xCorrectedPF(k,:) = correct(pf,yMeas(k)); % Filter updates and stores Particles[k|k], Weights[k|k]
    % The result is x[k|k]: Estimate of states at time k, utilizing
    % measurements up to time k. This estimate is the mean of all particles
    % because StateEstimationMethod was 'mean'.
    %
    % Now, predict particles at next time step. These are utilized in the
    % next correct command
    [~,tmp] = getStateEstimate(pf);
    COV(k,:)=diag(tmp);
    figure(1);

        plot(pf.Particles(1,:),pf.Particles(2,:),'.r')
        hold on;
        plot(pf.State(1),pf.State(2),'*c','linewidth',2)
        plot(yMeas(k,1),xTrue(k,2),'*g','linewidth',2)
        plot(xTrue(k,1),xTrue(k,2),'*m','linewidth',2)
        xlabel('Position [x]')
        ylabel('Velocity [dx/dt]')
        titlestr=strcat('Iteration=',num2str(k));
        title(titlestr)
        xlim([-2 2])
        ylim([-2 2])
        legend('Particles','Expected state','Y_m_e_a_s','X_true')
        drawnow
        hold off;
%         im=getframe(gcf);
%         writeVideo(vidfile, im);
        predict(pf); % Filter updates and stores Particles[k+1|k]
end
% close(vidfile)
%% Plot
figure();
subplot(2,1,1);
plot(timeVector,xTrue(:,1),timeVector,xCorrectedPF(:,1),timeVector,yMeas(:));
hold on
grid on
plot(timeVector,xCorrectedPF(:,1)-3*sqrt(COV(:,1)),'-m',timeVector,xCorrectedPF(:,1)+3*sqrt(COV(:,1)),'-m')
legend('True','Particlte filter estimate','Measured','3 \sigma')
% ylim([-1 1]);
xlim([0 2*pi])
ylabel('x_1');
subplot(2,1,2);
plot(timeVector,xTrue(:,2),timeVector,xCorrectedPF(:,2));
hold on
grid on
plot(timeVector,xCorrectedPF(:,2)-3*sqrt(COV(:,2)),'-m',timeVector,xCorrectedPF(:,2)+3*sqrt(COV(:,2)),'-m')
% ylim([-1.5 1.5]);
xlim([0 2*pi])
xlabel('Time [s]');
ylabel('dx_1/dt');



%% ODE
function[out]=DuffODE(t,X)

k=1;
eta=1000;
x=X(1);
dx=X(2);
d2x=-k*x-eta*x^3;

out=[dx;d2x];

end
function[t,out] = RK4(t,X,h,ODE)
k1=ODE(t,X)';
k2=ODE(t+h/2,X+h*k1/2)';
k3=ODE(t+h/2,X+h*k2/2)';
k4=ODE(t+h,X+h*k3)';
out=X+h*(k1+2*k2+2*k3+k4)/6;
t=t+h;
end