%%
clc;clear;close all;
%% MATLAB example PF from website
%10/30/2020
%Keith Poletti implementing
%https://www.mathworks.com/help/control/ug/nonlinear-state-estimation-using-unscented-kalman-filter.html
%% Initialize elements
R = 0.2;
T = 0.05; % [s] Filter sample time
timeVector = 0:T:5;
[~,xTrue]=ode45(@vdp1,timeVector,[2;0]);
rng(1); % Fix the random number generator for reproducible results
yTrue = xTrue(:,1);
yMeas = yTrue .* (1+sqrt(R)*randn(size(yTrue))); % sqrt(R): Standard deviation of noise

pf = particleFilter(@vdpParticleFilterStateFcn,@vdpExamplePFMeasurementLikelihoodFcn);
initialize(pf, 1000, [2;0], 0.01*eye(2));

%% Estimate
xCorrectedPF = zeros(size(xTrue));
for k=1:size(xTrue,1)
    % Use measurement y[k] to correct the particles for time k
    xCorrectedPF(k,:) = correct(pf,yMeas(k)); % Filter updates and stores Particles[k|k], Weights[k|k]
    % The result is x[k|k]: Estimate of states at time k, utilizing
    % measurements up to time k. This estimate is the mean of all particles
    % because StateEstimationMethod was 'mean'.
    %
    % Now, predict particles at next time step. These are utilized in the
    % next correct command
     predict(pf); % Filter updates and stores Particles[k+1|k]

    
    [~,tmp] = getStateEstimate(pf);
    COV(k,:)=diag(tmp);
        figure(1);

        plot(pf.Particles(1,:),pf.Particles(2,:),'.r')
        hold on;
        plot(pf.State(1),pf.State(2),'*c')
        plot(yMeas(k,1),xTrue(k,2),'*g')
        plot(xTrue(k,1),xTrue(k,2),'*m')
        xlabel('Position [x]')
        ylabel('Velocity [dx/dt]')
        titlestr=strcat('Iteration=',num2str(timeVector(end)/(timeVector(2)-timeVector(1))));
        title(titlestr)
        xlim([-5 5])
        ylim([-5 5])
        legend('Particles','Expected state','Y_m_e_a_s','X_true')
        drawnow
        hold off;
end
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
figure();
subplot(2,1,1);
plot(timeVector,-3*sqrt(COV(:,1)),'-m',timeVector,3*sqrt(COV(:,1)),'-m',timeVector,yMeas(:)-xTrue(:,1))
legend('True','Particlte filter estimate','Measured','3 \sigma')
% ylim([-1 1]);
xlim([0 2*pi])
ylabel('x_1');
subplot(2,1,2);
plot(timeVector,-3*sqrt(COV(:,2)),'-m',timeVector,3*sqrt(COV(:,2)),'-m',timeVector,yMeas(:)-xCorrectedPF(:,2))
legend('True','Particlte filter estimate','Measured','3 \sigma')
% ylim([-1 1]);
xlim([0 2*pi])
ylabel('dx_1/dt');

