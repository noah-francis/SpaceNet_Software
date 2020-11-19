function likelihood = OrbitalPFMeasurementLikelihoodFcn(particles,measurement,k)
% vdpExamplePFMeasurementLikelihoodFcn Example measurement likelihood function
%
% The measurement is the first state.
%
% likelihood = vdpParticleFilterMeasurementLikelihoodFcn(particles, measurement)
%
% Inputs:
%    particles - NumberOfStates-by-NumberOfParticles matrix that holds 
%                the particles
%
% Outputs:
%    likelihood - A vector with NumberOfParticles elements whose n-th
%                 element is the likelihood of the n-th particle
%
% See also extendedKalmanFilter, unscentedKalmanFilter

%   Copyright 2017 The MathWorks, Inc.

%#codegen

% The tag %#codegen must be included if you wish to generate code with 
% MATLAB Coder.
load ECI_Sensor_Locations.mat

c=299792458*10^(-3); %m/s
sig_t=100*10^(-9);
sig_r=1;
numParticles = size(particles,2);



% Validate the sensor measurement
numberOfMeasurements = 3; % Expected number of measurements
validateattributes(measurement, {'double'}, {'vector', 'numel', numberOfMeasurements}, ...
    'vdpExamplePFMeasurementLikelihoodFcn', 'measurement');

% The measurement is first state. Get all measurement hypotheses from particles
% for kk=1:numParticles
% %     [predictedMeasurement(1:3,kk), err(kk,:)] = TDOA(particles(1:3,kk),Sen_ECI(:,:,k)*10^(-3),c,sig_r*10^(-3),sig_t) ;
% end
predictedMeasurement=particles(1:3,:);%+c*sig_t * randn(3,size(particles,2));

% predictedMeasurement = predictedMeasurement./vecnorm(predictedMeasurement,2,1);
% 
% measurement=measurement/norm(measurement);
% err=mean(err);
% Assume the ratio of the error between predicted and actual measurements
% follow a Gaussian distribution with zero mean, variance 0.2
% sigma=cov(predictedMeasurement');
sigma = 5*sqrt([13.1379159113235,0.937401777094900,8.57397351437937]);
sigma=diag(sigma);
mu =0; % mean*
% sigma =0.1 *eye(3); % variance

% Use multivariate Gaussian probability density function, calculate
% likelihood of each particle

likelihood = zeros(numParticles,1);
% C = (2*pi)^(-numberOfMeasurements/2)*det(sigma) ^ (-0.5);
C = det(2 * pi * sigma) ^ (-0.5);

for kk=1:numParticles
    v = -(predictedMeasurement(:,kk)'-measurement);
   likelihood(kk) =  C *exp(-0.5 * (v / sigma * v') );
%     errorRatio = (predictedMeasurement(:,kk)'-measurement)./predictedMeasurement(:,kk)';
%     v = errorRatio-mu;
%     likelihood(kk) =C * exp(-0.5 * (v / sigma * v') );
end
% likelihood=max(likelihood)-(likelihood-max(likelihood));
end

function[XYZ,err]=TDOA(p_T,P,c,sig_r,sig_t)
% This is based on TDoA file from MATLAB file exchange



trials=10;
M=length(P);
in_est_error=200;
% Induce Noise into Sensor Position SEED For repeatablility
rng(40)
err_r=sig_r*randn(size(P));
err_t=sig_t*randn(M,1);
%     A= -2 * [P(:,2)'- P(:,1)', d(2);
%              P(:,3)'- P(:,1)', d(3); 
%              P(:,4)'- P(:,1)', d(4)];
%     b=[ d(2)^2-norm(P(:,2))^2+norm(P(:,1))^2;
%         d(3)^2-norm(P(:,3))^2+norm(P(:,1))^2;
%         d(4)^2-norm(P(:,4))^2+norm(P(:,1))^2];
% %     Ainv=A'*A;
% %     B=A'*b;
%     X=pinv(A)*b;



guess=P(1:3,1)+ P(1:3,1)/norm(P(1:3,1))*575;
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
dummy = repmat(p_T,1,M)-P;
toa = zeros(M,1);   %includes all toa information
for ii = 1:M
    toa(ii) = round(norm(dummy(:,ii))/c,7);    
end
tdoa = toa-toa(1); %tdoa(1)=[];

P=P+err_r;
tdoa = tdoa +err_t;


% err_tMeters=[0;c*err_t];
Pt=(sig_t)^2*eye(4);

    
%%% Taylor Series Expansion Solution
p_T_0 = guess + in_est_error*randn(3,1);    %initial estimate with some error (penalty term)
d = c*tdoa;
f = zeros(M,1);
del_f = zeros(M,3);


    

for ii=1:M
   f(ii)=norm(p_T_0-P(:,ii))-norm(p_T_0-P(:,1)); 
   del_f(ii,1) = (p_T_0(1)-P(1,ii))*norm(p_T_0-P(:,ii))^-1 - (p_T_0(1)-P(1,1))*norm(p_T_0-P(:,1))^-1;
   del_f(ii,2) = (p_T_0(2)-P(2,ii))*norm(p_T_0-P(:,ii))^-1 - (p_T_0(2)-P(2,1))*norm(p_T_0-P(:,1))^-1;
   del_f(ii,3) = (p_T_0(3)-P(3,ii))*norm(p_T_0-P(:,ii))^-1 - (p_T_0(3)-P(3,1))*norm(p_T_0-P(:,1))^-1;    

end
%    del_f(1,:)=(p_T_0 - P(:,1))'*norm(p_T_0 - P(:,1))^(-1) - p_T_0'/norm(p_T_0);
   x_nonlin = pinv(del_f)*(d-f)+p_T_0;
   guess=x_nonlin;
   X=x_nonlin;


%     X(1:3)=X(4)*X(1:3)/norm(X(1:3))+ X(1:3);
    

P=P-err_r;


% 
% rmse(k) = norm(p_T-X)^2;
end

XYZ=X;
err=p_T-X(1:3);
% if norm(err)>1e6
%   XYZ=NaN;
%   err=NaN;
%   DOP=NaN;
% end
end