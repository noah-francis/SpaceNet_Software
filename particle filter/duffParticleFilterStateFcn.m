function particles = duffParticleFilterStateFcn(particles) 
% vdpParticleFilterStateFcn Example state transition function for particle
%                           filter
%
% Discrete-time approximation to van der Pol ODEs for mu = 1. 
% Sample time is 0.05s.
%
% predictedParticles = vdpParticleFilterStateFcn(particles)
%
% Inputs:
%    particles - Particles at current time. Matrix with dimensions
%                [NumberOfStates NumberOfParticles] matrix
%
% Outputs:
%    predictedParticles - Predicted particles for the next time step
%
% See also particleFilter

%   Copyright 2017 The MathWorks, Inc.

%#codegen

% The tag %#codegen must be included if you wish to generate code with 
% MATLAB Coder.

[numberOfStates, numberOfParticles] = size(particles);
    
% Time-propagate each particle
%
% Euler integration of continuous-time dynamics x'=f(x) with sample time dt
dt=2*pi/50; % [s] Sample time
for kk=1:numberOfParticles
    particles(:,kk) = RK4(0,particles(:,kk),dt,@DuffODE)';
end

% Add Gaussian noise with variance 0.1 on each state variable
processNoise = 0.1^2*eye(numberOfStates);
particles = particles + processNoise * randn(size(particles));
end

function[out]=DuffODE(t,X)

k=1;
eta=1000;
x=X(1);
dx=X(2);
d2x=-k*x-eta*x^3;

out=[dx;d2x];

end
function[out] = RK4(t,X,h,ODE)
k1=ODE(t,X);
k2=ODE(t+h/2,X+h*k1/2);
k3=ODE(t+h/2,X+h*k2/2);
k4=ODE(t+h,X+h*k3);
out=X+h*(k1+2*k2+2*k3+k4)/6;

end