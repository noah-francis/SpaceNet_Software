function particles = OrbitalParticleFilterStateFcn(particles,DOP) 
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
c=299792458; %m/s
sig_t=420*10^(-9);
% Time-propagate each particle
%
% RK4 intergration technique
dt = 10; % [s] Sample time
for kk=1:numberOfParticles
    particles(:,kk) = rk4orbit(0,particles(:,kk)*10^(3),10)';
end
mean(particles,2)
% Add Gaussian noise with variance 0.025 on each state variable
processNoise=[150*ones(3,1);1*ones(3,1)];
% processNoise=cov(particles'*10^(-3));
% processNoise = sqrt(diag(processNoise));
processNoise = diag(processNoise);
particles = particles*10^(-3) + processNoise * randn(size(particles));
end

function sols = OrbitalODE(t,X)
    mu = 3.986004418*10^14;
    R=[X(1) X(2) X(3)];
    V=[X(4) X(5) X(6)];
    % Acceleration based on 
    A=-mu*R/(norm(R))^3;
    
    sols=[V';A'];

end
function[out] = RK4(t,X,h,ODE)
%Simple Runge-Kutta Method integrator
    k1 = ODE(t,X);
    k2 = ODE(t+h/2,X+h*k1/2);
    k3 = ODE(t+h/2,X+h*k2/2);
    k4 = ODE(t+h,X+h*k3);
    out = X+h/6*(k1+2*k2+2*k3+k4);

end

function X = rk4orbit(t,X,h)
% Numerical Solution, Runge-Kutta 4th Order
    k1 = forbit(t,X);
    k2 = forbit(t+h/2,X+k1*h/2);
    k3 = forbit(t+h/2,X+k2*h/2);
    k4 = forbit(t+h,X+k3*h);
% Step forward in time
X = X+(h/6)*(k1+2*k2+2*k3+k4);
end
function Xdot = forbit(t,X)
    mu = 3.986004418*10^14;    % Gravitational constant (m^3/s^2)
    r = X(1:3);                % Position (m)
    v = X(4:6);                % Velocity (ms^2)
    dr = v;
    dv = (-mu/(norm(r))^3).*r; % Newton's law of gravity
    Xdot = [dr; dv];
    
end