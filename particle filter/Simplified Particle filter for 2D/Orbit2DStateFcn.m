function particles = Orbit2DStateFcn(particles) 


[numberOfStates, numberOfParticles] = size(particles);

dt=10; % [s] Sample time
for kk=1:numberOfParticles
    particles(:,kk) = RK4(0,particles(:,kk),dt,@OrbitalODE2D);
end
% Add Gaussian noise with variance 0.1 on each state variable
processNoise = 0.15*eye(numberOfStates);
particles = particles + processNoise * randn( size(particles));
end
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
    k1 = ODE(t,X);
    k2 = ODE(t+h/2,X+h*k1/2);
    k3 = ODE(t+h/2,X+h*k2/2);
    k4 = ODE(t+h,X+h*k3);
    out = X+h/6*(k1+2*k2+2*k3+k4);

end