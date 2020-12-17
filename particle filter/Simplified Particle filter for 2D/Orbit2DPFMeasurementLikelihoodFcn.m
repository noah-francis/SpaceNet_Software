function likelihood = Orbit2DPFMeasurementLikelihoodFcn(particles,measurement)
c = 299792458 * 10^(-3); %(km/s)
R_E = 6378; %[km]
sig_t = 100 *10^(-9); % s
latS1 = 40 * pi/180; % [rad]
% place sensor unit at 35 deg
latS2 = 35 * pi/180; % [rad]
latS3 = 45 * pi/180; % [rad]
% define the vector locations of sensor units
rS1 = R_E * [ cos(latS1) ; sin(latS1)];
rS2 = R_E * [ cos(latS2) ; sin(latS2)];
rS3 = R_E * [ cos(latS3) ; sin(latS3)];

% Validate the sensor measurement
numberOfMeasurements = 3; % Expected number of measurements
validateattributes(measurement, {'double'}, {'vector', 'numel', numberOfMeasurements}, ...
    'vdpExamplePFMeasurementLikelihoodFcn', 'measurement');
% predict the measurement
predictedMeasurement = particles(1:2,:);
% calculate the distance to each particle
Dist2S1 = vecnorm( predictedMeasurement - rS1,2,1);
Dist2S2 = vecnorm( predictedMeasurement - rS2,2,1);
Dist2S3 = vecnorm( predictedMeasurement- rS3,2,1);

% calculate the time of arrivial
t1 = Dist2S1 / c;
t2 = Dist2S2 / c;
t3 = Dist2S3 / c;

predictedMeasurement = [ t1 - t1; t2 - t1; t3 - t1];
sigma = sig_t * eye(numberOfMeasurements); % variance
numParticles = size(particles,2);
likelihood = zeros(numParticles,1);
C = det(2*pi*sigma) ^ (-0.5);
for kk=1:numParticles
    if norm(particles(1:2,kk)) < R_E 
        likelihood(kk)=0;
    else
        v = -(predictedMeasurement(:,kk)-measurement');
        likelihood(kk) =  C * exp(-0.5 * (v' / sigma * v) );
    end
end

end