%{
Kyle D. Kemble
University of Colorado Boulder
Dept. Aerospace Engineering Sciences

ASEN-5050: Spaceflight Dynamics

Convert a State Vector (r, v) to Orbital Elements (a,e,i,RAAN,w,v)

Reference: Orbital Mechanics for Engineering Students by H. Curtis

StateVectorToOrbitElements.m

%}
function OrbitElements = StateVectorToOrbitElements(R,V)
%{
INPUTS:
R: Radius Vector in ECI [km]
V: Velocity Vector in ECI [km/s]

OUTPUTS:
Struct called "OrbitElements" with the format:
    'a'     Semimajor Axis [km]
    'e'     Eccentricity Magnitude
    'i'     Inclination [deg]
    'w'     Argument of Periapse [deg]
    'RAAN'  Right Ascension of the Ascending Node [deg]
    'v'     True Anomaly [deg]
    'Alt'   Altitudes of Apogee and Perigee Respectively [km]
    'FPA'	Flight Path Angle [deg]
    'OE'    Orbit Energy [MJ/kg]
%}
%% Universal Variables
r_earth = 6378.137; %km
u_earth = 398600; %km^3/s^2

%Inertial Axes
I = [1; 0; 0];
J = [0; 1; 0];
K = [0; 0; 1];

%% Determine Angular Momentum
H = cross(R,V); %km^2/s

%Orbit Energy
OE = norm(V)^2/2 - u_earth/norm(R); %km^2/s^2/kg => MJ/kg

%% Determination of Eccentricity
%Eccentricity Vector
e_vector = cross(V,H)/u_earth-R/norm(R);

%Scalor Eccentricity
e = norm(e_vector);

%% Semimajor Axis 
a = -u_earth/(norm(V)^2-(2*u_earth)/norm(R)); %km

%% Determine Inclination
i = acos(dot(K,H)/(norm(K)*norm(H))); %rad

%% Determine the Right Ascension of the Ascending Node
%Node Line Vector
N = cross(K,H);

%RAAN
if N(2) < 0
    RAAN = 2*pi - acos(dot(I,N)/(norm(I)*norm(N)));
else
    RAAN = acos(dot(I,N)/(norm(I)*norm(N)));
end

%% Determine The Argument of Periapse
if e_vector(3) < 0
    w = 2*pi - acos(dot(N,e_vector)/(norm(N)*norm(e_vector)));
else
    w = acos(dot(N,e_vector)/(norm(N)*norm(e_vector)));
end

%% Determine True Anomaly
if dot(R,V) < 0
    v = 2*pi - acos(dot(e_vector,R)/(norm(e_vector)*norm(R)));
else
    v = acos(dot(e_vector,R)/(norm(e_vector)*norm(R)));
end

%% Apogee & Perigee Altitudes
%Apogee
r_a = a*(1-e^2)/(1-e); %km
    h_a = r_a - r_earth; %km

%Perigee
r_p = a*(1-e^2)/(1+e); %km
    h_p = r_p - r_earth; %km

%% Flight Path Angle
phi_fpa = atan2(e*sin(v)/sqrt(1+2*e*cos(v)+e^2),...
               (1+e*cos(v))/sqrt(1+2*e*cos(v)+e^2));

%% Format Outputs
OrbitElements = struct(...
    'a',   a,...
    'e',   e,...
    'i',   i*180/pi,...
    'w',   w*180/pi,...
    'RAAN',RAAN*180/pi,...
    'v',   v*180/pi,...
    'Alt', [h_a,h_p],...
    'FPA', phi_fpa*180/pi,...
    'OE',  OE);