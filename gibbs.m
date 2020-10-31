function [v_vec, coplanar_val] = gibbs(r_vec, r_mag, mu)
%   Use Gibb's Method to produce a velocity vector from three consecutive
%   position vectors.
%   -----------------------------------------------------------------------
%
%   Inputs:
%   --------
%   'r_vec' -- array of three time-consecutive position vectors [km]
%   'r_mag' -- array of corresponding vector magnitides [km]
%   'mu'    -- Standard Gravitational Parameter of the body orbited 
%               [km^3 / s^2]
%
%   Returns:
%   --------
%   'v_vec'        -- velocity vector at the second time [km/s]
%   'coplanar_val' -- dot product of unit vectors that should be coplanar 
%                       for Gibb's Method to work properly [~] (closer to 
%                       zero is better)

% compute unit vectors (two are unused)
% C_12 = cross(r_vec(1), r_vec(2)) / norm(cross(r_vec(1), r_vec(2)));
C_23 = cross(r_vec(2), r_vec(3)) / norm(cross(r_vec(2), r_vec(3)));
% C_31 = cross(r_vec(3), r_vec(1)) / norm(cross(r_vec(3), r_vec(1)));

% check how coplanar the position vectors are
coplanar_val = dot((r_vec(1) / r_mag(1)), C_23);

% calculate N, D, and S
N = r_mag(1) * cross(r_vec(2), r_vec(3)) + r_mag(2) * cross(r_vec(3), r_vec(1)) + r_mag(3) * cross(r_vec(1), r_vec(2));
D = cross(r_vec(1), r_vec(2)) + cross(r_vec(2), r_vec(3)) + cross(r_vec(3), r_vec(1));
S = r_vec(1) * (r_mag(2) - r_mag(3)) + r_vec(2) * (r_mag(3) - r_mag(1)) + r_vec(3) * (r_mag(1) - r_mag(2));

% calculate velocity vector at second position vector
v_vec = sqrt(mu / (norm(N) * norm(D))) * ((cross(D, r_vec(2)) / r_mag(2)) + S);

end