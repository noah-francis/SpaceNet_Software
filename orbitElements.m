function [i, Omega, e, omega, M, n, theta, h] = orbitElements(state_vec, mu)
    %   Use a given state vector (x, y, z, xdot, ydot, zdot)' to
    %   compute Keplerian orbital elements.
    %   -----------------------------------------------------------------------
    %
    %   Inputs:
    %   --------  
    %   `mu` -- Standard Gravitational Parameter of the body orbited [km^3 / s^2]
    %   `state_vec` -- vector with position and velocity components in the ECEF frame [km] and [km / s]
    %
    %   Returns:
    %   --------
    %   `i`     -- inclination [deg]
    %   `Omega` -- right ascension of the ascending node [deg]
    %   `e`     -- eccentricity [~]
    %   `omega` -- argument of perigee [deg]
    %   `M`     -- mean anomaly [deg]
    %   `n`     -- mean motion [rev / day]
    %   `theta` -- true anomaly [deg]
    %   `h`     -- magnitude of the specific angular momentum [km^2 / s]

    % distance
    r = norm(state_vec(1:3));

    % speed
    v = norm(state_vec(4:6));

    % radial speed 
    vr = dot(state_vec(1:3), state_vec(4:6)) / r;
    % if (vr > 0): satellite is flying away from perigee
    % if (vr < 0): satellite is flying toward perigee

    % specific angular momentum
    h_vec = cross(state_vec(1:3), state_vec(4:6));
    % magnitude of specific angular momentum
    h = norm(h_vec);

    % inclination
    i = acosd(h_vec(3) / h);

    % node line vector
    N_vec = cross([0,0,1], h_vec);

    % node line magnitude
    N = norm(N_vec);

    % right ascension of the ascending node
    if (N_vec(2) >= 0)
        Omega = acosd(N_vec(1) / N);
    else
        Omega = 360 - acosd(N_vec(1) / N);
    end

    % eccentricity vector
    e_vec = ((v^2 - (mu / r)) * state_vec(1:3) - r * vr * state_vec(4:6)) / mu;
    % eccentricity
    e = norm(e_vec);

    % argument of perigee
    if (e_vec(3) >= 0)
        omega = acosd(dot(N_vec, e_vec) / (N * e));
    else
        omega = 360 - acosd(dot(N_vec, e_vec) / (N * e));
    end

    % true anomaly
    if (vr >= 0)
        theta = acosd(dot(e_vec, state_vec(1:3)) / (e * r));
    else
        theta = 360 - acosd(dot(e_vec, state_vec(1:3)) / (e * r));
    end

    % eccentric anomaly
    E = 2 * atan(sqrt((1 - e) / (1 + e)) * tand(theta / 2));

    % mean anomaly
    M = E - e * sin(E);%[rad]
    M = rad2deg(M);%[deg]

    % orbital period
    T = 2 * pi * h^3 / (mu^2 * (1 - e^2)^(3/2));%[s]

    % mean motion 
    n = 2 * pi / T;%[rad / s]
    n = n * 3600 * 24 / (2 * pi);%[rev / day]
            
            
end