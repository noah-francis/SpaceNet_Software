% "a random normal distribution of coordinates gives you a uniform 
%   distribution of directions"
% https://stackoverflow.com/a/9751925/14382550

% generate same plot from post (uniform spread of directions of magnitude 1)
figure(1)
v = randn(1000, 3); % normal distribution of coordinates
v = bsxfun(@rdivide, v, sqrt(sum(v.^2, 2))); % uniform distribution of directions
plot3(v(:,1), v(:,2), v(:,3), '.')
axis equal

% want the same thing, but with randomly distributed vector magnitudes
figure(2)
v = randn(1000, 3); % normal distribution of coordinates
v = bsxfun(@rdivide, v, sqrt(sum(v.^2, 2))); % uniform distribution of directions

err_mag_UNIFORM = 100 * rand(1000, 1); 
err_mag_NORMAL = 100 * randn(1000, 1);

err_mag = err_mag_UNIFORM;
% err_mag = err_mag_NORMAL;

v = err_mag .* v; % scaling each direction by a randomly distributed length (err_mag)
plot3(v(:,1), v(:,2), v(:,3), '.')
axis equal

%% My thinking is that we want a uniform spread of directions but a normal spread of magnitudes