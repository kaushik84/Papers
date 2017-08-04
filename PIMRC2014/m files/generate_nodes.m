function [x,y]=generate_nodes(R,K)
%BINOMIAL Generation of a Binomial Point Process with K Nodes on a Disk of Radius R.

%
% Random placement of transmitters:

% Inverse Transform Sampling
r = R * sqrt(rand(1,K));
phi = 2 * pi * rand(1,K);

% Coordinate Transformation
x = r.*cos(phi);
y = r.*sin(phi);

end