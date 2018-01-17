function [ ret ] = absorptionFF( x )
%ABSORPTIONFF Absorption coef for fluorescence material.
%   x coordinate, vectorized.
%   ret ~0.2
N = 500;
[X,Y] =meshgrid(linspace(-0.5,0.5, N));
X = X(:);
Y = Y(:);
ph = phantom('Modified Shepp-Logan', N);
F = scatteredInterpolant(X,Y,ph(:));
ret = 0.1 * F(x(1,:), x(2,:)) + 0.05;
% ret = 0.2 + 0.05 * x(1,:);


end

