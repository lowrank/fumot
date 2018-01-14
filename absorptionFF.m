function [ ret ] = absorptionFF( x )
%ABSORPTIONFF Absorption coef for fluorescence material.
%   x coordinate, vectorized.
%   ret ~0.2
ret = 0.2 + 0.15 * exp(-(x(1,:).^2 + x(2,:).^2) * 16);
% ret = 0.2 + 0.05 * x(1,:);


end

