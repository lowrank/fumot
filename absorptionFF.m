function [ ret ] = absorptionFF( x )
%ABSORPTIONFF Absorption coef for fluorescence material.
%   x coordinate, vectorized.
%   ret ~0.2
ret = 0.2 + 0.15 * exp(-((x(1,:)+0.2).^2 + (x(2,:)+0.2).^2) * 36) +...
    0.25 * exp(-((x(1,:)-0.2).^2 + (x(2,:)-0.2).^2) * 36);
% ret = 0.2 + 0.05 * x(1,:);


end

