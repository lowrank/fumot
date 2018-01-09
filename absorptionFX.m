function [ ret ] = absorptionFX( x )
%ABSORPTIONFX Absorption coef for excitation stage
%   x is coordinate, vectorized.
%   ret ~0.1
    ret = 0.1 + 0.05 * sin(4 * x(1,:)) .* sin(4 * x(2,:));
%     ret = 0.1 + 0 * x(1,:);
end

