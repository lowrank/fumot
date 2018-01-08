function [ ret ] = absorptionFM( x )
%ABSORPTIONFM Absorption coef for emission stage.
%  x coordinate, vectorized.
%  ret ~0.1
    ret = 0.1 + 0.02 * cos(4 * (x(1,:).^2 + x(2,:).^2));
end

