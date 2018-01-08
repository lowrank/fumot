function [ ret ] = diffusionFM( x )
%DIFFUSIONFM Diffusion coef for emission stage
%   x is coordinate, vectorized.
%   ret ~ 0.1
    ret = 0.1 + 0.02 * cos(2 * x(1,:)) .* cos(2 * x(2,:));
end

