function [ ret ] = diffusionFX( x )
%DIFFUSIONFX Diffusion coef for excitation stage
%   x is coordinate, vectorized. 
%   ret ~0.1
%     ret = 0.1 + 0.05 * (x(1,:).^2 + x(2,:).^2);
    ret = 0.1 + 0.00 * x(1,:);
end

