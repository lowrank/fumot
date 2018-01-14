function [ ret ] = LaplaceSqrtDiffusionFX( x )
%LAPLACEDIFFUSIONFX calculate \nabla\cdot \nabla Dx.
%   must modify according to Dx.
%     ret = (0 * x(1,:) - 0.25 * (0.05)^2) ./(0.1 + 0.05 * x(1,:)).^(3/2);
ret = 0 * x(1,:);
end

