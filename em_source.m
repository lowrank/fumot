function [ ret ] = em_source( x )
%EM_SOURCE emission source, only for auxillary function.
%   x : coordinate
%   ret: must be a C^2 fuction
    ret = 1 + 0 * x(1,:);
end

