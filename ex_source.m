function [ ret ] = ex_source( x )
%EX_SOURCE excitation source.
%   x:   coordinate
%   ret: must be a C^2 function.
    ret = 1 + 0 * x(1,:)';
end

