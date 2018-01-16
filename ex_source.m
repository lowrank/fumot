function [ ret ] = ex_source( x )
%EX_SOURCE excitation source.
%   x:   coordinate
%   ret: must be a C^2 function.
    ret =  exp(sqrt(6) * x(1,:)) + exp(-sqrt(3) * x(2,:));
end

