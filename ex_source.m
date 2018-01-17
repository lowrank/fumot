function [ ret ] = ex_source( x )
%EX_SOURCE excitation source.
%   x:   coordinate
%   ret: must be a C^2 function.
    ret =  exp(2* x(1,:)) + exp(-2 * x(2,:));
%     ret = 1 + 0 * x(1,:);
end

