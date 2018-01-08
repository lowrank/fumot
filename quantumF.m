function [ ret ] = quantumF( x )
%QUANTUMF Quantum Efficiency coef for fluorescence.
%   x coordinate, vectorized.
%   ret ~0.7
    ret = 0.7 + 0.1 * cos(2 * x(1,:)) .* cos(2 * x(2,:));
end

