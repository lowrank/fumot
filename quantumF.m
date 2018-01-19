function [ ret ] = quantumF( x )
%QUANTUMF Quantum Efficiency coef for fluorescence.
%   x coordinate, vectorized.
%   ret ~0.7
    ret = 0.7 + 0.1 * sin(6*pi * x(1,:)) .* sin(6 *pi* x(2,:));
end

