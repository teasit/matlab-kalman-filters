function y = quantize(y, resolution)
% QUANTIZE model the discrete resolution of a sensor
%
% Inputs:
%   y               signal to quantize (will lose precision)
%   resolution      the (sensor's) resolution causing quantization error.

y = resolution*round(y/resolution);
end

