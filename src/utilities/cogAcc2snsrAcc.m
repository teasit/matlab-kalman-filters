function a = cogAcc2snsrAcc(r, a, wz)
%SENSORPLACEMENTKINEMATICS Calculates the acceleration a sensor would
%measure from a known acceleration in the center of gravity.
%
% Inputs:
%   r: [2 1] vector with x and y positions relative to cog
%   a: cog-acceleration
%   wz: angular velocity around z-axis
a(1) = a(1) - wz^2*r(1);
a(2) = a(2) - wz^2*r(2);
end

