function a = addCentripetalAcc(a, v, wz)
%ADDCENTRIPETALACC Add centripetal acceleration to acceleration in vehicle-
%fixed coordinates.
%
% Inputs:
%   a: cog-acceleration WITHOUT centripetal parts
%   v: [2 1] vector with x and y velocities
%   wz: angular velocity around z-axis

a(1) = a(1) - wz*v(2);
a(2) = a(2) + wz*v(1);
end

