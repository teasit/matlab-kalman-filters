function p_delta = tanhParameterEstimationModel(p_delta,p0,maxRelOff)
%TANHPARAMETERESTIMATIONMODEL Implements a parameter estimation model
%that is bounded by specifying an initial guess of the parameter and a
%relative offset to that initial guess. The model can be applied in "joint
%parameter estimation" as described by Rudolph Van der Merwe (Dissertation
%2004). The output is the difference added to the nominal (initial
%estimate) value.
%
% Inputs:
%   p               Current parameter estimate.
%   p0              Initial parameter guess (nominal value).
%   maxRelOff       Percentage that specifies how far from the nominal
%                   value the model can deviate.
% Outputs:
%   p_delta         Difference added to the nominal value.

scaleFac = maxRelOff.*abs(p0);
p_delta = atan(p_delta./scaleFac);
end

