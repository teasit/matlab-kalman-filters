classdef (Abstract) SigmaPointsSelector < handle
    %SIGMAPOINTSSELECTOR Abstract interface class.
    methods (Abstract)
        %SELECT Must return Sigma Points and covariance weights.
        [sigmaPnts, wCov, wMean] = select
    end
end

