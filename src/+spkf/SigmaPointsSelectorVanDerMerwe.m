classdef SigmaPointsSelectorVanDerMerwe < spkf.SigmaPointsSelector
    %SIGMAPOINTSSELECTORVANDERMERWE Selection scheme according to R. Van der Merwe.
    %
    % Input arguments:
    %   (optional) alpha
    %   (optional) beta
    %   (optional) kappa
    
    properties(Access=public)
        %ALPHA Van der Merwe parameter for Sigma Point selection.
        %Typical range for parameter: 1E-3 < alpha <=1
        Alpha = 1E-3;
        
        %BETA Van der Merwe parameter for Sigma Point selection
        %Typical range for parameter: beta = 2 (if Gaussian priors)
        Beta = 2;
        
        %KAPPA Van der Merwe parameter for Sigma Point selection
        %Typical range for parameter: kappa = 0 (usually)
        Kappa = 0;
    end
    
    methods
        function obj = SigmaPointsSelectorVanDerMerwe(varargin)
            %SIGMAPOINTSSELECTORVANDERMERWE Constructs an instance.
            
            if nargin == 1
                obj.Alpha = varargin{1};
            elseif nargin == 2
                obj.Alpha = varargin{1};
                obj.Beta = varargin{2};
            elseif nargin == 3
                obj.Alpha = varargin{1};
                obj.Beta = varargin{2};
                obj.Kappa = varargin{3};
            elseif nargin > 3
                error("Too many input arguments")
            end
        end

        function [sigmaPnts, wCov, wMean] = select(obj, augmState, augmCovSqrt)
            %SELECT Selection of Sigma Points and weights.
            %
            %   Formula:
            %       sigma = [xa  xa+scale*sqrtPa  xa-scale*sqrtPa]
            %   Where:
            %       xa:     augmented state estimate
            %       scale:  scaling factor (depending on sigma point
            %               selection method)
            %       sqrtPa: square-root of augmented covariance matrix
            %   
            %   The resulting sigma points are given as a matrix of dimensions (n,2n+1)
            %   where n is the number of augmented states, not the (system) states.
            
            typename = class(augmState);
            
            L = numel(augmState);
            l = obj.Alpha^2*(L+obj.Kappa)-L;

            wMean = zeros(1,2, typename);
            wMean(1) = l/(L+l);
            wMean(2) = 1/(2*(L+l));

            wCov = wMean;
            wCov(1) = wCov(1) + (1-obj.Alpha^2+obj.Beta);
            
            augmCovSqrtScaled = sqrt(L+l)*augmCovSqrt;

            sigmaPnts = zeros(L,2*L+1, typename);
            sigmaPnts(:,2:L+1)   = +augmCovSqrtScaled;
            sigmaPnts(:,L+2:end) = -augmCovSqrtScaled;
            for i = 1:2*L+1
                sigmaPnts(:,i) = sigmaPnts(:,i) + augmState;
            end
        end
    end
end

