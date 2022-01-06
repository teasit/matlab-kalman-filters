classdef SquareRootUnscentedKalmanFilter < handle
    %SQUAREROOTUNSCENTEDKALMANFILTER Square-Root Unscented Kalman Filter
    %Implements the computationally more efficient version of the Unscented
    %Kalman Filter, as is defined in [1].
    %
    %   Inputs:
    %       f               state transition function
    %       h               measurement function
    %       x0              initial state estimate
    %       Px0             (optional) initial state error covariance
    %       Rv              (optional) process noise covariance
    %       Rn              (optional) measurement noise covariance
    %       sigmaSelector   (optional) sigma point selector
    %
    %   [1]:
    %       Merwe R, Wan E (2003) Sigma-Point Kalman Filters for Probabil-
    %       istic Inference in Dynamic State-Space Models.
    %       Proceedings of the Workshop on Advances in Machine Learning.
    
    properties(Access=public)
        %State (x) state estimate vector
        State (:,1)
    end
    properties(Access=public)
        %StateTransitionFcn (f) state transition function handle
        StateTransitionFcn
    end
    properties(Access=public)
        %MeasurementFcn (h) measurement function handle
        MeasurementFcn
    end
    properties(Access=public)
        %StateCovariance (Px) state estimate error covariance matrix
        StateCovariance (:,:)
    end
    properties(Access=public)
        %MeasurementNoiseCovariance (Rn) measurement noise covariance
        %function handle
        MeasurementNoiseCovariance (:,:)
    end
    properties(Access=public)
        %ProcessNoiseCovariance (Rv) process noise covariance function
        %handle
        ProcessNoiseCovariance (:,:)
    end
    properties(SetAccess=public, GetAccess=public)
        %SigmaPointsSelector Sigma Points selector object
        SigmaPointsSelector
    end

    properties(Access=protected, Hidden)
        pState (:,1) {isnumeric}
        pProcessNoiseMean (:,1) {isnumeric}
        pMeasurementNoiseMean (:,1) {isnumeric}
        pLengthState (1,1) {isnumeric}
        pLengthMeasurementNoise (1,1) {isnumeric}
        pLengthProcessNoise (1,1) {isnumeric}
        pSizeState (1,2) {isnumeric}
        pSizeMeasurement (1,2) {isnumeric}
        pSizeProcessNoise (1,2) {isnumeric}
        pSizeMeasurementNoise (1,2) {isnumeric}
        pStateCovSqrt (:,:) {isnumeric}
        pStateMeasurementCrossCov (:,:) {isnumeric}
        pMeasurementCovSqrt (:,:) {isnumeric}
        pMeasurementNoiseCovSqrt (:,:) {isnumeric}
        pProcessNoiseCovSqrt (:,:) {isnumeric}
    end
    
    properties(Constant=true, Hidden)
        DEFAULT_STATE_ERROR_COVARIANCE = 1;
        DEFAULT_MEASUREMENT_NOISE_COVARIANCE = 1;
        DEFAULT_PROCESS_NOISE_COVARIANCE = 1;
    end
    
    methods
        function set.StateTransitionFcn(obj, value)
            validateattributes(value,...
                {'function_handle'}, {'nonempty'},...
                mfilename(), 'StateTransitionFcn');
            obj.StateTransitionFcn = value;
        end
        function value = get.StateCovariance(obj)
            value = obj.pStateCovSqrt*obj.pStateCovSqrt';
        end
    end
    
    methods(Access=protected)
        function S = choleskyFactor(obj,wCov,sigmaPnts,mean)
            %CHOLESKYFACTOR Calculates Cholesky factor 'S' using QR-Decomposition.
            %Cholesky Update is necessary in case of negative weights.
            sigmaPntsDelta = sigmaPnts;
            sigmaPntsDelta = bsxfun(@minus, sigmaPntsDelta, mean);
            [~,S] = qr(sqrt(wCov(2))*(sigmaPntsDelta(:,2:end))',0);
            if sign(wCov(1)) == 1
                v = '+';
            else
                v = '-';
            end
            u = sqrt(complex(wCov(1)))*(sigmaPntsDelta(:,1));
            S = real(obj.cholUpdateFactor(complex(S),u,v));
        end
    end
    
    methods(Access=protected, Static)
        function S = cholUpdateFactor(S,U,v)
            %CHOLUPDATEFACTOR Rank 1 update to Cholesky factorization.
            %Wraps MathWorks implementation with some additional fallback
            %in case of failed cholesky update.
            if ~ischar(v)
                error("Input v must be a char of either '+' or '-'.")
            end
            for i = 1:size(U,2)
                u = U(:,i);
                if v == '+' % update
                    S = cholupdate(S,u,'+');
                elseif v == '-' % downdate
                    [S,p] = cholupdate(S,u,'-');
                    if p ~= 0
                        P = S'*S - u*u';
                        [~,Singulars,rightSingulars] = svd(P,'econ');
                        R = rightSingulars*sqrt(Singulars);
                        S = R';
                        if ~istriu(S)
                            [~,S] = qr(S);
                        end
                    end
                end
            end
            S = S';
        end
    end
    
    methods(Access=public)
        function [x, Px] = predict(obj, varargin)
            %PREDICT From the current state estimate and the state tran-
            %sition function, predict the next state estimate.
            %
            %   Inputs:
            %       u       (Optional) Control input.
            %
            %   Outputs:
            %       x       (Optional) Predicted state estimate.
            %       Px      (Optional) Predicted state error covariance.
            
            if nargin == 1
                u = [];
            elseif nargin == 2
                u = varargin{1};
            elseif nargin > 2
                error("Too many inputs.")
            end
            
            x_ = obj.State;
            v_ = obj.pProcessNoiseMean;
            n_ = obj.pMeasurementNoiseMean;
            
            lenX = length(x_);
            lenV = length(v_);
            
            % (3.195)
            augmState = [x_;v_;n_];
            augmCovSqrt = blkdiag(...
                obj.pStateCovSqrt,...
                obj.pProcessNoiseCovSqrt,...
                obj.pMeasurementNoiseCovSqrt);
            [sigmaPnts, wCov, wMean]...
                = obj.SigmaPointsSelector.select(augmState, augmCovSqrt);
            sigmaX = sigmaPnts(1:lenX(1),:);
            sigmaV = sigmaPnts(lenX(1)+1:lenX(1)+lenV(1),:);
            
            % Propagate Sigma Points through state transition function
            % (3.196)
            f = obj.StateTransitionFcn;
            for i = 1:size(sigmaX,2)
                if isempty(u)
                    sigmaX(:,i) = f(sigmaX(:,i), sigmaV(:,i));
                else
                    sigmaX(:,i) = f(sigmaX(:,i), sigmaV(:,i), u);
                end
            end
            
            % Weighted average of sigma points as new State Estimate
            % (3.197)
            for i = 1:size(sigmaX,1)
                x_(i) = wMean(1).*sigmaX(i,1) ...
                    + sum(wMean(2)*sigmaX(i,2:end),2);
            end
            
            % (3.198), (3.199)
            Sx_ = obj.choleskyFactor(wCov,sigmaX,x_);
            
            obj.State = x_;
            obj.pStateCovSqrt = Sx_;
            
            if nargout == 1
                x = x_;
            elseif nargout == 2
                x = x_;
                Px = obj.StateCovariance;
            elseif nargout > 2
                error('Too many output arguments')
            end
        end
        function [x, Px, res] = correct(obj, varargin)
            %CORRECT Updates the state estimate with new measurement.
            %
            %   Inputs:
            %       y       Current measurement.
            %       u       (Optional) Control input.
            %
            %   Outputs:
            %       x       (Optional) Corrected state estimate.
            %       Px      (Optional) Corrected state error covariance.
            %       res     (Optional) Residual (or innovations).
            
            if nargin == 1
                error(['A measurement "y" is required to perform a ' ...
                    'correction.'])
            elseif nargin == 2
                y = varargin{1};
                u = [];
            elseif nargin == 3
                y = varargin{1};
                u = varargin{2};
            elseif nargin > 3
                error("Too many inputs.")
            end
            
            x_ = obj.State;
            v_ = obj.pProcessNoiseMean;
            n_ = obj.pMeasurementNoiseMean;
            Sx_ = obj.pStateCovSqrt;
            
            typename = class(x_);
            
            lenX = length(x_);
            lenV = length(v_);
            lenY = length(y);
            
            % (3.195)
            augmState = [x_;v_;n_];
            augmCovSqrt = blkdiag(...
                obj.pStateCovSqrt,...
                obj.pProcessNoiseCovSqrt,...
                obj.pMeasurementNoiseCovSqrt);
            [sigmaPnts, wCov, wMean] = obj.SigmaPointsSelector.select(...
                augmState, augmCovSqrt);
            sigmaX = sigmaPnts(1:lenX(1),:);
            sigmaN = sigmaPnts(lenX(1)+lenV(1)+1:end,:);
            
            % (3.200)
            h = obj.MeasurementFcn;
            sigmaY = zeros(lenY,size(sigmaX,2), typename);
            for i = 1:size(sigmaX,2)
                if isempty(u)
                    sigmaY(:,i) = h(sigmaX(:,i), sigmaN(:,i));
                else
                    sigmaY(:,i) = h(sigmaX(:,i), sigmaN(:,i), u);
                end
            end
            
            % (3.201)
            if coder.target('MATLAB')
                y_ = wMean(1)*sigmaY(:,1) ...
                    + sum(wMean(2)*sigmaY(:,2:end),2);
            else
                % https://de.mathworks.com/matlabcentral/answers/458934
                y_ = wMean(1)*sigmaY(:,1) + ...
                    sum(bsxfun(@times,wMean(2),sigmaY(:,2:end)),2);
            end
            
            deltaX = zeros(size(sigmaX), typename);
            for i = 1:size(deltaX,2)
                deltaX(:,i) = sigmaX(:,i) - x_;
            end
            deltaY = zeros(size(sigmaY), typename);
            for i = 1:size(deltaY,2)
                deltaY(:,i) = sigmaY(:,i) - y_;
            end
            
            % (3.202), (3.203)
            Sy_ = obj.choleskyFactor(wCov,sigmaY,y_);
            
            % (3.204)
            Pxy_ = deltaX;
            for i = 1:size(Pxy_,1)
                Pxy_(i,1) = wCov(1)*Pxy_(i,1);
                Pxy_(i,2:end) = wCov(2)*Pxy_(i,2:end);
            end
            Pxy_ = Pxy_*deltaY';
            
            % (3.205)
            K_ = (Pxy_/Sy_')/Sy_;
            
            % (3.206)
            residual = y-y_;
            x_ = x_ + K_*(residual);
            
            % (3.207)
            U_ = K_*Sy_;
            
            % (3.208)
            Sx_ = obj.cholUpdateFactor(Sx_',U_,'-');
            
            obj.State = x_;
            obj.pStateCovSqrt = Sx_;
            obj.pMeasurementCovSqrt = Sy_;
            obj.pStateMeasurementCrossCov = Pxy_;
            
            if nargout == 1
                x = obj.State;
            elseif nargout == 2
                x = obj.State;
                Px = obj.StateCovariance;
            elseif nargout == 3
                x = obj.State;
                Px = obj.StateCovariance;
                res = residual;
            elseif nargout > 3
                error('Too many output arguments')
            end
        end
        function reset(obj, x0, Px0)
            %RESET Resets instance to initial values.
            
            %Parse initial state
            if ~isempty(x0)
                if iscolumn(x0)
                    obj.State = x0;
                else
                    error("State has to be a column vector!")
                end
            else
                error("No initial state specified.")
            end
            
            sizeState = size(x0);
            
            %Parse state error covariance
            if isempty(Px0)
                obj.pStateCovSqrt = diag(ones(sizeState))...
                    *obj.DEFAULT_STATE_ERROR_COVARIANCE;
            elseif issymmetric(Px0) && length(Px0) == sizeState(1)
                obj.pStateCovSqrt = chol(Px0, 'lower');
            elseif isscalar(Px0)
                obj.pStateCovSqrt = diag(...
                    ones(1,sizeState(1))*sqrt(Px0));
            else
                error(strcat("State Covariance has to be symmetric ", ...
                    "[n,n] matrix where 'n' is the number of rows of ", ...
                    "the state."))
            end
        end
        function obj = SquareRootUnscentedKalmanFilter(f, h, x0, Px0, Rv, Rn)
            %SQUAREROOTUNSCENTEDKALMANFILTER Constructs SR-UKF instance.
            
            %State transition function
            typename = class(x0);
            if ~isa(f, 'function_handle')
                error(["The state transition function has to be a "
                    "function handle!"])
            else
                obj.StateTransitionFcn = f;
            end
            
            %Measurement function
            if ~isa(h, 'function_handle')
                error(["The measurement function has to be a function "
                    "handle!"])
            else
                obj.MeasurementFcn = h;
            end
            
            %Initial state
            if ~isempty(x0)
                if iscolumn(x0)
                    obj.State = x0;
                else
                    error("State has to be a column vector!")
                end
            else
                error("No initial state specified.")
            end
            sizeState = size(x0);
            obj.pLengthState = length(x0);
            obj.pSizeState = sizeState;
            
            %State error covariance
            if isempty(Px0)
                obj.pStateCovSqrt = diag(ones(sizeState))...
                    *obj.DEFAULT_STATE_ERROR_COVARIANCE;
            elseif issymmetric(Px0) && length(Px0) == sizeState(1)
                obj.pStateCovSqrt = chol(Px0, 'lower');
            elseif isscalar(Px0)
                obj.pStateCovSqrt = diag(...
                    ones(1,sizeState(1))*sqrt(Px0));
            else
                error(strcat("State Covariance has to be symmetric ", ...
                    "[n,n] matrix where 'n' is the number of rows of ", ...
                    "the state."))
            end
            
            %Process noise covariance
            if isempty(Rv)
                processNoiseCovariance = ...
                    obj.DEFAULT_PROCESS_NOISE_COVARIANCE;
                obj.ProcessNoiseCovariance = processNoiseCovariance;
            elseif issymmetric(Rv)
                processNoiseCovariance = Rv;
                obj.ProcessNoiseCovariance = processNoiseCovariance;
            else
                error("Process noise covariance is not symmetric.")
            end
            lengthProcessNoise = length(processNoiseCovariance);
            sizeProcessNoise = [lengthProcessNoise 1];
            processNoiseMean = zeros(sizeProcessNoise, typename);
            obj.pLengthProcessNoise = lengthProcessNoise;
            obj.pSizeProcessNoise = sizeProcessNoise;
            obj.pProcessNoiseMean = processNoiseMean;
            obj.pProcessNoiseCovSqrt = sqrt(processNoiseCovariance);
            
            %Measurement noise covariance
            if isempty(Rn)
                measurementNoiseCovariance = ...
                    obj.DEFAULT_MEASUREMENT_NOISE_COVARIANCE;
                obj.MeasurementNoiseCovariance = ...
                    measurementNoiseCovariance;
            elseif issymmetric(Rn)
                measurementNoiseCovariance = Rn;
                obj.MeasurementNoiseCovariance = ...
                    measurementNoiseCovariance;
            else
                error("Measurement noise covariance is not symmetric.")
            end
            lengthMeasurementNoise = length(measurementNoiseCovariance);
            sizeMeasurementNoise = [lengthMeasurementNoise 1];
            measurementNoiseMean = zeros(sizeMeasurementNoise, typename);
            obj.pLengthMeasurementNoise = lengthMeasurementNoise;
            obj.pSizeMeasurementNoise = sizeMeasurementNoise;
            obj.pMeasurementNoiseMean = measurementNoiseMean;
            obj.pMeasurementNoiseCovSqrt = ...
                sqrt(measurementNoiseCovariance);
            
            obj.pSizeMeasurement = size(...
                h(x0, zeros(sizeMeasurementNoise)));
            
            obj.SigmaPointsSelector = ...
                spkf.SigmaPointsSelectorVanDerMerwe();
        end
    end
end

