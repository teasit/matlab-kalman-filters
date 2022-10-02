classdef HinfSteadyStateFilter < handle
    %HINFSTEADYSTATEFILER Steady-State H-Infinity Filter.
    
    properties (SetAccess = protected)
        Gain
        State
        StateCovariance {mustBeSymmetric(StateCovariance, 'P')}
        ProcessCovariance {mustBeSymmetric(ProcessCovariance, 'Q')}
        MeasurementCovariance {mustBeSymmetric(MeasurementCovariance, 'R')}
        EstimationOutput {mustBeSymmetric(EstimationOutput, 'L')}
        EstimationCovariance {mustBeSymmetric(EstimationCovariance, 'S')}
        PerformanceBound (1,1) {mustBeNonnegative}
        System ss {mustBeDiscreteStateSpace} = ss([],[],[],[],1)
    end
    
    methods (Access = protected)
        function determineGainMatrix(filter)
            sys = filter.System;
            A = sys.A;
            B = sys.B;
            C = sys.C;
            Q = filter.ProcessCovariance;
            R = filter.MeasurementCovariance;
            L = filter.EstimationOutput;
            S = filter.EstimationCovariance;
            S_ = L'*S*L;
            theta = filter.PerformanceBound;
            I = eye(size(B,1));
            
            P = idare(A', I, Q, inv(C'/R*C-theta*S_));
            filter.StateCovariance = P;
        end
    end
    methods
        function filter = HinfSteadyStateFilter(sys,Q,R,x0,theta,L,S)
            arguments
                sys ss {mustBeDiscreteStateSpace}
                Q  = []
                R  = []
                x0 = []
                theta = []
                L  = []
                S  = []
            end
            nx = size(sys.A,1);
            ny = size(sys.C,1);
            
            if isempty(Q)
                Q = eye(nx);
            end
            if isempty(R)
                R = eye(ny);
            end
            if isempty(x0)
                x0 = zeros(nx,1);
            end
            if isempty(theta)
                theta = eps;
            end
            if isempty(L)
                L = eye(nx);
            end
            if isempty(S)
                S = eye(nx);
            end
            
            filter.System = sys;
            filter.ProcessCovariance = Q;
            filter.MeasurementCovariance = R;
            filter.State = x0;
            filter.PerformanceBound = theta;
            filter.EstimationOutput = L;
            filter.EstimationCovariance = S;
            
            determineGainMatrix(filter)
        end
    end
end

function mustBeSymmetric(A, name)
[m,n] = size(A);
assert(m==n, 'Matrix ''%s'' must be symmetric.', name)
end

function mustBeDiscreteStateSpace(sys)
assert(sys.Ts, 'System must be a discrete state-space model.')
end
