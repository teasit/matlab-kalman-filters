function f = eulerDiscretization(firstOrderFcn, dt)
%EULERDISCRETIZATION Discretizes first order differential equation function
%that has the form: dxdt = f(x) using Euler integration. The output then is
%a state transition function of the form: x = f(x) = x + dxdt*dt. In case
%the differential equation takes additional inputs (e.g. control inputs),
%these can also be passed through the discretization function.

% persistent DELTA_T; if isempty(DELTA_T); DELTA_T = dt; end
f = @(x, varargin) x + firstOrderFcn(x, varargin{:})*dt;
end

