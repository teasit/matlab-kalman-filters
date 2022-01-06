function [t,ax] = plotEstimatedVariableWithCovariance(tRef,xRef,tEst,xEst,cEst,varargin)
t_ = tiledlayout(2,1);
ax_(1) = nexttile;
plot(tRef, xRef, 'Color', 'k', 'LineWidth', 1)
hold on; grid on
% plot(tEst, xEst, 'Color', 'b', 'Marker', 'x',...
%     'MarkerSize', 2, 'LineStyle', 'None')
for i = 1:size(xEst,3)
    plot(tEst, xEst(:,:,i), 'Color', 'b', 'LineStyle', '-', 'LineWidth', 1)
end
if nargin > 5
    xMeasured = varargin{:};
    plot(tEst, xMeasured, 'Color', 'g', 'LineStyle', '-', 'LineWidth', 1)
    legend(ax_(1), "Referenz", "SRUKF", "Messwert")
else
    legend(ax_(1), "Referenz", "SRUKF")
end
ax_(2) = nexttile;
hold on; grid on
cEst = squeeze(cEst);
for i = 1:size(cEst,2)
    plot(tEst, cEst(:,i), 'Color', 'b')
end
% legend(ax_(2), "Covariance")

if nargout == 1
    t = t_;
elseif nargout == 2
    t = t_;
    ax = ax_;
end

end

