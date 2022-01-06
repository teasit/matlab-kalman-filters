function [t,ax] = plotEstimatedVariables(tRef,xRef,tEst,varargin)
xEst = varargin;
numVars = numel(xEst);
if isvector(xEst{1})
    numTiles = 1;
else
    numTiles = size(xEst{1},2);
end
if numTiles == 1
    numTilesX = 1;
    numTilesY = 1;
elseif mod(numTiles,2)
    numTilesX = (numTiles+1)/2;
    numTilesY = (numTiles-1)/2;
else
    numTilesX = numTiles/2;
    numTilesY = numTiles/2;
end

t_ = tiledlayout(numTilesX,numTilesY);

colors = {'b','g'};
linestyle = {'-','-'};
linewidth = [1 1];

for i = 1:numTiles
    ax_(i) = nexttile;
    plot(tRef, xRef(:,i), 'Color', 'k', 'LineWidth', 1)
    hold on; grid on
    for j = 1:numVars
        xEstj = xEst{j};
%         plot(tEst, xEstj(:,i), 'Color', colors{j}, 'Marker', 'x',...
%             'MarkerSize', 2, 'LineStyle', 'None');
        plot(tEst, xEstj(:,i), 'Color', colors{j}, 'LineStyle', linestyle{j}, 'LineWidth', linewidth(j))
    end
end

linkaxes(ax_)

if nargout == 1
    t = t_;
elseif nargout == 2
    t = t_;
    ax = ax_;
end

end

