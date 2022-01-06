function error = rmse(xRef,xEst)
%RMSE Calculates the Root Mean Square Error for predictions compared to a
%reference dataset.
xRefDwnSmpld = xRef(1:length(xRef)/length(xEst):end);
error = sqrt(mean((xRefDwnSmpld - xEst).^2));
end

