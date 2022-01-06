function errorNorm = nrmse(xRef,xEst)
%NRMSE Calculates the Normalized Root Mean Square Error for predictions
%compared to a reference dataset.
error = estimation.rmse(xRef,xEst);
errorNorm = error/(max(xRef)-min(xRef));
end

