function p = tanhParameterEstimationEvaluation(p_delta,p0,maxRelOff)
scaleFac = maxRelOff.*abs(p0);
p = p0 + scaleFac.*tanh(p_delta);
end

