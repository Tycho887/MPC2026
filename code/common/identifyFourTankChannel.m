function G = identifyFourTankChannel(t, normalizedOutput, stepTime, order, useToolbox)
%IDENTIFYFOURTANKCHANNEL Fit a continuous first/second-order zero-free channel.
if useToolbox
    data = iddata(normalizedOutput,double(t>=stepTime),t(2)-t(1));
    options = tfestOptions('Display','off');
    fitted = tfest(data,order,0,options);
    G = tf(fitted.Numerator,fitted.Denominator);
else
    mask = t>=stepTime;
    G = fitContinuousStepModel(t(mask)-stepTime,normalizedOutput(mask),order);
end
end
