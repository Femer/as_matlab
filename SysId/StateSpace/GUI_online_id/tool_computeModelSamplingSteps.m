function samplingStep = tool_computeModelSamplingSteps(model, dataValidation)

%compute the sampling time of the mode, that is, every how many steps
%the input of the validation data must be fed into the model.
%the minimum value is necessary 1

%pay attention: model.Dt is in seconds, meanTs in milliseconds, so convert
%meanTs in seconds!
samplingStep = floor(model.Dt / (dataValidation.meanTs / 1e3));

%safety check
if(samplingStep < 1)
    samplingStep = 1;
end

end

