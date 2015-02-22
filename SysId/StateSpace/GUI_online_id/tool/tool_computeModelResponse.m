function [tModel, yModel, wModel] = tool_computeModelResponse(model, data, ...
                                                              timeSampleRealSys, ...
                                                              samplingStep)

rudderVal = data.dataYaw.InputData;
yawVal = data.dataYaw.OutputData;
yawRateVal = data.dataYawRate.OutputData;

A = model.A;
B = model.B;

numbSamples = length(rudderVal);


%x = zeros(2, numbSamples);

%start with same initial condition of the validation data
x(1, 1) = yawRateVal(1);
x(2, 1) = yawVal(1);
tModel(1) = data.time_sysId(1);
count = 1;

%compute model evolution
for k = 1 : numbSamples - 1
    %update the model state with the validation data every timeSampleRealSys steps
    if(mod(k, timeSampleRealSys) == 0 && timeSampleRealSys ~= -1)
        x(:, count) = [yawRateVal(k); yawVal(k)];
    end
    %compute model response every samplingStep steps
    if(mod(k, samplingStep) == 0)
        x(:, count+1) = A * x(:, count) + B * rudderVal(k);
        tModel(count+1) = data.time_sysId(k);
        count = count + 1;
    end
end

%extract yawRate and yaw
wModel = x(1, :);
yModel = x(2, :);

%return row vectors
yModel = yModel';
wModel = wModel';
end

