function [yModel, wModel] = tool_computeModelResponse(model, data, timeSampleRealSys)

rudderVal = data.dataYaw.InputData;
yawVal = data.dataYaw.OutputData;
yawRateVal = data.dataYawRate.OutputData;

A = model.A;
B = model.B;

numbSamples = length(rudderVal);


x = zeros(2, numbSamples);

%start with same initial condition
x(1, 1) = yawRateVal(1);
x(2, 1) = yawVal(1);

%compute model evolution
for k = 1 : numbSamples - 1
    %read the ral state very timeSampleRealSys steps
    if(mod(k, timeSampleRealSys) == 0)
        x(:, k) = [yawRateVal(k); yawVal(k)];
    end
    x(:, k+1) = A * x(:, k) + B * rudderVal(k);
end

%extract yawRate and yaw
wModel = x(1, :);
yModel = x(2, :);

% %limit yModel to be in [-pi, pi]
% indexOverPi = yModel > pi;
% indexLessMinusPi = yModel < - pi;
% 
% yModel(indexOverPi) = yModel(indexOverPi) - 2 * pi;
% yModel(indexLessMinusPi) = yModel(indexLessMinusPi) + 2 * pi;

%return row vectors
yModel = yModel';
wModel = wModel';
end

