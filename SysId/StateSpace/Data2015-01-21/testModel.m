clc;
clear;
close all;

%load step data
extractStepsFromHelmsman;

%tool
addpath('../tools/');

%% Id

%select which seq yu want to use for the identification
nameSeqId = 'tack5B';

eval(['seqId = stepTacks.' num2str(nameSeqId) ';']);

%least square solution
[a, b, Dt] = tool_computeBestScalarValues(seqId);

%% validation

%select which seq yu want to use for the validation
nameSeqVal = 'tack5B';
eval(['seqVal = stepTacks.' num2str(nameSeqVal) ';']);

rudderVal = seqVal.dataYaw.InputData;
yawVal = seqVal.dataYaw.OutputData;
yawRateVal = seqVal.dataYawRate.OutputData;

numbSamples = length(rudderVal);

wModel = zeros(numbSamples, 1);
yModel = zeros(numbSamples, 1);

%start with same initial condition
wModel(1) = yawRateVal(1);
yModel(1) = yawVal(1);
yModel(2) = yawVal(2);

%compute model output
for k = 1 : numbSamples - 1
    wModel(k+1) = a * wModel(k) + b * rudderVal(k);
    
    if(k >= 2)
        yModel(k+1) = yModel(k) + Dt * wModel(k);
    end
end

%compare validation data and model output
lW = 1.3;
figure;

subplot(2,1,1);
plot(yawRateVal .* 180 / pi, 'b', 'LineWidth', lW);
hold on;
plot(wModel .* 180 / pi, 'm--', 'LineWidth', lW);
grid on;
title('yawRate');

subplot(2,1,2);
plot(yawVal .* 180 / pi, 'b', 'LineWidth', lW);
hold on;
plot(yModel .* 180 / pi, 'm--', 'LineWidth', lW);
grid on;
title('yaw');



