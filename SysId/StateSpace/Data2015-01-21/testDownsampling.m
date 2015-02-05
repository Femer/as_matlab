clc;
clear;
close all;


extractStepsFromHelmsman;
%%
groupNSamples = 5;

dataFull = stepTacks.tack5B;

[dataDownSampled] = tool_downsampleData(dataFull, groupNSamples);

fullO = dataFull.dataYaw.OutputData;
downO = dataDownSampled.dataYaw.OutputData;

fullR = dataFull.dataYaw.InputData;
downR = dataDownSampled.dataYaw.InputData;

fullT = dataFull.bgud_time;
downT = dataDownSampled.bgud_time;

figure;
subplot(2,1,1);
plot(fullT, fullO, 'b+');
hold on;
plot(downT, downO, 'm--*');

subplot(2,1,2);
plot(fullT, fullR, 'b+');
hold on;
plot(downT, downR, 'm--*');