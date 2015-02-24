clc;
clear;
close all;

load('debugValues');

weights = [1,3,1,35];

sim_MpcVsLqr(realModel, lqrModel, mpcModel, predHor_steps, ...
                      weights, deltas, constraints, typeTack);
                  

data = mpcData;

% figure;
% subplot(321);
% plot(data.xHatSim(1,:));
% subplot(322);
% plot(data.xHatEst(1,2:end));
% 
% subplot(323);
% plot(data.xHatSim(2,:));
% subplot(324);
% plot(data.xHatEst(2,2:end));
% 
% subplot(325);
% %plot(data.rudCmd);
% plot(data.xHatSim(3,:));
% subplot(326);
% plot(data.xHatEst(3,:));