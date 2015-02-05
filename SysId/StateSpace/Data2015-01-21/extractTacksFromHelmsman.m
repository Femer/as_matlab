%create data from tack helmsmen

clc;
clear;

showPlot = 0;

% %test unifor mean sample time for every sequence, in uSec
% meanTs = 9855.6988;

data.att = load('test1Att');
data.rc = load('test1RC');
data.bgud = load('test1BGUD');
data.gps = load('test1GPS');
data.qgc2 = load('test1QGC2');
data.wsai = load('test1WSAI');


data.est = load('test1EST0');
data.gpos = load('test1GPOS');

%add path for extra helper functions
addpath('../tools/');

%% identify tack1, TYPE 0, no velocity after tack
% start = ((71.6 * 60) + 12.9)* 1e6;
% stop = start + (5.3 * 1e6);
% 
% tack1 = tool_selectData(start, stop, data);
% 
% tack = tack1;
% dataYawRate = iddata(tack.yawRate_interp, tack.rudder);
% dataYawRate.OutputName = 'yawRate';
% dataYawRate.InputName = 'rudder';
% dataYaw = iddata(tack.yaw_interp, tack.rudder);
% dataYaw.OutputName = 'yaw';
% dataYaw.InputName = 'rudder';
% stepTacks.tack1.dataYawRate = dataYawRate;
% stepTacks.tack1.dataYaw = dataYaw;
% stepTacks.tack1.bgud_time = tack.bgud_time;
% 
% if(showPlot)
%     tool_plot(tack1, 'tack1');
%     figure;
%     plot(dataYawRate);
% end

%% identify tack3, TYPE 1. rud_45 = 0.85; alpha_hat = 45;
%qiute good, only some overshoot after tack
start = ((79.65 * 60) + 15.8)* 1e6;
stop = start + (2.8 * 1e6);

tack3 = tool_selectData(start, stop, data);

stepTacks.tack3 = tool_createIdData(tack3);

if(showPlot)
    tool_plot(tack3, 'tack3');
    figure;
    plot(dataYawRate);
end

%% identify tack4, TYPE 1, double one rud_45 = 0.85; alpha_hat = 45;
%double tack, not so bad but we needed some more velocity after the first
%tack
start = ((82.3 * 60) + 16.5)* 1e6;
stop = start + (3.2 * 1e6);

tack4 = tool_selectData(start, stop, data);

stepTacks.tack4 = tool_createIdData(tack4);

if(showPlot)
    tool_plot(tack4, 'tack4');
    figure;
    plot(dataYawRate);
end

%% identify tack5, TYPE 1, double one rud_45 = 0.8; alpha_hat = 40;
%not bad, maybe the best type 1 tack
start = ((87.8 * 60) + 20.49)* 1e6;
stop = start + (4.1 * 1e6);

tack5 = tool_selectData(start, stop, data);

stepTacks.tack5 = tool_createIdData(tack5);

if(showPlot)
    tool_plot(tack5, 'tack5');
    figure;
    plot(dataYawRate);
end

%% identify tack5B, TYPE 1, double one rud_45 = 0.8; alpha_hat = 40;
%not bad, maybe the best type 1 tack
start = ((87.8 * 60) + 35.9)* 1e6;
stop = start + (3.51 * 1e6);

tack5B = tool_selectData(start, stop, data);

stepTacks.tack5B = tool_createIdData(tack5B);

if(showPlot)
    tool_plot(tack5B, 'tack5B');
    figure;
    plot(dataYawRate);
end
%% identify tack8, TYPE 0 rud_45 = 0.8;
%not enough velocity after the tack
start = ((103.2 * 60) + 17.5)* 1e6;
stop = start +(6.66 * 1e6);

tack8 = tool_selectData(start, stop, data);

stepTacks.tack8 = tool_createIdData(tack8);

if(showPlot)
    tool_plot(tack8, 'tack8');
    figure;
    plot(dataYawRate);
end
%% identify tack6, ERROR IN TWD MEAN DIRECTION CHANGED ON THE GO
start = ((94 * 60) + 30)* 1e6;
stop = start +(6.26 * 1e6);

tack6 = tool_selectData(start, stop, data);

stepTacks.tack6 = tool_createIdData(tack6);

if(showPlot)
    tool_plot(tack6, 'tack6');
    figure;
    plot(dataYawRate);
end
%% identify tack7, TYPE 0 ERROR IN TWD MEAN DIRECTION CHANGED ON THE GO
% start = ((101 * 60) + 10)* 1e6;
% stop = start +(6.87 * 1e6);
% 
% tack7 = tool_selectData(start, stop, data);
% 
% tack = tack7;
% dataYawRate = iddata(tack.yawRate_interp, tack.rudder);
% dataYawRate.OutputName = 'yawRate';
% dataYawRate.InputName = 'rudder';
% dataYaw = iddata(tack.yaw_interp, tack.rudder);
% dataYaw.OutputName = 'yaw';
% dataYaw.InputName = 'rudder';
% stepTacks.tack7.dataYawRate = dataYawRate;
% stepTacks.tack7.dataYaw = dataYaw;
% stepTacks.tack7.bgud_time = tack.bgud_time;
% 
% if(showPlot)
%     tool_plot(tack7, 'tack7');
%     figure;
%     plot(dataYawRate);
% end
%%
clearvars data showPlot start stop dataYaw dataYawRate 
clearvars -regexp ^tack