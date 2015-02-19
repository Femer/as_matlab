
clc;
clear;
%close all;

%2 models: the first one is used to simulate the boat, it is the "real"
%model. The first model is selected using typeOfModel variable.
%The second model is loaded from idModel2015-02-18 ad it is the model used
%by the KF, the MPC and LQR definition. It is called the uncertModel

%type of model you want to use to simulate the boat, identified from
%2015-01-21, assuming it is the "real" one
realModelType = 'capital';%little (a,b) or capital (A,B)

%choose the factor you want to use to downsample the model used by the MPC
factorSampleTime = 10;

%load real model, load uncertModel (for LQR) and load uncert downsampled
%model (for MPC)
[realModel, uncertModel, uncertMPCModel] = loadModels(realModelType, factorSampleTime);

%load tunable parameters
tunablePramas;

%build extended model
[realModelExt, uncertModelExt, uncertMPCModelExt] = extentendModels(realModel, ...
                                                                    uncertModel, ...
                                                                    uncertMPCModel);
        
%initial state and guess for simulation
loadInitParams;


%% Build MPC using model downsampled
mpcController = buildMPC(uncertMPCModelExt, predHor, Q, R, ...
                         tack, rudderMax, rudderVelocity);


%% Build LQR using normally sampled model
display('Building LQR');
[K_LQR, ~, ~] = dlqr(uncertModelExt.A, uncertModelExt.B, Q, R); 

%% Sim MPC and LQR 

% ---- MPC ----
[xHatSimMPC, xHatEstMPC, rudMPC, timeRunMPC] = simMPC(  mpcController, ...
                                                        realModelExt, uncertModelExt, ...
                                                        xHatSimMPC1, guessP1_1, guessX1Hat,...
                                                        rudderBeforeTack, factorSampleTime, ...
                                                        N, measNoise, useRealState,...
                                                        convarianceStr, rudderMax);
% ---- LQR ----
[xHatSimLQR, xHatEstLQR, rudLQR] = simLQR(  K_LQR, realModelExt, ...
                                            uncertModelExt, xHatSimLQR1,...
                                            guessP1_1, guessX1Hat,...
                                            rudderBeforeTack, N, ...
                                            measNoise, useRealState, ...
                                            rudderVelocity, convarianceStr, ...
                                            rudderMax);

%% plot
% ---- plot comparison ----

%compute simulation time, starting from time 0
time = (0:N) .* realModel.Dt; 

%time when rudder input by MPC was computed, converted in seconds
timeRunMPC = timeRunMPC .* realModel.Dt;

% %compute yaw overshoot limit
% if(strcmp(tack, 'p2s'))
%     limitYaw = -absDeltaYaw;
% else
%     limitYaw = absDeltaYaw;
% end

figure;

set(gcf,'name', ...
    ['MPC (yalmip), steps prediction Horizon: ' num2str(predHor) ...
    '; type of model: ' realModelType], ...
    'numbertitle', 'off');

lW0 = 1.3;

if useRealState
    leg = { 'w real MPC', ...
            'w real LQR', ...
            '\psi real MPC', ...
            '\psi real LQR', ...
            'rud MPC',...
            'rud LQR'};
else
    leg = { {'w real MPC', 'w by KF'}, ...
            {'w real LQR', 'w by KF'}, ...
            {'\psi real MPC', '\psi by KF'}, ...
            {'\psi real LQR', '\psi by KF'}, ...
            'rud MPC',...
            'rud LQR'};
end

yawReal = [xHatSimMPC(yawIndex, :);
           xHatSimLQR(yawIndex, :)];
       
yawRateReal = [xHatSimMPC(yawRateIndex, :);
               xHatSimLQR(yawRateIndex, :)];
           
%rudder for the real system, NOT for the extended one    
rudder{1} = rudMPC;
rudder{2} = rudLQR;
      
rudderTime{1} = timeRunMPC;
rudderTime{2} = time;
    
%state estimated by the time vayring Kalman filter
yawEstVector = [xHatEstMPC(yawIndex,:);
                xHatEstLQR(yawIndex,:)];
            
%state estimated by the time vayring Kalman filter
yawRateEstVector = [xHatEstMPC(yawRateIndex,:);
                    xHatEstLQR(yawRateIndex,:)];

for i = 1 : 2
   index = i;
   
   h1 = subplot(3, 2, index);
   
   plot(time(2:end), yawRateReal(i, :) .* 180 / pi, ...
       'LineWidth', 1.9, 'Color', [88 25 225] ./ 255);
   hold on;
   
   if useRealState == 0
       plot(time(2:end-1), yawRateEstVector(i, :) .* 180 / pi, 'c-.', ...
           'LineWidth', 1.9, 'Color', [245 86 1] ./ 255);
   end
   grid on;
   ylabel('[deg/s]');
   xlabel('Time [sec]');
   xlim([time(1) time(end)]);
   legend(leg{index});
   legend(leg{index}, 'Location', 'southeast');
   
   
   index = index + 2;
   
   h2 = subplot(3, 2, index);
   plot(time(2:end), yawReal(i, :) .* 180 / pi, ...
       'LineWidth', 1.9, 'Color', [88 25 225] ./ 255);
   hold on;
   %plot([time(2) time(end)], [limitYaw limitYaw] .* 180 / pi, 'r-.', 'LineWidth', lW0);
   
   if useRealState == 0
       plot(time(2:end-1), yawEstVector(i, :) .* 180 / pi, 'c-.', ...
           'LineWidth', 1.9, 'Color', [245 86 1] ./ 255);
   end
   grid on;
   ylabel('[deg]');
   xlabel('Time [sec]');
   xlim([time(1) time(end)]);
   legend(leg{index});
      
   
   index = index + 2;
   
   h3 = subplot(3, 2, index);
   
   timeR = rudderTime{i};
   plot(timeR, rudder{i}, 'm--*', 'LineWidth', 1.4);
   hold on;
   plot([timeR(1) timeR(end)], [rudderMax rudderMax], 'r-.', 'LineWidth', lW0);
   plot([timeR(1) timeR(end)], [-rudderMax -rudderMax], 'r-.', 'LineWidth', lW0);
   grid on;
   ylabel('Rudder [cmd]');
   xlabel('Time [sec]');
   xlim([timeR(1) timeR(end)]);
   legend(leg{index}, 'Location', 'southeast');
   
   linkaxes([h1, h2, h3], 'x');
   
   %h3.XLim = [timeR(1), timeR(end)];
end
