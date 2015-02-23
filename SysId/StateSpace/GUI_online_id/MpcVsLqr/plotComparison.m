function plotComparison(lqrData, mpcData, params, deltas)

figure;

%usefull index
yawRateIndex = 1;
yawIndex = 2;
%another useful index for the extended state
lastRudderIndex = 3;

set(gcf,'name', ...
    ['MPC (forces), steps prediction Horizon: '], ...
    'numbertitle', 'off');

%line width
lW0 = 1.3;

%legend
leg = { {'w real MPC', 'w by KF'}, ...
        {'w real LQR', 'w by KF'}, ...
        {'\psi real MPC', '\psi by KF'}, ...
        {'\psi real LQR', '\psi by KF'}, ...
        'rud MPC',...
        'rud LQR'};



yawReal = [mpcData.xHatSim(yawIndex, :);
          lqrData.xHatSim(yawIndex, :)];

yawRateReal = [mpcData.xHatSim(yawRateIndex, :);
               lqrData.xHatSim(yawRateIndex, :)];

%rudder for the real system, NOT for the extended one
rudder{1} = mpcData.rudCmd;
rudder{2} = lqrData.rudCmd;

controllerTimeStr{1} = mpcData.time;
controllerTimeStr{2} = lqrData.time;

%state estimated by the time varying Kalman filter
yawEstVector = [mpcData.xHatEst(yawIndex,:);
                lqrData.xHatEst(yawIndex,:)];

%state estimated by the time vayring Kalman filter
yawRateEstVector = [mpcData.xHatEst(yawRateIndex,:);
                    lqrData.xHatEst(yawRateIndex,:)];
%full simulation time vector
time = (0: params.N - 1) * params.realModelDt;

for i = 1 : 2
    index = i;
    
    controllerTime = controllerTimeStr{i} * params.realModelDt;
    
    h1 = subplot(3, 2, index);
        
    plot(time, yawRateReal(i, :) .* 180 / pi, ...
        'LineWidth', 1.9, 'Color', [88 25 225] ./ 255);
    hold on;
    
    plot(controllerTime, yawRateEstVector(i, :) .* 180 / pi, 'c-.', ...
        'LineWidth', 1.9, 'Color', [245 86 1] ./ 255);
    
    plot([controllerTime(1) controllerTime(end)], [deltas(1) deltas(1)], 'c--', 'LineWidth', 1.9);
    plot([controllerTime(1) controllerTime(end)], -[deltas(1)  deltas(1)], 'c--', 'LineWidth', 1.9);
    grid on;
    ylabel('[deg/s]');
    xlabel('Time [sec]');
    xlim([time(1) time(end)]);
    legend(leg{index});
    legend(leg{index}, 'Location', 'southeast');
    
    
    index = index + 2;
    
    h2 = subplot(3, 2, index);
    plot(time, yawReal(i, :) .* 180 / pi, ...
        'LineWidth', 1.9, 'Color', [88 25 225] ./ 255);
    hold on;
    
    plot(controllerTime, yawEstVector(i, :) .* 180 / pi, 'c-.', ...
        'LineWidth', 1.9, 'Color', [245 86 1] ./ 255);
    
    plot([controllerTime(1) controllerTime(end)], [deltas(2) deltas(2)], 'c--', 'LineWidth', 1.9);
    plot([controllerTime(1) controllerTime(end)], -[deltas(2)  deltas(2)], 'c--', 'LineWidth', 1.9);
    
    grid on;
    ylabel('[deg]');
    xlabel('Time [sec]');
    xlim([time(1) time(end)]);
    legend(leg{index});
    
    
    index = index + 2;
    
    h3 = subplot(3, 2, index);
    
    plot(controllerTime, rudder{i}, 'm--*', 'LineWidth', 1.4);
    hold on;
    plot([controllerTime(1) controllerTime(end)], ...
         [params.rudderMax params.rudderMax], 'r-.', 'LineWidth', lW0);
    plot([controllerTime(1) controllerTime(end)], ...
         [-params.rudderMax -params.rudderMax], 'r-.', 'LineWidth', lW0);
    
    plot([controllerTime(1) controllerTime(end)], [deltas(3) deltas(3)], 'c--', 'LineWidth', 1.9);
    plot([controllerTime(1) controllerTime(end)], -[deltas(3)  deltas(3)], 'c--', 'LineWidth', 1.9);
    
    grid on;
    ylabel('Rudder [cmd]');
    xlabel('Time [sec]');
    xlim([controllerTime(1) controllerTime(end)]);
    legend(leg{index}, 'Location', 'southeast');
    
    linkaxes([h1, h2, h3], 'x');
    
    %set(h3, 'XLim', [controllerTime(1), controllerTime(end)]);
end


end

