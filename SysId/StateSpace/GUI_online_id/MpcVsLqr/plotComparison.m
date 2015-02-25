function plotComparison(lqrData, mpcData, params, deltas, ...
                        mpcParams, predHor_steps, weights, ...
                        addNoise, mpcModelDt)


figure;

%usefull index
yawRateIndex = 1;
yawIndex = 2;
%another useful index for the extended state
lastRudderIndex = 3;

set(gcf,'name', ...
    ['MPC (forces) vs LQR. Prediction Horizon: ' num2str(predHor_steps) ...
    ' steps (' num2str(predHor_steps * mpcModelDt) ...
    ' sec). Weights:' num2str(weights)], ...
    'numbertitle', 'off');

%line width
lW0 = 1.3;

%legend, see if we added noise in the simulation
if(addNoise == 1)
    leg = { {'w real MPC', 'w by KF'}, ...
            {'w real LQR', 'w by KF'}, ...
            {'\psi real MPC', '\psi by KF'}, ...
            {'\psi real LQR', '\psi by KF'}, ...
            'rud MPC',...
            'rud LQR'};
else
    leg = { 'w real MPC', ...
            'w real LQR', ...
            '\psi real MPC', ...
            '\psi real LQR', ...
            'rud MPC',...
            'rud LQR'};
end

%real yaw and yawRate
yawReal = [mpcData.xHatSim(yawIndex, :);
          lqrData.xHatSim(yawIndex, :)];

yawRateReal = [mpcData.xHatSim(yawRateIndex, :);
               lqrData.xHatSim(yawRateIndex, :)];

%rudder for the real system, NOT for the extended one
rudder{1} = mpcData.rudCmd;
rudder{2} = lqrData.rudCmd;

rudderLimit(1) = mpcParams.upperBound(2);
rudderLimit(2) = params.rudderMax;

controllerTimeStr{1} = mpcData.time;
controllerTimeStr{2} = lqrData.time;

%state estimated by the time varying Kalman filter
%take data from index 2 since at index1 the tack was not started yet
yawEstVector = [mpcData.xHatEst(yawIndex,2:end);
                lqrData.xHatEst(yawIndex,2:end)];

%state estimated by the time vayring Kalman filter
yawRateEstVector = [mpcData.xHatEst(yawRateIndex,2:end);
                    lqrData.xHatEst(yawRateIndex,2:end)];
%full simulation time vector
time = (0: params.N - 1) * params.realModelDt;

for i = 1 : 2
    index = i;
    
    controllerTime = controllerTimeStr{i} * params.realModelDt;
    
    h1 = subplot(3, 2, index);
            
    plot(time, yawRateReal(i, :) .* 180 / pi, ...
        'LineWidth', 1.9, 'Color', [88 25 225] ./ 255);
    hold on;
    
    %see if we added noise in the simulation
    if(addNoise == 1)
        plot(controllerTime(2:end), yawRateEstVector(i, :) .* 180 / pi, 'c-.', ...
            'LineWidth', 1.9, 'Color', [245 86 1] ./ 255);
    end
    
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
    %see if we added noise in the simulation
    if(addNoise == 1)
        plot(controllerTime(2:end), yawEstVector(i, :) .* 180 / pi, 'c-.', ...
            'LineWidth', 1.9, 'Color', [245 86 1] ./ 255);
    end
    
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
         [rudderLimit(i) rudderLimit(i)], 'r-.', 'LineWidth', lW0);
    plot([controllerTime(1) controllerTime(end)], ...
         [-rudderLimit(i) -rudderLimit(i)], 'r-.', 'LineWidth', lW0);
    
    plot([controllerTime(1) controllerTime(end)], [deltas(3) deltas(3)], 'c--', 'LineWidth', 1.9);
    plot([controllerTime(1) controllerTime(end)], -[deltas(3)  deltas(3)], 'c--', 'LineWidth', 1.9);
    
    grid on;
    ylabel('Rudder [cmd]');
    xlabel('Time [sec]');
    xlim([controllerTime(1) controllerTime(end)]);
    legend(leg{index}, 'Location', 'southeast');
    
    linkaxes([h1, h2, h3], 'x');
    
    set(h3, 'XLim', [controllerTime(1), controllerTime(end)]);
end


end

