function tool_validateModel(handles, model, validationLog, stateUpdateTime)

%create iddata from validationLog
meanTs = mean(diff(validationLog.time)); %in milliseconds

dataId.dataYawRate = iddata(validationLog.yawRate, ...
                            validationLog.rudder,...
                            meanTs);
                        
dataId.dataYaw = iddata(validationLog.yaw, ...
                        validationLog.rudder,...
                        meanTs);
                    
dataId.time_sysId = validationLog.time;
dataId.meanTs = meanTs;

%convert stateUpdateTime from seconds to validation data time
if(stateUpdateTime ~= -1)
    stateUpdateTime = round(stateUpdateTime / (meanTs / 1e3));%remember: meanTs is in mSec 
end

%compute the sampling time of the mode, that is, every how many steps
%the input of the validation data must be fed into the model.
%the minimum value is necessary 1
samplingStep = tool_computeModelSamplingSteps(model, dataId);

%compute model response
[tModel, yModel, wModel] = tool_computeModelResponse(model, dataId, ...
                                                     stateUpdateTime, ...
                                                     samplingStep);
%plot real response Vs Model response

%clear axes
cla(handles.a_yawRate);
cla(handles.a_yaw);
cla(handles.a_rudder);

legend(handles.a_yawRate,'hide');
legend(handles.a_yaw,'hide');
legend(handles.a_rudder,'hide');
%time in sec
time_sec = validationLog.time ./ 1e3;
tModel = tModel ./ 1e3;

lW = 1.7;

%plot yawRate in deg/s
plot(handles.a_yawRate, time_sec, validationLog.yawRate .* 180 / pi, ...
    'b--*', 'LineWidth', lW);
hold(handles.a_yawRate, 'on');
plot(handles.a_yawRate, tModel, wModel .* 180 / pi, ...
    'Color', [245 86 1] ./ 255, 'LineStyle', '--', ...
    'Marker', '*', 'LineWidth', lW);
legend(handles.a_yawRate, 'real', 'model');

grid(handles.a_yawRate);
xlabel(handles.a_yawRate, 'Time [s]');
ylabel(handles.a_yawRate, 'w [deg/s]');

%plot yaw in deg
plot(handles.a_yaw, time_sec, validationLog.yaw .* 180 / pi, ...
    'b--*', 'LineWidth', lW);
hold(handles.a_yaw, 'on');
plot(handles.a_yaw, tModel, yModel .* 180 / pi, ...
    'Color', [245 86 1] ./ 255, 'LineStyle', '--', ...
    'Marker', '*', 'LineWidth', lW);
legend(handles.a_yaw, 'real', 'model');

grid(handles.a_yaw);
xlabel(handles.a_yaw, 'Time [s]');
ylabel(handles.a_yaw, '\psi [deg]');

%plot rudder
plot(handles.a_rudder, time_sec, validationLog.rudder, ...
    'm--*', 'LineWidth', lW);
grid(handles.a_rudder);
xlabel(handles.a_rudder, 'Time [s]');
ylabel(handles.a_rudder, 'rud [cmd]');

%link axes when zoom or move
linkaxes([handles.a_yawRate, handles.a_yaw, handles.a_rudder], 'x');
xlim(handles.a_yawRate, [time_sec(1) time_sec(end)]);
end

