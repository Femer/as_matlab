function tool_validateModel(handles, model, validationLog, resampleFactor)

%create iddata from validationLog
meanTs = mean(diff(validationLog.time)); %in milliseconds

dataId.dataYawRate = iddata(validationLog.yawRate, ...
                            validationLog.rudder,...
                            meanTs);
                        
dataId.dataYaw = iddata(validationLog.yaw, ...
                        validationLog.rudder,...
                        meanTs);
                    
dataId.time_sysId = validationLog.time;

%convert resampleFactor from seconds to validation data time
if(resampleFactor ~= -1)
    resampleFactor = round(resampleFactor / (meanTs / 1e3));%remember: meanTs is in mSec 
end
%compute model response
[yModel, wModel] = tool_computeModelResponse(model, dataId, resampleFactor);

%plot real response Vs Model response

%clear axes
cla(handles.a_yawRate);
cla(handles.a_yaw);
cla(handles.a_rudder);
%time in sec
time_sec = validationLog.time ./ 1e3;

lW = 1.7;

%plot yawRate in deg/s
plot(handles.a_yawRate, time_sec, validationLog.yawRate .* 180 / pi, ...
    'b--*', 'LineWidth', lW);
hold(handles.a_yawRate, 'on');
plot(handles.a_yawRate, time_sec, wModel .* 180 / pi, ...
    'c--*', 'LineWidth', lW);
legend(handles.a_yawRate, 'real', 'model');

grid(handles.a_yawRate);
xlabel(handles.a_yawRate, 'Time [s]');
ylabel(handles.a_yawRate, 'w [deg/s]');

%plot yaw in deg
plot(handles.a_yaw, time_sec, validationLog.yaw .* 180 / pi, ...
    'b--*', 'LineWidth', lW);
hold(handles.a_yaw, 'on');
plot(handles.a_yaw, time_sec, yModel .* 180 / pi, ...
    'c--*', 'LineWidth', lW);
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

