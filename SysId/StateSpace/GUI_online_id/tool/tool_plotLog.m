function tool_plotLog(handles, logStr)

%time in sec
time_sec = logStr.time ./ 1e3;

lW = 1.7;

%clear axes
cla(handles.a_yawRate);
cla(handles.a_yaw);
cla(handles.a_rudder);
%plot yawRate in deg/s
plot(handles.a_yawRate, time_sec, logStr.yawRate .* 180 / pi, ...
    'b--*', 'LineWidth', lW);
grid(handles.a_yawRate);
xlabel(handles.a_yawRate, 'Time [s]');
ylabel(handles.a_yawRate, 'w [deg/s]');

%plot yaw in deg
plot(handles.a_yaw, time_sec, logStr.yaw .* 180 / pi, ...
    'b--*', 'LineWidth', lW);
grid(handles.a_yaw);
xlabel(handles.a_yaw, 'Time [s]');
ylabel(handles.a_yaw, '\psi [deg]');

%plot rudder
plot(handles.a_rudder, time_sec, logStr.rudder, ...
    'm--*', 'LineWidth', lW);
grid(handles.a_rudder);
xlabel(handles.a_rudder, 'Time [s]');
ylabel(handles.a_rudder, 'rud [cmd]');

%link axes when zoom or move
linkaxes([handles.a_yawRate, handles.a_yaw, handles.a_rudder], 'x');
xlim(handles.a_yawRate, [time_sec(1) time_sec(end)]);

end

