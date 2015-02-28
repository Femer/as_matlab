function tool_plotFiltFiltLog(handles, filtfiltLog)

%take the selected log struct
[error, logStr, ~] = tool_getSelectedLog(handles);

if(error == 0)
    %clear axes
    cla(handles.a_yawRate);
    cla(handles.a_yaw);
    cla(handles.a_rudder);

    legend(handles.a_yawRate,'hide');
    legend(handles.a_yaw,'hide');
    legend(handles.a_rudder,'hide');

    %time in sec
    time_sec = logStr.time ./ 1e3;

    lW = 1.7;

    %plot yawRate in deg/s
    plot(handles.a_yawRate, time_sec, logStr.yawRate .* 180 / pi, ...
        'b--*', 'LineWidth', lW);
    hold(handles.a_yawRate, 'on');
    plot(handles.a_yawRate, time_sec, filtfiltLog.yawRate .* 180 / pi, ...
        'Color', [178 130 3] ./ 255, 'LineStyle', '--', ...
        'Marker', 'o', 'LineWidth', lW);
    legend(handles.a_yawRate, 'real', 'filtfilt');

    grid(handles.a_yawRate);
    xlabel(handles.a_yawRate, 'Time [s]');
    ylabel(handles.a_yawRate, 'w [deg/s]');

    %plot yaw in deg
    plot(handles.a_yaw, time_sec, logStr.yaw .* 180 / pi, ...
        'b--*', 'LineWidth', lW);
    hold(handles.a_yaw, 'on');
    plot(handles.a_yaw, time_sec, filtfiltLog.yaw .* 180 / pi, ...
        'Color', [178 130 3] ./ 255, 'LineStyle', '--', ...
        'Marker', 'o', 'LineWidth', lW);
    legend(handles.a_yaw, 'real', 'filtfilt');

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

end

