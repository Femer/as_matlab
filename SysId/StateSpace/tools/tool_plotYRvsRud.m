function tool_plotYRvsRud(tStart, tStop, dataStruct, nameFig)

    
    

    %select index
    att_index = tool_selectIndexOnTime(tStart, tStop, dataStruct.att.time_att);
    rc_index = tool_selectIndexOnTime(tStart, tStop, dataStruct.rc.time_rc);
    
    %take data
    att_time = dataStruct.att.time_att(att_index) ./ 1e6;
    yawRate = dataStruct.att.ATT_YawRate(att_index);
    
    rc_time = dataStruct.rc.time_rc(rc_index) ./ 1e6;
    rudder = dataStruct.rc.RC_Ch0(rc_index);
    mode = dataStruct.rc.RC_Ch4(rc_index);
    
    [u, v, ~, est0_time] = tool_extractVelocites(tStart, tStop, dataStruct.est0);
    est0_time = est0_time ./ 1e6;
    %plot yawRate vs rudder command
    lineW = 1.2;
    
    figure;
    set(gcf,'name', nameFig,'numbertitle','off');
    
    subplot(2,1,1);
    plot(att_time, yawRate, 'b--', 'LineWidth', lineW);
    hold on;
    plot(rc_time, rudder, 'm--', 'LineWidth', lineW);
    plot(rc_time, mode, 'r');
    legend('yawRate', 'rud', 'mode');
    grid on;
    
    subplot(2,1,2);
    plot(est0_time, u, 'c');
    hold on;
    plot(est0_time, v, 'm');
    legend('u', 'v');
    grid on;
    
end