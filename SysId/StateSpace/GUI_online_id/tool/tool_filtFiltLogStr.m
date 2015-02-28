function [error, filtFiltLog, nameFiltLog] = tool_filtFiltLogStr(handles)

filtFiltLog = [];
nameFiltLog = [];

%take the selected log struct
[error, logStr, nameLog] = tool_getSelectedLog(handles);

if(error == 0)
    %compute the butterworth filter using data set by the user
    cutOffFreq = str2double(get(handles.e_cutOffFreq, 'String'));
    filtOrder = str2double(get(handles.e_orderFilt, 'String'));

    %compute normalized cut off frequency
    normalizedCutF = (cutOffFreq * 2) / tool_computeAvgSamplingF(logStr.time);
    
    %butterworth filt design
    [b, a] = butter(filtOrder, normalizedCutF, 'low');
    
    %filtfilt yawRate and yaw
    filtFiltLog.yawRate = filtfilt(b, a, logStr.yawRate);
    filtFiltLog.yaw = filtfilt(b, a, logStr.yaw);
    
    %copy time and rudder as they are
    filtFiltLog.time = logStr.time;
    filtFiltLog.rudder = logStr.rudder;
    
    %name of the filtered log, if cutOffFreq is not an integere, use the
    %letter 'p' instead of the .
    strCutFreq = num2str(cutOffFreq);
    indexPoint = regexp(strCutFreq, '\.');
    if(not(isempty(indexPoint)))
        strCutFreq(indexPoint) = 'p';
    end
    
    nameFiltLog = [nameLog '_filt_cf' strCutFreq ...
                   '_ord' num2str(filtOrder)];
end

end

