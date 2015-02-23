function [error, logStr] = tool_getSelectedLog(handles)

logStr = [];
%take the selected log and use it as validation data
contents = cellstr(get(handles.p_logList,'String')); %returns p_logList contents as cell array
indexPopUp = get(handles.p_logList,'Value');
selectedLog = contents{indexPopUp}; %returns selected item from p_logList

%check
if(indexPopUp == 1)
    msgbox('Please select a log to use as validation data', 'Error','error');
    error = 1;
else
    eval(['logStr = handles.logs.' selectedLog ';']);
    %everything ok
    error = 0;
end


end

