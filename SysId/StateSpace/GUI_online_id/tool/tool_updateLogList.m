function tool_updateLogList(handles)

%debug
assignin('base', 'handles', handles);
    
%add every log in handles.logs into the String field
logList = fieldnames(handles.logs);
%add strin 'log list'
logList = ['log list'; logList];

set(handles.p_logList, 'String', logList);

end

