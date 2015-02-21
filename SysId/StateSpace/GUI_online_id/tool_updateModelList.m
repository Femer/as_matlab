function tool_updateModelList(handles)

%add every model in handles.idModels into the String field
modelList = fieldnames(handles.idModels);
%add string 'log list'
modelList = ['idetentified models' modelList];

set(handles.p_idModels, 'String', modelList);

end

