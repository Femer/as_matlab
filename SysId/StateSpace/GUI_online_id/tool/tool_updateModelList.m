function tool_updateModelList(handles)

%add every model in handles.idModels into the String field
modelList = fieldnames(handles.idModels);

%string for p_idModels
set(handles.p_idModels, 'String', ['identified models'; modelList]);

%string for p_realModel
set(handles.p_realModel, 'String', ['"real" model in simulation'; modelList]);

%string for p_lqrModel
set(handles.p_lqrModel, 'String', ['LQR model'; modelList]);

%string for p_mpcModel
set(handles.p_mpcModel, 'String', ['MPC model'; modelList]);

end

