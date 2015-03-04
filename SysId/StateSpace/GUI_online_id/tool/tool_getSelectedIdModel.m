function [error, modelSelected, indexModel, nameModel] = tool_getSelectedIdModel(handles)

modelSelected = [];
%see which model has been selected
contents = cellstr(get(handles.p_idModels,'String'));
indexModel = get(handles.p_idModels,'Value');
nameModel = contents{indexModel};
%make sure indexModel ~= 1
if(indexModel ~= 1) 
    eval(['modelSelected = handles.idModels.' nameModel ';']);
    error = 0;
else
    msgbox('Please select an identified model', 'Error','error');
    error = 1;
end

end

