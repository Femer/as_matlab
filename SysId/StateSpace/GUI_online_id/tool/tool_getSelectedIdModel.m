function [error, modelSelected] = tool_getSelectedIdModel(handles)

%see which model has been selected
contents = cellstr(get(handles.p_idModels,'String'));
nameModel = contents{get(handles.p_idModels,'Value')};
%make sure nameModel ~= identified model
if(strcmp(nameModel, 'identified models') ~= 1) 
    eval(['modelSelected = handles.idModels.' nameModel ';']);
    error = 0;
else
    msgbox('Please select an identified model', 'Error','error');
    error = 1;
end


end

