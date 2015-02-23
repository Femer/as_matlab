function [error, modelSelected] = tool_getRealModel(handles)

modelSelected = [];
%see which model has been selected
contents = cellstr(get(handles.p_realModel,'String'));
indexModel = get(handles.p_realModel,'Value');
nameModel = contents{indexModel};
%make sure indexModel ~= 1
if(indexModel ~= 1) 
    eval(['modelSelected = handles.idModels.' nameModel ';']);
    error = 0;
else
    msgbox('Please select a real model to use in the simulation', 'Error','error');
    error = 1;
end

end

