function [error, lqrModel, nameLqr, mpcModel, nameMpc] = tool_getLqrMpcModels(handles)

lqrModel = [];
mpcModel = [];
nameLqr = [];
nameMpc = [];
%hope for the best
error = 0;

%see which model has been selected
contents = cellstr(get(handles.p_lqrModel,'String'));
indexLqr = get(handles.p_lqrModel,'Value');
nameLqr = contents{indexLqr};

%make sure indexModel ~= 1
if(indexLqr ~= 1) 
    eval(['lqrModel = handles.idModels.' nameLqr ';']);
else
    msgbox('Please select a model to design the LQR controller', 'Error','error');
    error = 1;
    return;
end

%see which model has been selected
contents = cellstr(get(handles.p_mpcModel,'String'));
indexMpc = get(handles.p_mpcModel,'Value');
nameMpc = contents{indexMpc};

%make sure indexModel ~= 1
if(indexMpc ~= 1) 
    eval(['mpcModel = handles.idModels.' nameMpc ';']);
else
    msgbox('Please select a model to design the MPC controller', 'Error','error');
    error = 1;
    return;
end

end

