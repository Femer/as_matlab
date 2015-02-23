function predHor_s = tool_getPredHorSec(handles, model)

%based on the prediction horizon, in steps, and the model sampling time,
%compute the prediction horizon in seconds
contents = cellstr(get(handles.p_predHorizon,'String')); 
strPredHor = contents{get(handles.p_predHorizon,'Value')};

predHor_s = str2double(strPredHor) * model.Dt;

end

