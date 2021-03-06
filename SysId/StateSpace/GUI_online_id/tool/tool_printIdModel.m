function tool_printIdModel(handles, model, indexModel)


%print matrices A and B. Print sample time of the model in milliseconds!
s = {'A', ...
    num2str(model.A(1,:)), ...
    num2str(model.A(2,:)), ...
    ' ',...
    'B', ...
    num2str(model.B(1,:)), ...
    num2str(model.B(2,:)), ...
    ' ',...
    ['Dt: ' num2str(model.Dt * 1e3) ' mS']};%in milliseconds


set(handles.t_modelSelected, 'String', s);

%based on the prediction horizon, in steps, compute the prediction horizon
%in seconds
predHor_s = tool_getPredHorSec(handles, model);

%show it on screen
set(handles.t_predHorizon, 'String', ...
    ['prediction horizon: ' num2str(predHor_s) ' [sec].']);

%put the same model selected even for the MPC and LQR model
set(handles.p_lqrModel, 'Value', indexModel);
set(handles.p_mpcModel, 'Value', indexModel);

end

