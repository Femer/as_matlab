function extendedModel = tool_extendModel(model)

extendedModel.A = [model.A,                      model.B;
                   zeros(1, length(model.A)),    1];
    
extendedModel.B = [model.B;
                   1];
               
extendedModel.C = eye(3);

extendedModel.Dt = model.Dt;

end

