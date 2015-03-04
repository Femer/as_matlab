function idModels = tool_loadDefaultModels()

defaultModels = load('defaultIdentifiedModels');

%default_model_1 identified from tack8 from test on 21/01/2015 using d2d sampling =
%10 and black box model
idModels.default_model_1 = defaultModels.model1;

%default_model_2 identified from tack5B from test on 21/01/2015 using d2d sampling =
%10 and grey box model
idModels.default_model_2 = defaultModels.model2;

%default_model_3 identified from cross18 from test on 18/02/2015 using d2d sampling =
%10 and black box model
idModels.default_model_3 = defaultModels.model3;
end

