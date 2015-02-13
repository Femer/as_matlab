function [newModel] = tool_changeModelSampleTime(model, factor)

% newModel.Dt = model.Dt * factor;
% newModel.A = (model.A) ^ factor;
% 
% powerA = zeros(size(model.A));
% 
% for i = 0 : factor-1
%    powerA = powerA + powerA ^i;
% end
% 
% newModel.B = powerA * model.B;

appS = ss(model.A, model.B, [], [], 1);
appN = d2d(appS, factor);

newModel.A = appN.a;
newModel.B = appN.b;
newModel.Dt = model.Dt * factor;

end

