function [realModel, uncertModel, uncertMPCModel] = loadModels(uncertModelType, factorSampleTime)

%model used to simulate the real boat
realModel = load('idModel2015-02-18');
realModel = realModel.model;

display('---------- Info ----------');

if(strcmp(uncertModelType, 'little'))
    %load identified model with a and b
    load('linModelScalar');
    uncertModel = linModelScalar;
    display('Model with a and b.');
    %select wich  model you want to use
    nameModel = 'tack8';%7
else
    %load identified model with A and B
    load('linModelFull');
    uncertModel = linModelFull;
    display('Model with A and B.');
    %select wich  model you want to use
    nameModel = 'tack6';
end

%take the linear model to build the LQR and the MPC
eval(['uncertModel = uncertModel.' nameModel ';']);

%change the sample time of the uncertModel, this undersampled model will be used
%be the MPC
oldDt = uncertModel.Dt;
uncertMPCModel = tool_changeModelSampleTime(uncertModel, factorSampleTime);
display(['Model for the MPC, sample time changed by a factor of ' num2str(factorSampleTime) ...
    ': from ' num2str(oldDt) ' to ' num2str(uncertMPCModel.Dt) ' [sec].']);

end

