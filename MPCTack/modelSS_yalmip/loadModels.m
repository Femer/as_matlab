function [realModel, uncertModel, uncertMPCModel] = loadModels(realModelType, factorSampleTime)

%model used in the kalman filter and in the MPC and LQR definition
uncertModel = load('idModel2015-02-18');
uncertModel = uncertModel.model;

display('---------- Info ----------');

if(strcmp(realModelType, 'little'))
    %load identified model with a and b
    load('linModelScalar');
    realModel = linModelScalar;
    display('Model with a and b.');
    %select wich  model you want to use
    nameModel = 'tack8';%7
else
    %load identified model with A and B
    load('linModelFull');
    realModel = linModelFull;
    display('Model with A and B.');
    %select wich  model you want to use
    nameModel = 'tack6';
end

%take the linear model to use in simulation
eval(['realModel = realModel.' nameModel ';']);

%change the sample time of the uncertModel, this undersampled model will be used
%be the MPC
oldDt = uncertModel.Dt;
uncertMPCModel = tool_changeModelSampleTime(uncertModel, factorSampleTime);
display(['Model for the MPC, sample time changed by a factor of ' num2str(factorSampleTime) ...
    ': from ' num2str(oldDt) ' to ' num2str(uncertMPCModel.Dt) ' [sec].']);

end

