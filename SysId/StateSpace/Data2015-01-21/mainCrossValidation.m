clc;
clear;
%close all;

%load step data
%extractStepsFromHelmsman;

%load full tack data
extractTacksFromHelmsman;

%load controlled course
%extractControlledCourse;

%load full tack from test 2015-01-15
%load('stepTacks2015-01-15');

%tool
addpath('../tools/');

%estimated model with a and B or with A and B
typeOfModel = 'little'; %little or capital

%choose if you want to increase the sample time of the model you selected
factorSampleTime = 10;

if(strcmp(typeOfModel, 'little'))
    %load identified model with scalar value for a and b
    load('linModelScalar');
    linModels = linModelScalar;
    display('Model estimated with a and b');
    %select the sequence from which you want to load the identified model
    nameSeqId = 'tack8';
else
    load('linModelFull');
    linModels = linModelFull;
    display('Model estimated with A and B');
    %select the sequence from which you want to load the identified model
    nameSeqId = 'tack6';
end


%%

%do you want to plot yawRate or yaw ?
choosenOutput = 'yaw';

display(['Model estimated with data from ' nameSeqId '; estimated using ' typeOfModel '.']);
eval(['model = linModels.' nameSeqId ';']);

%change the sample time of the model 
if(factorSampleTime > 1)
    oldDt = model.Dt;
    model = tool_changeModelSampleTime(model, factorSampleTime);
    display(['Model sample time changed by a factor of ' num2str(factorSampleTime) ...
        ': from ' num2str(oldDt) ' to ' num2str(model.Dt) ' [sec].']);
end

%use every sequence in steptacks for validation
seqNames = fieldnames(stepTacks);
seqNumb = length(seqNames);

tool_crossValidation(stepTacks, model, nameSeqId, choosenOutput, factorSampleTime);







