clc;
clear;
%close all;

%load step data
%extractStepsFromHelmsman;

%load full tack data
%extractTacksFromHelmsman;

%load controlled course
%extractControlledCourse;

%load cross validation data from test on 2015-02-10
load('dataCrossValidation2015-02-10');

%tool
addpath('../tools/');

%Update the model state with the real state every secResampleState sec
secResampleState = 1;

%estimated model with a and B or with A and B
typeOfModel = 'capital'; %little or capital

%choose if you want to increase the sample time of the model you selected
factorSampleTime = 1;

if(strcmp(typeOfModel, 'little'))
    %load identified model with scalar value for a and b
    load('linModelScalar');
    linModels = linModelScalar;
    display('Model estimated with a and b');
    %select the sequence from which you want to load the identified model
    nameSeqId = 'tack8';
    dateSysId = '21-01-2015';
else
    load('linModelFull');
    linModels = linModelFull;
    display('Model estimated with A and B');
    %select the sequence from which you want to load the identified model
    nameSeqId = 'tack6';
    dateSysId = '21-01-2015';
end


%%

%do you want to plot yawRate or yaw ?
choosenOutput = 'yawRate';

display(['Model estimated with data from ' nameSeqId '; estimated using ' typeOfModel ' model.']);
eval(['model = linModels.' nameSeqId ';']);

%change the sample time of the model 
if(factorSampleTime > 1)
    oldDt = model.Dt;
    model = tool_changeModelSampleTime(model, factorSampleTime);
    display(['Model sample time changed by a factor of ' num2str(factorSampleTime) ...
        ': from ' num2str(oldDt) ' to ' num2str(model.Dt) ' [sec].']);
end

timeSampleRealSys = round(secResampleState / model.Dt);
display(['Model state updated with real state every ' ...
        num2str(secResampleState) ' sec.']);

%use every sequence in steptacks for validation
seqNames = fieldnames(stepTacks);
seqNumb = length(seqNames);

tool_crossValidation(stepTacks, model, nameSeqId, choosenOutput, ...
                     factorSampleTime, dateSysId, timeSampleRealSys);







