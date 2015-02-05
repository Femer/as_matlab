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
typeOfModel = 'capital'; %little or capital

if(strcmp(typeOfModel, 'little'))
    %load identified model with scalar value for a and b
    load('linModelScalar');
    linModel = linModelScalar;
    display('Model estimated with a and b');
else
    load('linModelFull');
    linModel = linModelFull;
    display('Model estimated with A and B');
end

%%

%downsample data
groupNSamples = 1;

%do you want to plot yawRate or yaw ?
choosenOutput = 'yaw';

%select the sequence from which you want to load the identified model
nameSeqId = 'tack6';
display(['Model estimated with data from ' nameSeqId '; estimated using ' typeOfModel]);
eval(['model = linModel.' nameSeqId ';']);

%use every sequence in steptacks for validation
seqNames = fieldnames(stepTacks);
seqNumb = length(seqNames);

tool_crossValidation(stepTacks, model, nameSeqId, choosenOutput, groupNSamples);







