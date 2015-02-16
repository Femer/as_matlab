clc;
clear;
close all;

addpath('../tools/');

%load step data from takcs
load('dataIDStepTacks2015-02-10');

%% ID

%for every tack in stepTacks, identify ARX models
tackNames = fieldnames(stepTacks);
tackNumb = length(tackNames);

for i = 1 : tackNumb
   
    eval(['seqId = stepTacks.' tackNames{i} ';']);

    %least square solution
    [A, B, Dt] = tool_computeBestFullAB(seqId);
    model.A = A;
    model.B = B;
    model.Dt = Dt;
    
    eval(['linModelFull.' tackNames{i} ' = model;']);
end

%% save estimated transfer functions
save('linModelFull', 'linModelFull');
