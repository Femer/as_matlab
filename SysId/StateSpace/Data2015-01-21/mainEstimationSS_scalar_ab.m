clc;
clear;
close all;

%load step data from takcs
extractStepsFromHelmsman;

%% ID


%for every takc in stepTacks, identify ARX models
tackNames = fieldnames(stepTacks);
tackNumb = length(tackNames);

for i = 1 : tackNumb
   
    eval(['seqId = stepTacks.' tackNames{i} ';']);

    %least square solution
    [a, b, Dt] = tool_computeBestScalarValues(seqId);
    model.A = [a,   0;
               Dt,  1];
           
    model.B = [b;
               0];
           
    model.Dt = Dt;
    
    eval(['linModelScalar.' tackNames{i} ' = model;']);
end

%% save estimated transfer functions
save('linModelScalar', 'linModelScalar');
