clc;
clear;
close all;

%load validation data from test on 10-02-2015
load('dataCrossValidation2015-02-10');

%tool
addpath('../tools/');

%estimated model with a and B or with A and B
typeOfModel = 'capital'; %little or capital

display('----------------------------');
if(strcmp(typeOfModel, 'little'))
    %load identified model with scalar value for a and b
    load('linModelScalar');
    linModel = linModelScalar;
    display('Model estimated with scalar a and b');
else
    load('linModelFull');
    linModel = linModelFull;
    display('Model estimated with full A and B');
end
display(['----------------------------' 10]);

nameModels = fieldnames(linModel);
numberModels = length(nameModels);

numericModelFit = zeros(numberModels, length(fieldnames(stepTacks)));

for i = 1 : numberModels
    eval(['model = linModel.' nameModels{i} ';']);
    numericModelFit(i, :) = tool_printFitPercent(stepTacks, model, nameModels{i});
end

%compute avg
outputAvgFit = {'models'; 'Avg fit'};

for i = 1 : numberModels
    outputAvgFit{1, i+1} = nameModels{i};
    avgVal = sum(numericModelFit(i, :)) / length(numericModelFit(i, :));
    outputAvgFit{2, i+1} = num2str(avgVal, '%2.2f');
end

display(outputAvgFit);
