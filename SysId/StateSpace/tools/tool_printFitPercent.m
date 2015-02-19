function [numericModelFit] = tool_printFitPercent(listValidationData, model, nameModel, timeSampleRealSys)



nameValidationData = fieldnames(listValidationData);
numberValidationData = length(nameValidationData);


numericModelFit = [];

modelsFit = {'valData'; 'fit'};

%for every validation data in listValidationData
for indexValData = 1 : numberValidationData
    
    eval(['valData = listValidationData.' ...
          nameValidationData{indexValData} ';']);
           
    realOut = valData.dataYaw.OutputData;

      
    modelsFit{1, indexValData+1} = nameValidationData{indexValData};
    %simulate model response
    [modelYaw, ~] = tool_computeModelResponse(model, valData, timeSampleRealSys);
       
    %limit modelsOutput to be in [-pi, pi]
    indexOverPi = modelYaw > pi;
    indexLessMinusPi = modelYaw < - pi;
    
    modelYaw(indexOverPi) = modelYaw(indexOverPi) - 2 * pi;
    modelYaw(indexLessMinusPi) = modelYaw(indexLessMinusPi) + 2 * pi;
    
    fitPerc = 100 * (1 - (norm(realOut - modelYaw)) / (norm(realOut - mean(realOut))));
    
    modelsFit{2, indexValData + 1} = num2str(fitPerc, '%2.2f');
    
    numericModelFit(indexValData) = fitPerc;

end


%print
display([nameModel ' used to fit the model.']);
display(modelsFit);


end

