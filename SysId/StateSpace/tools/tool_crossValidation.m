function [] = tool_crossValidation(validationData, model, nameModel, output, groupNSamples)

nameValData = fieldnames(validationData);
numberValData = length(nameValData);


%simulate system responde for every validation data
for indexData = 1 : numberValData
   eval(['valData = validationData.' nameValData{indexData} ';']); 
   
   inputModel = valData;
   %if groupNSamples > 1 that means that time sample of the model has been changed
   %using tool_changeModel.
   %To simulate the model response, use tool_downsampleData to extract downsampled data
   %with the same groupNSamples, but now DO NOT use a mean value within
   %every sequence made by groupNSamples samples.
   if(groupNSamples > 1)
       inputModel = tool_downsampleData(valData, groupNSamples);
   end
   [yModel, wModel] = tool_computeModelResponse(model, inputModel);
   
   tempOut.yModel = yModel;
   tempOut.wModel = wModel;
   tempOut.bgud_time = inputModel.bgud_time;
   
   %save output data and time
   eval(['modelOutput.' nameValData{indexData} ' = tempOut;']);
end

%plot
numRow = ceil(numberValData / 2);
numCol = 2;

if(strcmp(output, 'yaw'))
    leg = {'\psi real', '\psi model'};
else
    leg = {'w real', 'w model'};
end



figure;
set(gcf,'name', ...
    ['Model estimated from: ' nameModel], ...
    'numbertitle', 'off');


for i = 1 : numRow
    %plot results 
    index = (i-1) * 2 + 1;
    
    if(strcmp(output, 'yaw'))
        eval(['preictedOutput = modelOutput.' nameValData{index} '.yModel;']);
        eval(['realOut = validationData.' nameValData{index} '.dataYaw.OutputData;']);
    else
        eval(['preictedOutput = modelOutput.' nameValData{index} '.wModel;']);
        eval(['realOut = validationData.' nameValData{index} '.dataYawRate.OutputData;']);
    end
    
    %time vectors
    eval(['timeModelData = modelOutput.' nameValData{index} '.bgud_time;']);
    eval(['timeRealData = validationData.' nameValData{index} '.bgud_time;']);
    timeModelData = (timeModelData - timeModelData(1)) ./ 1e6;
    timeRealData = (timeRealData - timeRealData(1)) ./ 1e6;
    
    subplot(numRow, numCol, index);
    plot(timeRealData, realOut .* 180 / pi, 'b');
    hold on;
    plot(timeModelData, preictedOutput .* 180 / pi, 'm--', 'LineWidth', 1.3);
    
    grid on;
    legend(leg);
    title(nameValData{index});
    xlabel('Time [sec]');
    
    index = index + 1;
    if(index <= numberValData)
        
        if(strcmp(output, 'yaw'))
            eval(['preictedOutput = modelOutput.' nameValData{index} '.yModel;']);
            eval(['realOut = validationData.' nameValData{index} '.dataYaw.OutputData;']);
        else
            eval(['preictedOutput = modelOutput.' nameValData{index} '.wModel;']);
            eval(['realOut = validationData.' nameValData{index} '.dataYawRate.OutputData;']);
        end
        %time vectors
        eval(['timeModelData = modelOutput.' nameValData{index} '.bgud_time;']);
        eval(['timeRealData = validationData.' nameValData{index} '.bgud_time;']);
        timeModelData = (timeModelData - timeModelData(1)) ./ 1e6;
        timeRealData = (timeRealData - timeRealData(1)) ./ 1e6;
    
        subplot(numRow, numCol, index);
        plot(timeRealData, realOut .* 180 / pi, 'b');
        hold on;
        plot(timeModelData, preictedOutput .* 180 / pi, 'm--', 'LineWidth', 1.3);

        grid on;
        legend(leg);
        title(nameValData{index});
        xlabel('Time [sec]');
    end
end

end

