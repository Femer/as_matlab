function [] = tool_crossValidation(validationData, model, nameModel, ...
                                   output, groupNSamples, dateSysId, timeSampleRealSys)

nameValData = fieldnames(validationData);
numberTotalValData = length(nameValData);


%simulate system responde for every validation data
for indexData = 1 : numberTotalValData
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
   [yModel, wModel] = tool_computeModelResponse(model, inputModel, timeSampleRealSys);
   
   tempOut.yModel = yModel;
   tempOut.wModel = wModel;
   tempOut.time_sysId = inputModel.time_sysId;
   
   %save output data and time
   eval(['modelOutput.' nameValData{indexData} ' = tempOut;']);
end



if(strcmp(output, 'yaw'))
    leg = {'\psi real', '\psi model'};
else
    leg = {'w real', 'w model'};
end

%in each figure use no more than 3 rows
totalRows = ceil(numberTotalValData / 2);
numbFig = ceil(totalRows / 3);

%take the index to access to the full data, regardless of which
%figure we are now plotting
indexData = 1;

for indexFig = 1 : numbFig

    figure;
    set(gcf,'name', ...
        ['Model estimated from: ' nameModel ' (' dateSysId '); figure ' num2str(indexFig)], ...
        'numbertitle', 'off');
    
    %number of validation data to plot in this figure
    if((numberTotalValData - ((indexFig - 1) * 6)) >= 6)
        valDataFig = 6;
    else
        valDataFig = numberTotalValData - ((indexFig - 1) * 6);
    end
    
    %in this figure, how many rows do I have to make?
    numRow = ceil(valDataFig / 2);
    numCol = 2;
    
    indexSubplot = 1;
    
    for i = 1 : numRow
        
        if(strcmp(output, 'yaw'))
            eval(['preictedOutput = modelOutput.' nameValData{indexData} '.yModel;']);
            eval(['realOut = validationData.' nameValData{indexData} '.dataYaw.OutputData;']);
        else
            eval(['preictedOutput = modelOutput.' nameValData{indexData} '.wModel;']);
            eval(['realOut = validationData.' nameValData{indexData} '.dataYawRate.OutputData;']);
        end
        
        %time vectors
        eval(['timeModelData = modelOutput.' nameValData{indexData} '.time_sysId;']);
        eval(['timeRealData = validationData.' nameValData{indexData} '.time_sysId;']);
        timeModelData = (timeModelData - timeModelData(1)) ./ 1e6;
        timeRealData = (timeRealData - timeRealData(1)) ./ 1e6;
        
        subplot(numRow, numCol, indexSubplot);
        plot(timeRealData, realOut .* 180 / pi, 'b', 'LineWidth', 1.8);
        hold on;
        plot(timeModelData, preictedOutput .* 180 / pi, 'm--', 'LineWidth', 1.9);
        
        grid on;
        legend(leg);
        title(nameValData{indexData});
        xlabel('Time [sec]');
        
        indexData = indexData + 1;
        indexSubplot = indexSubplot + 1;
        if(indexSubplot <= valDataFig)
            
            if(strcmp(output, 'yaw'))
                eval(['preictedOutput = modelOutput.' nameValData{indexData} '.yModel;']);
                eval(['realOut = validationData.' nameValData{indexData} '.dataYaw.OutputData;']);
            else
                eval(['preictedOutput = modelOutput.' nameValData{indexData} '.wModel;']);
                eval(['realOut = validationData.' nameValData{indexData} '.dataYawRate.OutputData;']);
            end
            %time vectors
            eval(['timeModelData = modelOutput.' nameValData{indexData} '.time_sysId;']);
            eval(['timeRealData = validationData.' nameValData{indexData} '.time_sysId;']);
            timeModelData = (timeModelData - timeModelData(1)) ./ 1e6;
            timeRealData = (timeRealData - timeRealData(1)) ./ 1e6;
            
            subplot(numRow, numCol, indexSubplot);
            plot(timeRealData, realOut .* 180 / pi, 'b', 'LineWidth', 1.8);
            hold on;
            plot(timeModelData, preictedOutput .* 180 / pi, 'm--', 'LineWidth', 1.9);
            
            grid on;
            legend(leg);
            title(nameValData{indexData});
            xlabel('Time [sec]');
            
            indexData = indexData + 1;
            indexSubplot = indexSubplot + 1;
        end
    end

end

end

