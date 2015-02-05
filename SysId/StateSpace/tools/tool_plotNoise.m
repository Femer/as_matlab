function [] = tool_plotNoise(validationData, model, nameModel, output)

nameValData = fieldnames(validationData);
numberValData = length(nameValData);


%simulate system responde for every validation data
for indexData = 1 : numberValData
   eval(['valData = validationData.' nameValData{indexData} ';']); 
   
   [yModel, wModel] = tool_computeScalarModelResponse(model, valData);
   
   tempOut.yModel = yModel;
   tempOut.wModel = wModel;
   
   
   eval(['modelOutput.' nameValData{indexData} ' = tempOut;']);
end

%plot
numRow = ceil(numberValData / 2);
numCol = 2;

if(strcmp(output, 'yaw'))
    leg = {'\psi noise'};
else
    leg = {'w noise'};
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
    noise = realOut - preictedOutput;
    
    subplot(numRow, numCol, index);
    plot(noise .* 180 / pi, 'm--', 'LineWidth', 1.3);
    
    grid on;
    legend(leg);
    title(nameValData{index});
    
    index = index + 1;
    if(index <= numberValData)
        
        if(strcmp(output, 'yaw'))
            eval(['preictedOutput = modelOutput.' nameValData{index} '.yModel;']);
            eval(['realOut = validationData.' nameValData{index} '.dataYaw.OutputData;']);
        else
            eval(['preictedOutput = modelOutput.' nameValData{index} '.wModel;']);
            eval(['realOut = validationData.' nameValData{index} '.dataYawRate.OutputData;']);
        end
        noise = realOut - preictedOutput;
    
        subplot(numRow, numCol, index);
        plot(noise .* 180 / pi, 'm--', 'LineWidth', 1.3);

        grid on;
        legend(leg);
        title(nameValData{index});
    end
end

end

