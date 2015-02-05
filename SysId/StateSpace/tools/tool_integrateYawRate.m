function [h_fig] = tool_integrateYawRate(yawValidationDataId, time, listEstTFs, ...
                                    nameList, nameValidationData, useJustNumb)

nameTFs = fieldnames(listEstTFs);
numberTFs = length(nameTFs);

modelsOutput = zeros(numberTFs, length(yawValidationDataId.OutputData));

%simulate every tf in listEstTFs
for indexTF = 1 : numberTFs
   eval(['myTf = listEstTFs.' nameTFs{indexTF} ';']); 
   
   myTfOut = sim(myTf, yawValidationDataId);
   
   modelsOutput(indexTF, :) = (myTfOut.OutputData)';
end

realOut = yawValidationDataId.OutputData;

%convert time in seconds
time = time ./ 1e6;
rel_time = (time - time(1));

%integrate modelsOutput using cumtrapz and time and add initlia yaw
for i = 1 : size(modelsOutput, 1)
    modelsOutput(i, :) = cumtrapz(time, modelsOutput(i, :)) + realOut(1);
    
    %limit modelsOutput to be in [-pi, pi]
    indexOverPi = modelsOutput(i, :) > pi;
    indexLessMinusPi = modelsOutput(i, :) < - pi;
    
    modelsOutput(i, indexOverPi) = modelsOutput(i, indexOverPi) - 2 * pi;
    modelsOutput(i, indexLessMinusPi) = modelsOutput(i, indexLessMinusPi) + 2 * pi;
end

%plot
numRow = ceil(numberTFs / 2);
numCol = 2;


h_fig = figure;
set(gcf,'name', ...
    ['Tfs estimated from: ' nameList '; validation data from ' nameValidationData], ...
    'numbertitle', 'off');

leg = {[yawValidationDataId.OutputName{1} ' real'],...
       [yawValidationDataId.OutputName{1} ' model']};
if(isempty(useJustNumb))
    for i = 1 : numRow
        %plot results 
        index = (i-1) * 2 + 1;
        subplot(numRow, numCol, index);
        plot(rel_time, realOut .* 180 / pi, 'b');
        hold on;
        plot(rel_time, modelsOutput(index,:) .* 180 / pi, 'm--', 'LineWidth', 1.3);

        grid on;
        legend(leg, 'Location', 'northwest');
        title(nameTFs{index});
        xlabel('Time [sec]');
        ylabel('deg');

        index = index + 1;
        if(index <= numberTFs)
            subplot(numRow, numCol, index);
            plot(rel_time, realOut .* 180 / pi, 'b');
            hold on;
            plot(rel_time, modelsOutput(index,:) .* 180 / pi, 'm--', 'LineWidth', 1.3);

            grid on;
            legend(leg, 'Location', 'northwest');
            title(nameTFs{index});
            xlabel('Time [sec]');
            ylabel('deg');
        end
    end
else
    plot(rel_time, realOut .* 180 / pi, 'b', 'LineWidth', 1.45);
    hold on;
    plot(rel_time, modelsOutput(useJustNumb,:) .* 180 / pi, 'm--', 'LineWidth', 1.5);
    
    grid on;
    legend(leg, 'Location', 'northwest');
    xlabel('Time [sec]');
    ylabel('deg');
end

end

