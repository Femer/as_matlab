function [retVal, model] = tool_idModel(typeModel, timeSelected, logStr)

retVal = 0;

%convert timeSelected from seconds to milliseconds
timeSelected = timeSelected .* 1e3;

%cut logStr based on timeSelected(1) and timeSelected(2)
index = tool_selectIndexOnTime(timeSelected(1), timeSelected(2), logStr.time);

%debug
assignin('base', 'logStr', logStr);
assignin('base', 'timeSelected', timeSelected);
assignin('base', 'index', index);
    
time = logStr.time(index);
meanTs = mean(diff(time)); %in milliseconds

dataId.dataYawRate = iddata(logStr.yawRate(index), ...
                            logStr.rudder(index),...
                            meanTs);
                        
dataId.dataYaw = iddata(logStr.yaw(index), ...
                        logStr.rudder(index),...
                        meanTs);
dataId.time_sysId = time;

%see which type of model we have to use
if typeModel == 1
    [A, B, Dt] = tool_computeBlackBoxModel(dataId);
else
    [aStar, bStar, Dt] = tool_computeGreyBoxModel(dataId);
    A = [aStar,  0;
         Dt,     1];
     
    B = [bStar;
         0];
end

model.A = A;
model.B = B;
model.Dt = Dt;

end

