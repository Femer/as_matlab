function [index] = tool_selectIndexOnTime(tStart,tStop, time)

index = time >= tStart & time <= tStop;
end

