function [structOut] = tool_createIdData(tack)

%compute mean sample time in the data, in uSec
meanTs = mean(diff(tack.bgud_time));
%create data for yawRate identification
dataYawRate = iddata(tack.yawRate_interp, tack.rudder, meanTs);
dataYawRate.OutputName = 'yawRate';
dataYawRate.InputName = 'rudder';
dataYawRate.InputUnit = 'Pixhawk_cmd';
dataYawRate.OutputUnit = 'rad/s';
dataYawRate.TimeUnit = 'microseconds';

%create data for yaw validation
dataYaw = iddata(tack.yaw_interp, tack.rudder, meanTs);
dataYaw.OutputName = 'yaw';
dataYaw.InputName = 'rudder';
dataYaw.InputUnit = 'Pixhawk_cmd';
dataYaw.OutputUnit = 'rad';
dataYaw.TimeUnit = 'microseconds';

structOut.dataYawRate = dataYawRate;
structOut.dataYaw = dataYaw;
structOut.bgud_time = tack.bgud_time;

end

