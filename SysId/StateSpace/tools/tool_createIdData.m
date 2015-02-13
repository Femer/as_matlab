function [structOut] = tool_createIdData(tack)

%compute mean sample time in the data, in uSec
meanTs = mean(diff(tack.time_sysId));
%create data for yawRate identification
dataYawRate = iddata(tack.yawRate_sysId, tack.rudder_sysId, meanTs);
dataYawRate.OutputName = 'yawRate';
dataYawRate.InputName = 'rudder';
dataYawRate.InputUnit = 'Pixhawk_cmd';
dataYawRate.OutputUnit = 'rad/s';
dataYawRate.TimeUnit = 'microseconds';

%create data for yaw validation
dataYaw = iddata(tack.yaw_sysId, tack.rudder_sysId, meanTs);
dataYaw.OutputName = 'yaw';
dataYaw.InputName = 'rudder';
dataYaw.InputUnit = 'Pixhawk_cmd';
dataYaw.OutputUnit = 'rad';
dataYaw.TimeUnit = 'microseconds';

structOut.dataYawRate = dataYawRate;
structOut.dataYaw = dataYaw;
structOut.time_sysId = tack.time_sysId;

end

