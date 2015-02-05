function [dataDownSampled] = tool_downsampleData(data, groupNSamples)

index = 1;

N = length(data.dataYaw.OutputData);

%downsample yaw, yawRate, rudder and time
yaw = data.dataYaw.OutputData;
yawRate = data.dataYawRate.OutputData;
rudder = data.dataYawRate.InputData;
time = data.bgud_time;

count = 1;

yawDown = [];
yawRateDown = [];
rudderDown = [];
timeDown = [];

while(index <= N - groupNSamples + 1)
   %take data every groupNSamples steps
   
   stop = index + groupNSamples - 1;
   
   %take avg of every samples between [index, stop], except for time
   yawDown(count) = yaw(index);
   yawRateDown(count) = yawRate(index);
   rudderDown(count) = rudder(index);
   timeDown(count) = time(index);
   
   count = count + 1;
   
   index = index + groupNSamples;
end

%compute new mean sampling time
meanTs = mean(diff(timeDown));

%create data for yawRate identification
dataYawRate = iddata(yawRateDown', rudderDown', meanTs);
dataYawRate.OutputName = 'yawRate';
dataYawRate.InputName = 'rudder';
dataYawRate.InputUnit = 'Pixhawk_cmd';
dataYawRate.OutputUnit = 'rad/s';
dataYawRate.TimeUnit = 'microseconds';

%create data for yaw validation
dataYaw = iddata(yawDown', rudderDown', meanTs);
dataYaw.OutputName = 'yaw';
dataYaw.InputName = 'rudder';
dataYaw.InputUnit = 'Pixhawk_cmd';
dataYaw.OutputUnit = 'rad';
dataYaw.TimeUnit = 'microseconds';

dataDownSampled.dataYawRate = dataYawRate;
dataDownSampled.dataYaw = dataYaw;
dataDownSampled.bgud_time = timeDown;

end

