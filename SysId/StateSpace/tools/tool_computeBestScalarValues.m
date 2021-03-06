function [aStar, bStar, Dt] = tool_computeBestScalarValues(dataId)

%take yaw rate
W = dataId.dataYawRate.OutputData;
%take yaw angle
Y = dataId.dataYaw.OutputData;
%take rudder commands
U = dataId.dataYawRate.InputData;

%backward compatibility
if(isfield(dataId, 'bgud_time'))
    time = dataId.bgud_time;
else
    time = dataId.time_sysId;
end
%compute average delta dt over the sampled measurements
Dt = mean(diff(time)); %mean deltaT in uSec
%convert Dt from uSec to sec
Dt = Dt / 1e6;

%number of measurements
N = length(W);

%create matrix for pseudoinverse computation
A = [   W(1 : N-1),     U(1 : N-1);
        W(1 : N-2),     U(1 : N-2)];
    
b = [   W(2 : N);
        (1 / Dt) * (Y(3 : N) - Y(2 : N-1))];
    
%compute pseudoinverse
x = A \ b;

%take result
aStar = x(1);
bStar = x(2);
end

