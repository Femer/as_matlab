function [A, B, Dt] = tool_computeBestFullAB(dataId)

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
%convert Dt from milliSec to sec
Dt = Dt / 1e3;

%number of measurements
N = length(W);

%create matrix for pseudoinverse computation
M1 = [   W(1 : N-1),     Y(1 : N-1),    U(1 : N-1)];
    
b1 = W(2 : N);

M2 = M1;
    
b2 = Y(2 : N);
    
%compute pseudoinverse
x1 = M1 \ b1;
x2 = M2 \ b2;

%take result
A = [x1(1),     x1(2);
     x2(1),     x2(2)];
 
B = [   x1(3);
        x2(3)];
end

