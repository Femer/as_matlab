function initParam = loadInitParams(realModel, typeTack)

%% tunable 
%total simulation time, in seconds
tF = 10;

%gaussian noise on measurements
varYawRate = 10 * pi / 180;
varYaw = 2 * pi / 180;
varRudder = 0.001;

%covariance matrices for Kalman filter
noiseModel = 0.0005;

%covariance matrix P{k-1}
noiseP11 = 0.5;

%rudder value before starting the tack
rudderBeforeTack = 0; %between -1 and 1

%absolute value of alpha star
absAlphaNew = 45 * pi / 180;

%physical constraints
rudderMax = 1;
rudderVel_cmd_sec = 4;

%% UNtunable

%constraints
initParam.rudderMax = rudderMax;
initParam.rudderVel_cmd_sec = rudderVel_cmd_sec;
%noise statistics
initParam.convarianceStr.R1 = noiseModel * eye(3); %noise on the model
initParam.convarianceStr.R2 = blkdiag(varYawRate, varYaw, varRudder); %noisy measurements

%compute real yaw0 based in type of tack
if(typeTack == 1)
   %port 2 starboard haul
   yaw0 = 2 * absAlphaNew;
else
    %starboard 2 port haul
   yaw0 = -2 * absAlphaNew;
end

%total simulation steps fo the real model to reach tF
initParam.N = round(tF / realModel.Dt);
initParam.realModelDt = realModel.Dt;
%mesurements noise
initParam.measNoise = [sqrt(varYawRate) * randn(1, initParam.N);
                       sqrt(varYaw) * randn(1, initParam.N);
                       sqrt(varRudder) * randn(1, initParam.N)];

initParam.varYawRate = varYawRate;
initParam.varYaw = varYaw;
initParam.varRudder = varRudder;

%guess on the initial state of the KF
% initParam.guessX1Hat = [  0 + sqrt(varYawRate) * randn();
%                           yaw0 +  sqrt(varYaw) * randn();
%                           rudderBeforeTack +  sqrt(varRudder) * randn()];

%guess on the variance matrix of the KF
initParam.guessP1_1 = noiseP11 * eye(3);

%start of the realModel used to simulate the boat in the C.L. MPC
initParam.xHatReal0 = [ 0;
                        yaw0;
                        rudderBeforeTack];

end