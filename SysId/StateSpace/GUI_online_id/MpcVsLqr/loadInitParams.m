function initParam = loadInitParams(realModel, typeTack)

%% tunable 
%total simulation time, in seconds
tF = 10;

%gaussian noise on measurements
varYawRate = 10 * pi / 180;
varYaw = 2 * pi / 180;
varRudder = 0.01;

%covariance matrices for Kalman filter
noiseModel = 0.0005;

%covariance matrix P{k-1}
noiseP11 = 1;

%rudder value before starting the tack
rudderBeforeTack = 0; %between -1 and 1

%absolute value of alpha star
absAlphaNew = 45 * pi / 180;

%% UNtunable

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

%mesurements noise
initParam.measNoise = [sqrt(varYawRate) * randn(1, N);
                       sqrt(varYaw) * randn(1, N);
                       sqrt(varRudder) * randn(1, N)];

%total simulation steps fo the real model to reach tF
initParam.N = round(tF / realModel.Dt);

%guess on the initial state of the KF
initParam.guessX1Hat = [  0 + sqrt(varYawRate) * randn();
                          yaw0 +  sqrt(varYaw) * randn();
                          rudderBeforeTack +  sqrt(varRudder) * randn()];

%guess on the variance matrix of the KF
initParam.guessP1_1 = noiseP11 * blkdiag(eye(2), 0);

%start of the realModel used to simulate the boat in the C.L. MPC
initParam.xHatReal0 = [ 0;
                        yaw0;
                        rudderBeforeTack];

end