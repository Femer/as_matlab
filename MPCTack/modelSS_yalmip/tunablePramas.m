%total simulation time, in seconds
tF = 10;

%prediction horizon of the MPC
predHor = 20;

%simulate using real state or estimated state by K.F. ?
useRealState = 0;

% Weights 
qYawRate = 1;
qYaw = 3;
rU = 1;%0.5
sChattering = 35;

%constraints
rudderMax = 0.75;
rudderVelocity = 1 / 0.5; %command/sec: rudder can go full right to full lest in 0.5 sec

%delta values to define the band near the origin
deltaYR = 7;
deltaY = 10;
deltaRud = 0.15;

%gaussian noise on measurements
varYawRate = 10 * pi / 180;
varYaw = 2 * pi / 180;
varRudder = 0.01;

%covariance matrices for Kalman filter
convarianceStr.R1 = 0.0005 * eye(3); %noise on the model
convarianceStr.R2 = blkdiag(varYawRate, varYaw, varRudder); %noisy measurements

%initial haul
tack = 'p2s'; %p2s s2p
              
%absDeltaYaw = 10 * pi / 180;

%% other param based on tunable parameters

absAlphaNew = 45 * pi / 180;

%matrices
Q = blkdiag(qYawRate, qYaw, rU);
R = sChattering;

%rudder value before starting the tack
rudderBeforeTack = 0; %between -1 and 1

%take the sample time of the real model, in seconds
meanTsSec = realModel.Dt;

%take the sample time of the model used by MPC or LQR
meanTsSecDown = uncertMPCModel.Dt;

%total simulation steps to reach tF
N = round(tF / meanTsSec);

%convert rudderVelocity from command/sec to command/simulationStep

%rudder velocity based on simulation time of the model used by MPC
rudderVelMPC = rudderVelocity * meanTsSecDown;

%rudder velocity in simulation time
rudderVelSim = rudderVelocity * meanTsSec;

%mesurements noise
measNoise = [sqrt(varYawRate) * randn(1, N);
             sqrt(varYaw) * randn(1, N);
             sqrt(varRudder) * randn(1, N)];
         
 display(['Horizon MPC: ' num2str(predHor * uncertMPCModel.Dt) ...
     ' [sec]; number of prediction steps: ' num2str(predHor) '.']);
 display('--------------------------');
 
 %compute yawReference based on type of tack
if(strcmp(tack, 'p2s'))
    alpha0 = 45 .* pi / 180;
    yawRef = -absAlphaNew - alpha0;
else
    alpha0 = -45 .* pi / 180;
    yawRef = +absAlphaNew - alpha0;
end

%usefull index
yawRateIndex = 1;
yawIndex = 2;
%another useful index for the extended state
lastRudderIndex = 3;