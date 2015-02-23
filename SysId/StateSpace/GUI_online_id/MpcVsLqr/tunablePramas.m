

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