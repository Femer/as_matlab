% Based on weights specified, compute LQR gain and MPC Hessian H matrix

clc;
clear;
close all;

typeOfModel = 'little';%little (a,b) or capital (A,B)

%change sampling time of the model used by the MPC
factorSamplingTimeMPC = 10;

% Weights 
qYawRate = 0.00001;
qYaw = 10;
rU = 0.001;
sChattering = 5;

%constraints
rudderMax = 0.9;%0.9;
rudderVelocity = 1.8 / 0.5; %command/sec: rudder can go full right to full lest in 0.5 sec

%matrices
Q = blkdiag(qYawRate, qYaw, rU);
R = sChattering;


if(strcmp(typeOfModel, 'little'))
    %load identified model with a and b
    load('linModelScalar');
    linModel = linModelScalar;
    nameModelUsed = 'a and b';
    %select wich  model you want to use
    nameModel = 'tack8';
else
    %load identified model with A and B
    load('linModelFull');
    linModel = linModelFull;
    nameModelUsed = 'A and B';
    %select wich  model you want to use
    nameModel = 'tack6';
end


%take the linear model to use in simulation
eval(['model = linModel.' nameModel ';']);

% extended model xHat = [yawRate_k, yaw_k, rudder_{k-1}],
% uHat = [rudder_{k} - rudder{k-1}];
% we use the extended model to achieve the same cost function in the
% LQR and in the MPC. Since the LQR can have a cost function of the form
% x' * Q * x + u' * R * u, that tries to bring the system to the origin

%LQR: use normal system sampling time
AExt = [model.A,                      model.B;
        zeros(1, length(model.A)),    1];
    
BExt = [model.B;
        1];
    
%compute gain matrix for the LQR and the solution of the discrete riccati
%equation
[K_LQR, ~, ~] = dlqr(AExt, BExt, Q, R); 

%MPC: change model sampling time by a factor of factorSamplingTimeMPC
modelMPC = tool_changeModelSampleTime(model, factorSamplingTimeMPC);

AExt = [modelMPC.A,                      modelMPC.B;
        zeros(1, length(modelMPC.A)),    1];
    
BExt = [modelMPC.B;
        1];
%discrete riccati solution
M = dare(AExt, BExt, Q, R);

%MPC forces Hesian matrix
H = [sChattering; qYawRate; qYaw; rU];

%MPC forces Hesian matrix for final cost
H_final = blkdiag(0, M);


%every simulation step lasts meanTsSec seconds.
rudderVelSim = rudderVelocity * modelMPC.Dt;
%MPC lower and upper bound
lowerBound = [-rudderVelSim; -rudderMax];
upperBound = [rudderVelSim; rudderMax];


%print matrices and info
display('---------- Info ----------');
display(['Model with estimated ' nameModelUsed '.']);
display(['LQR model sampling time: ' num2str(model.Dt) ' [sec].']);
display(['MPC model sampling time: ' num2str(modelMPC.Dt) ' [sec].']);
display('--------------------------');
printmat(K_LQR, 'LQR gain', '', 'ASO_LQR_K1 ASO_LQR_K2 ASO_LQR_K3');

printmat(H', 'MPC H', '', 'ASO_MPC_H1 ASO_MPC_H2 ASO_MPC_H3 ASO_MPC_H4');

printmat(lowerBound', 'MPC lower bound', '', 'ASO_MPC_LB1 ASO_MPC_LB2');

printmat(upperBound', 'MPC uppper bound', '', 'ASO_MPC_UB1 ASO_MPC_UB2');

printmat(H_final(2, 2:end), 'MPC H final (2, 2:end)', '', ...
        'ASO_MPC_HF22 ASO_MPC_HF23 ASO_MPC_HF24');
    
printmat(H_final(3, 2:end), 'MPC H final (3, 2:end)', '', ...
        'ASO_MPC_HF32 ASO_MPC_HF33 ASO_MPC_HF34');
    
printmat(H_final(4, 2:end), 'MPC H final (4, 2:end)', '', ...
        'ASO_MPC_HF42 ASO_MPC_HF43 ASO_MPC_HF44');
