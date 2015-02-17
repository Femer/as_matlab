% Based on weights specified, compute LQR gain and MPC Hessian H matrix
% and write a txt file cointaining all the parameters needed
%to use in QGC to send the settings to the boat.

clc;
clear;
close all;

typeOfModel = 'capital';%little (a,b) or capital (A,B)

%change sampling time of the model used by the MPC
factorSamplingTimeMPC = 10;

% Weights 
qYawRate = 1;
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

%text string for the .txt file
fileID = fopen('as_optimal_control_param.txt','w');
%LQR gain
fprintf(fileID, '1\t50\tASO_LQR_K1\t%0.10f\t9\n', K_LQR(1));
fprintf(fileID, '1\t50\tASO_LQR_K2\t%0.10f\t9\n', K_LQR(2));
fprintf(fileID, '1\t50\tASO_LQR_K3\t%0.10f\t9\n', K_LQR(3));
% 
% %H vector
fprintf(fileID, '1\t50\tASO_MPC_H1\t%0.10f\t9\n', H(1));
fprintf(fileID, '1\t50\tASO_MPC_H2\t%0.10f\t9\n', H(2));
fprintf(fileID, '1\t50\tASO_MPC_H3\t%0.10f\t9\n', H(3));
fprintf(fileID, '1\t50\tASO_MPC_H4\t%0.10f\t9\n', H(4));

%Hfinal
for i = 2 : 4
   for j = 2 : 4
       fprintf(fileID, '1\t50\tASO_MPC_HF%d%d\t%0.10f\t9\n', i, j, H_final(i,j));
   end
end

%lower bound
fprintf(fileID, '1\t50\tASO_MPC_LB1\t%0.10f\t9\n', lowerBound(1));
fprintf(fileID, '1\t50\tASO_MPC_LB2\t%0.10f\t9\n', lowerBound(2));

%upper bound
fprintf(fileID, '1\t50\tASO_MPC_UB1\t%0.10f\t9\n', upperBound(1));
fprintf(fileID, '1\t50\tASO_MPC_UB2\t%0.10f\t9\n', upperBound(2));

fclose(fileID);

%C definiction of matrix AExt
strAExt = 'const float mpc_AExt[3][3] = {';

for i = 1 : 3
   strAExt = [strAExt '{'];
   for j = 1 : 3
       strAExt = [strAExt num2str(AExt(i,j), '%.10f') 'f, '];  
   end
   strAExt(end-1) = '}';
   strAExt(end) = [];
   strAExt = [strAExt ',' 10];
end
strAExt(end-1) = [];

strAExt = [strAExt '};'];

%print matrices and info
display('---------- Info ----------');
display(['Model with estimated ' nameModelUsed '.']);
display(strAExt);
display(['LQR model sampling time: ' num2str(model.Dt) ...
        ' [sec] == ' num2str(round(model.Dt * 1e6)) ' [uSec].']);
display(['MPC model sampling time: ' num2str(modelMPC.Dt) ...
        ' [sec] == ' num2str(round(modelMPC.Dt * 1e6)) ' [uSec].']);
display('--------------------------');

s = {'ASO_LQR_K1' 'ASO_LQR_K2' 'ASO_LQR_K3'};
s1 = sprintf('%12s\t%12s\t%12s\t\n', s{:});
s2 = sprintf('%.10f\t%.10f\t%.10f\n', K_LQR);
display(horzcat(s1,s2));

%printmat(H', 'MPC H', '', 'ASO_MPC_H1 ASO_MPC_H2 ASO_MPC_H3 ASO_MPC_H4');
s = {'ASO_MPC_H1' 'ASO_MPC_H2' 'ASO_MPC_H3', 'ASO_MPC_H4'};
s1 = sprintf('%12s\t%12s\t%12s\t%12s\t\n', s{:});
s2 = sprintf('%.10f\t%.10f\t%.10f\t%.10f\n', H');
display(horzcat(s1,s2));


%printmat(lowerBound', 'MPC lower bound', '', 'ASO_MPC_LB1 ASO_MPC_LB2');
s = {'ASO_MPC_LB1' 'ASO_MPC_LB2'};
s1 = sprintf('%12s\t%12s\t\n', s{:});
s2 = sprintf('%.10f\t%.10f\n', lowerBound');
display(horzcat(s1,s2));

%printmat(upperBound', 'MPC uppper bound', '', 'ASO_MPC_UB1 ASO_MPC_UB2');
s = {'ASO_MPC_UB1' 'ASO_MPC_UB2'};
s1 = sprintf('%12s\t%12s\t\n', s{:});
s2 = sprintf('%.10f\t%.10f\n', upperBound');
display(horzcat(s1,s2));


% printmat(H_final(2, 2:end), 'MPC H final (2, 2:end)', '', ...
%         'ASO_MPC_HF22 ASO_MPC_HF23 ASO_MPC_HF24');
s = {'ASO_MPC_HF22' 'ASO_MPC_HF23' 'ASO_MPC_HF24'};
s1 = sprintf('%12s\t%12s\t%12s\t\n', s{:});
s2 = sprintf('%.10f\t%.10f\t%.10f\n', H_final(2, 2:end));
display(horzcat(s1,s2));


% printmat(H_final(3, 2:end), 'MPC H final (3, 2:end)', '', ...
%         'ASO_MPC_HF32 ASO_MPC_HF33 ASO_MPC_HF34');
s = {'ASO_MPC_HF32' 'ASO_MPC_HF33' 'ASO_MPC_HF34'};
s1 = sprintf('%12s\t%12s\t%12s\t\n', s{:});
s2 = sprintf('%.10f\t%.10f\t%.10f\n', H_final(3, 2:end));
display(horzcat(s1,s2));

% printmat(H_final(4, 2:end), 'MPC H final (4, 2:end)', '', ...
%         'ASO_MPC_HF42 ASO_MPC_HF43 ASO_MPC_HF44');
s = {'ASO_MPC_HF42' 'ASO_MPC_HF43' 'ASO_MPC_HF44'};
s1 = sprintf('%12s\t%12s\t%12s\t\n', s{:});
s2 = sprintf('%.10f\t%.10f\t%.10f\n', H_final(4, 2:end));
display(horzcat(s1,s2));    
