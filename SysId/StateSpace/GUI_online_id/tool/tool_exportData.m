function  tool_exportData(lqrModel, nameLqr, mpcModel, nameMpc, ...
                          weights, deltas, constraints, predHor_steps)


% Weights 
qYawRate = weights(1);
qYaw = weights(2);
rU = weights(3);
sChattering = weights(4);

%sampling time of the model used for the LQR, in uSec
lqrSamplingTime = round(lqrModel.Dt * 1e6);
mpcSamplingTime = round(mpcModel.Dt * 1e6);

%constraints
rudderMax = constraints(1);% cmd
rudderVelocity = constraints(2); %cmd/s

%deltas
deltaYaw = deltas(1); %in deg!
deltaRudder = deltas(2);


%build cost matrices
Q = blkdiag(qYawRate, qYaw, rU);
R = sChattering;

%build extended state and extended model
lqrExtMod = tool_extendModel(lqrModel);
mpcExtMod = tool_extendModel(mpcModel);
    
%compute gain matrix for the LQR 
[K_LQR, ~, ~] = dlqr(lqrExtMod.A, lqrExtMod.B, Q, R); 

%MPC forces Hesian matrix
H = [sChattering; qYawRate; qYaw; rU];

%MPC forces Hesian matrix for final cost using dare solution
M = dare(mpcExtMod.A, mpcExtMod.B, Q, R);
H_final = blkdiag(0, M);

%convert rudderVel from cmd/s using the sample time of the model
rudderVelSim = rudderVelocity * mpcModel.Dt;
%MPC lower and upper bound
lowerBound = [-rudderVelSim; -rudderMax];
upperBound = [rudderVelSim; rudderMax];

%write txt file
%text string for the .txt file
nameFile = ['lqr_' nameLqr '_mpc_' nameMpc];
fileID = fopen([nameFile '.txt'], 'w');
%LQR gain
fprintf(fileID, '1\t50\tASO_LQR_K1\t%0.10f\t9\n', K_LQR(1));
fprintf(fileID, '1\t50\tASO_LQR_K2\t%0.10f\t9\n', K_LQR(2));
fprintf(fileID, '1\t50\tASO_LQR_K3\t%0.10f\t9\n', K_LQR(3));

%H vector
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

%delta values
fprintf(fileID, '1\t50\tASO_DLT_Y_D\t%0.10f\t9\n', deltaYaw);
fprintf(fileID, '1\t50\tASO_DLT_RD_CM\t%0.10f\t9\n', deltaRudder);

%sampling time
fprintf(fileID, '1\t50\tASO_SPL_LQR_US\t%d\t6\n', lqrSamplingTime);
fprintf(fileID, '1\t50\tASO_SPL_MPC_US\t%d\t6\n', mpcSamplingTime);

%export model matrix A and B, for the MPC
fprintf(fileID, '1\t50\tASO_MPC_A11\t%0.10f\t9\n', mpcModel.A(1,1));
fprintf(fileID, '1\t50\tASO_MPC_A12\t%0.10f\t9\n', mpcModel.A(1,2));
fprintf(fileID, '1\t50\tASO_MPC_A21\t%0.10f\t9\n', mpcModel.A(2,1));
fprintf(fileID, '1\t50\tASO_MPC_A22\t%0.10f\t9\n', mpcModel.A(2,2));

fprintf(fileID, '1\t50\tASO_MPC_B1\t%0.10f\t9\n', mpcModel.B(1));
fprintf(fileID, '1\t50\tASO_MPC_B2\t%0.10f\t9\n', mpcModel.B(2));

%steps for the prediction horizon
fprintf(fileID, '1\t50\tASO_PRED_HOR\t%d\t6\n', predHor_steps);

fclose(fileID);

end

