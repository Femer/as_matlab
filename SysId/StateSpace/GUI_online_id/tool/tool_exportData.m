function  tool_exportData(model, nameModel, ...
                                   weights, deltas, constraints)


% Weights 
qYawRate = weights(1);
qYaw = weights(2);
rU = weights(3);
sChattering = weights(4);

%sampling time of the model, in uSec
samplingTime = round(model.Dt * 1e6);

%constraints
rudderMax = constraints(1);%
rudderVelocity = constraints(2);

%deltas
deltaYawRate = deltas(1); %in deg!
deltaYaw = deltas(2); %in deg!
deltaRudder = deltas(3);


%build cost matrices
Q = blkdiag(qYawRate, qYaw, rU);
R = sChattering;

%build extended state and extended model
AExt = [model.A,                      model.B;
        zeros(1, length(model.A)),    1];
    
BExt = [model.B;
        1];
    
%compute gain matrix for the LQR and the solution of the discrete riccati
%equation
[K_LQR, M, ~] = dlqr(AExt, BExt, Q, R); 

%MPC forces Hesian matrix
H = [sChattering; qYawRate; qYaw; rU];

%MPC forces Hesian matrix for final cost
H_final = blkdiag(0, M);

%convert rudderVel from cmd/s using the sample time of the model
rudderVelSim = rudderVelocity * model.Dt;
%MPC lower and upper bound
lowerBound = [-rudderVelSim; -rudderMax];
upperBound = [rudderVelSim; rudderMax];

%write txt file
%text string for the .txt file
fileID = fopen([nameModel '.txt'], 'w');
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

%delta values
fprintf(fileID, '1\t50\tASO_DLT_YR_D\t%0.10f\t9\n', deltaYawRate);
fprintf(fileID, '1\t50\tASO_DLT_Y_D\t%0.10f\t9\n', deltaYaw);
fprintf(fileID, '1\t50\tASO_DLT_RD_CM\t%0.10f\t9\n', deltaRudder);

%sampling time
fprintf(fileID, '1\t50\tASO_SPL_T_US\t%d\t9\n', samplingTime);
%TODO adjust king of param since ASO_SPL_T_US will bel an INT32

%TODO export model matrix A and B

fclose(fileID);

end

