function [xHatSimLQR, xHatEstLQR, rudLQR] = simLQR( K_LQR, realModel, ...
                                                uncertModel, xHatSimLQR1,...
                                                guessP1_1, guessX1Hat,...
                                                rudderBeforeTack, N, ...
                                                measNoise, useRealState, ...
                                                rudderVel_cmd_sec, convarianceStr, ...
                                                rudderMax)

display('Computing LQR response');

%usefull index
yawRateIndex = 1;
yawIndex = 2;
%another useful index for the extended state
lastRudderIndex = 3;

%rudder velocity in simulation time
rudderVelSim = rudderVel_cmd_sec * realModel.Dt;

%state evolution of the real system    
nx = size(realModel.A, 1);
xHatSimLQR = zeros(nx, N);

%init
xHatSimLQR(:, 1) = xHatSimLQR1;
rudHatLQR = zeros(1, N);

%initial value of covariance prediction error, estimated step and first
%control input
P_k1_k1 = guessP1_1;
xHatEst_k1_k1 = guessX1Hat;

xHatEstLQR = zeros(nx, N-1);
xHatEstLQR(:, 1) = guessX1Hat;

%rudder value before tacking
u_k1 = rudderBeforeTack;

for k = 1 : N-1
   %we are starting now the step k, prediction phase
   [K_k, P_k_k] = kfPrediction(uncertModel, convarianceStr, P_k1_k1);
   
   %now we read the corrupted measurements
   meas_k = xHatSimLQR(:, k) + measNoise(:, k);
   
   %update step to predict the real state
   [xEst_k_k] = kfUpdate(uncertModel, K_k, meas_k, u_k1, xHatEst_k1_k1);
   
   %save predicted step for later plots
   if(k == 1)
        xHatEstLQR(:, 1) = guessX1Hat;
   else
        xHatEstLQR(:, k) = xEst_k_k;
   end
   
   if useRealState
       %compute LQR control input using real state
       feedbackState = xHatSimLQR(:, k);
   else
       %compute LQR control input using meas
       feedbackState = xEst_k_k;
   end
   rudHatLQR(k) = -K_LQR * feedbackState;

   %constrin max value and velocity
   rudHatLQR(k) = rudderSaturation(rudHatLQR(k), xHatSimLQR(lastRudderIndex, k),...
                                   rudderMax, rudderVel_cmd_sec, realModel.Dt);
      
   %update system dynamic
   xHatSimLQR(:, k+1) = realModel.A * xHatSimLQR(:, k) + realModel.B * rudHatLQR(k);
   
   %save kalman filter variables at the end of step k
   P_k1_k1 = P_k_k;
   u_k1 = rudHatLQR(k);
   xHatEst_k1_k1 = xEst_k_k;
end

%from uHatMPC compute rudder sequence for the normal system (not the
%extended one)
rudLQR = xHatSimLQR(lastRudderIndex, 1 : N-1);

   
%rudLQR = cumsum([rudderBeforeTack, rudHatLQR]);

end

