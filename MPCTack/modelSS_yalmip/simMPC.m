function [xHatSimMPC, xHatEstMPC, rudMPC, timeRunMPC] = simMPC(mpcController, ...
                                                    realModel, uncertModel, ...
                                                    xHatSimMPC1, guessP1_1, guessX1Hat,...
                                                    rudderBeforeTack, factorSampleTime, ...
                                                    N, measNoise, useRealState, ...
                                                    convarianceStr, rudderMax)
% compute the 'real' system evolution using the normally sampled model called realModel.
% run the MPC every 'factorSampleTime' 
display('Computing MPC response');

%usefull index
yawRateIndex = 1;
yawIndex = 2;
%another useful index for the extended state
lastRudderIndex = 3;

nx = size(uncertModel.A, 1);
%state evolution of the real system     
xHatSimMPC = zeros(nx, N);

%init
xHatSimMPC(:, 1) = xHatSimMPC1;
rudHatMPC = [];

%initial value of covariance prediction error
P_k1_k1 = guessP1_1;
xHatEst_k1_k1 = guessX1Hat;

xHatEstMPC = zeros(nx, N-1);

%rudder value before tacking
u_k1 = rudderBeforeTack;

%simulation steps when MPC run
timeRunMPC = [];

%index to acces rudHatMPC and timeRunMPC
indexRunMPC = 1;

%input to the extended system
uHat = 0;

for k = 1 : N-1
    
   %we are starting now the step k, prediction phase using uncertModel
   [K_k, P_k_k] = kfPrediction(uncertModel, convarianceStr, P_k1_k1);
   
   %now we read the corrupted measurements adding noise to the 'real' state
   meas_k = xHatSimMPC(:, k) + measNoise(:, k);
   
   %update step to predict the real state at step k with info up to k using uncertModel
   [xEst_k_k] = kfUpdate(uncertModel, K_k, meas_k, u_k1, xHatEst_k1_k1);

   %save predicted step for later plots
   if(k == 1)
       xHatEstMPC(:, 1) = guessX1Hat;
   else
       xHatEstMPC(:, k) = xEst_k_k;
   end
   
   %compute MPC control every timeComputeMPCSim steps
   if(mod(k, factorSampleTime) == 0)
       
       if useRealState
           %compute new optimal control using real state
           rudHatMPC(indexRunMPC) = mpcController{xHatSimMPC(:, k)};
       else
           %compute new optimal control using meas
           rudHatMPC(indexRunMPC) = mpcController{xEst_k_k};
       end
              
       %save simulation step when a new optimal control was comptuted
       timeRunMPC(1, indexRunMPC) = k;
       %update uHat
       uHat = rudHatMPC(indexRunMPC);
       indexRunMPC = indexRunMPC + 1;
   else
       %no new MPC was computed, do not change the rudder input, in the
       %extended state model this means that uHat is 0
       uHat = 0;
   end
   
   %update system dynamic using the last rudder input computed
   xHatSimMPC(:, k+1) = realModel.A * xHatSimMPC(:, k) + realModel.B * uHat;
   
   %save kalman filter variables at the end of step k
   P_k1_k1 = P_k_k;
   u_k1 = uHat;
   xHatEst_k1_k1 = xEst_k_k;
end


%from uHatMPC compute rudder sequence for the normal system (not the
%extended one)
rudMPC = cumsum([rudderBeforeTack, rudHatMPC]);
%at time 0, rudder was equal to rudderBeforeTack
timeRunMPC = [0, timeRunMPC];
end

