function dataOut = simController(controller, typeController,...
                                 realModel, controllerModel, params, ...
                                 mpcParams)

%usefull index
yawRateIndex = 1;
yawIndex = 2;
%another useful index for the extended state
lastRudderIndex = 3;

%typecontroller possible values
indexMpcCtr = 1;
indexLqrCtr = 2;

%sample time of the controller w.r.t. realModel timing
controllerSampleTime = round(controllerModel.Dt / realModel.Dt);

%state evolution of the real system    
nx = size(realModel.A, 1);
xHatSim = zeros(nx, params.N);

%uHat = zeros(1, params.N);

%initial value of covariance prediction error, estimated step and first
%control input
P_k1_k1 = params.guessP1_1;
%P_k_k = P_k1_k1;
%xHatEst_k1_k1 = params.guessX1Hat;

%xHatEst = zeros(nx, params.N - 1);
%xHatEst(:, 1) = xHatEst_k1_k1;
%xEst_k_k = xHatEst_k1_k1;
%rudder value before tacking
u_k1 = params.xHatReal0(lastRudderIndex);

%time k = 1, the state situation BEFORE starting the tack

%real system
xHatSim(:, 1) = params.xHatReal0;

%controller data
time(1) = 1;
%realRudder for the real system, not for the extended one!
realRud(1) = xHatSim(lastRudderIndex, 1);

count = 2;

for k = 1 : params.N - 1
    
   %see if we have to compute a new optimal control action
   if(mod(k, controllerSampleTime + 1) == 0)
       
       %is this the very first time we compute  the optimal action?
       if(k == controllerSampleTime + 1)
           %pick a random value of the state for the KF near the state of
           %the "real" system
           xHatEst_k1_k1 = xHatSim(:, k-1) + [ sqrt(params.varYawRate) * randn();
                                             sqrt(params.varYaw) * randn();
                                             sqrt(params.varRudder) * randn()];
       end
       
       %save time
       time(count) = k;
       
       %we are starting now the step k, prediction phase
       [K_k, P_k_k] = kfPrediction(controllerModel, params.convarianceStr, P_k1_k1);
              
       %now we read the corrupted measurements
       meas_k = xHatSim(:, k) + params.measNoise(:, k);
       
       %update step to predict the real state
       xEst_k_k = kfUpdate(controllerModel, K_k, meas_k, u_k1, xHatEst_k1_k1);
       
       %Sinche the real rudder does not affect the state of the real boat,
       %here we use the real rudder command without noise on it
       xEst_k_k(3) = xHatSim(3, k);
       
       %save predicted step for later plots
       xHatEst(:, count) = xEst_k_k;
              
       %compute new optimal action
       if(typeController == indexMpcCtr)
           %MPC controller
           mpcParams.minusAExt_times_x0 = -controllerModel.A * xHatEst(:, count); 
           
           %MPC solve has to be used
           [solverout, exitflag, info] = controller(mpcParams);
                     
           if(exitflag == 1)
               uHat(k) = solverout.u0;
           else
               %display(info);
               display(['Some at iteration ' num2str(k) ' problem in solver, recovery rud cmd']);
               uHat(k) = u_k1;
           end
       else
           %LQR controller with negative feedback
           uHat(k) = - controller * xHatEst(:, count);
       end
              
       %constrain max rudder value and velocity
       %real previousRudder
       prevRealRud = realRud(count-1);

       %saturation and velocity constraints
      uHat(k) = rudderSaturation(uHat(k), xHatSim(3, k),...
                                 params.rudderMax, params.rudderVel_cmd_sec, ...
                                 controllerModel.Dt);
       %real rudder
       realRud(count) = uHat(k) + prevRealRud;
       
       count = count + 1;
       
       %save kalman filter variables at the end of step k
       P_k1_k1 = P_k_k;
       xHatEst_k1_k1 = xEst_k_k;
   else
       %do not compute a new optimal action, use the last computed one,
       %taht is, use uHat = 0
       uHat(k) = 0;
   end
   
   %update system dynamic
   xHatSim(:, k+1) = realModel.A * xHatSim(:, k) + realModel.B * uHat(k);
   
   %save last uHat given to the extended "real" model
   u_k1 = uHat(k);
end

%save output data
dataOut.time = time;
dataOut.xHatSim = xHatSim;
dataOut.xHatEst = xHatEst;
dataOut.typeController = typeController;

%from uHatcompute rudder sequence for the normal system (not the
%extended one)
dataOut.rudCmd = realRud;%xHatEst(lastRudderIndex, :);

end

