function dataOut = simController(jollyParam, typeController,...
                                 realModel, controllerModel, params, ...
                                 mpcParams)
%debug
typeController

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

%init "real" system state
xHatSim(:, 1) = params.xHatReal0;
uHat = zeros(1, params.N);

%initial value of covariance prediction error, estimated step and first
%control input
P_k1_k1 = params.guessP1_1;
P_k_k = P_k1_k1;
xHatEst_k1_k1 = params.guessX1Hat;

xHatEst = zeros(nx, params.N - 1);
xHatEst(:, 1) = xHatEst_k1_k1;
xEst_k_k = xHatEst_k1_k1;
%rudder value before tacking
u_k1 = params.xHatReal0(lastRudderIndex);

count = 1;

for k = 1 : params.N - 1
    
   %see if we have to compute a new optimal control action
   if(mod(k, controllerSampleTime) == 0)
       %save time
       time(count) = k;
       %we are starting now the step k, prediction phase
       [K_k, P_k_k] = kfPrediction(controllerModel, params.convarianceStr, P_k1_k1);
              
       %now we read the corrupted measurements
       meas_k = xHatSim(:, k) + params.measNoise(:, k);
       
       %update step to predict the real state
       xEst_k_k = kfUpdate(controllerModel, K_k, meas_k, u_k1, xHatEst_k1_k1);
       
       %save predicted step for later plots
       if(k == 1)
           xHatEst(:, 1) = params.guessX1Hat;
       else
           xHatEst(:, count) = xEst_k_k;
       end
              
       %compute new optimal action
       if(typeController == indexMpcCtr)
           %MPC controller
           mpcParams.minusAExt_times_x0 = -controllerModel.A * xHatEst(:, count); 
           
           %based on predHor_steps ( == jollyParams) value see which MPC solve has to be used
            if(jollyParam == 10)
                [solverout, exitflag, info] = mpc_boatTack_h10(mpcParams);
            else
                error('error!');
            end

           if(exitflag == 1)
               uHat(:, k) = solverout.u0;
           else
               display(info);
               display('Some problem in solver, recovery rud cmd');
               uHat(k) = u_k1;
           end
       else
           %LQR controller with negative feedback
           uHat(k) = - jollyParam * xHatEst(:, count);
       end
              
       %constrain max rudder value and velocity
       uHat(k) = rudderSaturation(uHat(k), xHatSim(lastRudderIndex, k),...
                                  params.rudderMax, params.rudderVel_cmd_sec, ...
                                  realModel.Dt);
       
       count = count + 1;
   else
       %do not compute a new optimal action, use the last computed one
       uHat(k) = u_k1;
   end
   
   %update system dynamic
   xHatSim(:, k+1) = realModel.A * xHatSim(:, k) + realModel.B * uHat(k);
   
   %save kalman filter variables at the end of step k
   P_k1_k1 = P_k_k;
   u_k1 = uHat(k);
   xHatEst_k1_k1 = xEst_k_k;
end

%save output data
dataOut.time = time;
dataOut.xHatSim = xHatSim;
dataOut.xHatEst = xHatEst;
dataOut.typeController = typeController;

%from uHatcompute rudder sequence for the normal system (not the
%extended one)
dataOut.rudCmd = xHatSim(lastRudderIndex, 1 : params.N - 1);

end

