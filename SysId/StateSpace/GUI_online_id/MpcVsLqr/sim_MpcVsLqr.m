function sim_MpcVsLqr(realModel, lqrModel, mpcModel, predHor_steps, ...
                      weights, deltas, constraints, typeTack, addNoise)

%mpc subfolders
p = genpath('MpcVsLqr/mpc_boatTack_h10');
addpath(p);

p = genpath('MpcVsLqr/mpc_boatTack_h15');
addpath(p);

p = genpath('MpcVsLqr/mpc_boatTack_h20');
addpath(p);

p = genpath('MpcVsLqr/mpc_boatTack_h25');
addpath(p);

p = genpath('MpcVsLqr/mpc_boatTack_h30');
addpath(p);
% Weights 
qYawRate = weights(1);
qYaw = weights(2);
rU = weights(3);
sChattering = weights(4);

%constraints
rudderMax = constraints(1);
rudderVelocity = constraints(2); %cmd / s

%take sampling time of the real model and of the lqr and mpc models
realDt_s = realModel.Dt;
lqrDt_s = lqrModel.Dt;
mpcDt_s = mpcModel.Dt;

%lqrDt_s and mpcDt_s must be >= idDt_s
if(lqrDt_s < realDt_s || mpcDt_s < realDt_s)
   msgbox({'The system assumed as the real one', ...
           'must have a Dt smaller than the ones used', ...
           'to compute the MPC and the LQR'}, 'Error','error'); 
   return;
end

%build cost matrices
Q = blkdiag(qYawRate, qYaw, rU);
R = sChattering;
%build hessina matrix for MPC
H = [sChattering; qYawRate; qYaw; rU];

%build extended state and extended model
realExtMod = tool_extendModel(realModel);
lqrExtMod = tool_extendModel(lqrModel);
mpcExtMod = tool_extendModel(mpcModel);

[nx, nu] = size(realExtMod.B);
    
%compute gain matrix for the LQR 
[K_LQR, ~, ~] = dlqr(lqrExtMod.A, lqrExtMod.B, Q, R);

%compute final cost of the MPC
M = dare(mpcExtMod.A, mpcExtMod.B, Q, R);
H_final = blkdiag(0, M);

%before starting the two simulations, load usefull paramters
initParam = loadInitParams(realModel, typeTack);

%simulaions

%typecontroller possible values
indexMpcCtr = 1;
indexLqrCtr = 2;

%LQR, no mpcParams
lqrData = simController(K_LQR, indexLqrCtr,...
                        realExtMod, lqrExtMod, initParam, ...
                        [], addNoise);

%MPC 

%which MPC horizon do we have?
if(predHor_steps == 10)
    mpcHandler = @mpc_boatTack_h10;
elseif(predHor_steps == 15)
    mpcHandler = @mpc_boatTack_h15;
elseif(predHor_steps == 20)
    mpcHandler = @mpc_boatTack_h20;
elseif(predHor_steps == 25)
    mpcHandler = @mpc_boatTack_h25;
elseif(predHor_steps == 30)
    mpcHandler = @mpc_boatTack_h30;
else
    error('No valid MPC solver for this prediction horizon!');
end

%collect mpcParams
mpcParams.Hessians = H;
mpcParams.HessiansFinal = H_final;
%veolicty from cmd/s to cmd/samplingMpcModel
velMpcSamplingTime = rudderVelocity * mpcModel.Dt;
mpcParams.lowerBound = [-velMpcSamplingTime; -rudderMax];
mpcParams.upperBound = [velMpcSamplingTime; rudderMax];
mpcParams.C = [zeros(nx, nu), mpcExtMod.A];
mpcParams.D = [mpcExtMod.B, -eye(nx)];


mpcData = simController(mpcHandler, indexMpcCtr,...
                       realExtMod, mpcExtMod, initParam, ...
                       mpcParams, addNoise);
          
%plot comparison
plotComparison(lqrData, mpcData, initParam, deltas, ...
               mpcParams, predHor_steps, weights, addNoise, mpcModel.Dt);

end

