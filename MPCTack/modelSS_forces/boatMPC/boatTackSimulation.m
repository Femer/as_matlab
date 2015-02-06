clc;
clear;
%close all;

typeOfModel = 'little';%little (a,b) or capital (A,B)

if(strcmp(typeOfModel, 'little'))
    %load identified model with a and b
    load('linModelScalar');
    linModel = linModelScalar;
    display('Model with a and b.');
    %select wich  model you want to use
    nameModel = 'tack8';
else
    %load identified model with A and B
    load('linModelFull');
    linModel = linModelFull;
    display('Model with A and B.');
    %select wich  model you want to use
    nameModel = 'tack6';
end

%% tunable parameters

%take the linear model to use in MPC
eval(['model = linModel.' nameModel ';']);

% Weights 
qYawRate = 0.00001;
qYaw = 10;
rU = 0.001;
sChattering = 5;

%matrices
Q = blkdiag(qYawRate, qYaw, rU);
R = sChattering;


%rudder value before starting the tack
rudderBeforeTack = 0; %between -0.9 and 0.9

%take the sample time of the selected model, in seconds
meanTsSec = model.Dt;

%total simulation time, in seconds
tF = 5;
%total simulation steps to reach tF
N = round(tF / meanTsSec);
          
%initial conditions
tack = 'p2s'; %p2s s2p
absAlphaNew = 45 * pi / 180;

%gaussian noise on measurements
varYawRate = 2 * pi / 180;
varYaw = 10 * pi / 180;
varRudder = 0.000001;
       
measNoise = [sqrt(varYawRate) * randn(1, N);
             sqrt(varYaw) * randn(1, N);
             sqrt(varRudder) * randn(1, N)];
         
%covariance matrices for Kalman filter
convarianceStr.R1 = 0 * eye(3); %perfect model
convarianceStr.R2 = blkdiag(varYawRate, varYaw, varRudder); %noisy measurements

%constraints
rudderMax = 0.9;%0.9;
rudderVelocity = 1.8 / 0.5; %command/sec: rudder can go full right to full lest in 0.5 sec
absDeltaYaw = 10 * pi / 180;

%convert rudderVelocity from command/sec to command/simulationStep

%every simulation step lasts meanTsSec seconds.
rudderVelSim = rudderVelocity * meanTsSec;


% extended model xHat = [yawRate_k, yaw_k, rudder_{k-1}],
% uHat = [rudder_{k} - rudder{k-1}];
% we use the extended model to achieve the same cost function in the
% LQR and in the MPC. Since the LQR can have a cost function of the form
% x' * Q * x + u' * R * u, that tries to bring the system to the origin, we
% start with a yaw angle = -yawRef and we bring the system to the origin.
AExt = [model.A,                      model.B;
        zeros(1, length(model.A)),    1];
    
BExt = [model.B;
        1];
    
CExt = blkdiag(1, 1, 1);

%use extended state space model
model.A = AExt;
model.B = BExt;
model.C = CExt;

%use extended state space model
[nx, nu] = size(BExt);
        
%compute yawReference based on type of tack
if(strcmp(tack, 'p2s'))
    alpha0 = 45 .* pi / 180;
    yawRef = -absAlphaNew - alpha0;
else
    alpha0 = -45 .* pi / 180;
    yawRef = +absAlphaNew - alpha0;
end

%remeber that you have extended the state of the mdoel
xHatRef = [ 0;
            0;
            0];
%guess on the initial state of the KF
guessX1Hat = [  2 * pi / 180;
                -yawRef + (5 * pi / 180);
                rudderBeforeTack];
guessP1_1 = blkdiag(0.2 * eye(2), 0);

%usefull index
yawRateIndex = 1;
yawIndex = 2;
%another useful index for the extended state
lastRudderIndex = 3;

%Weight matrix for final cost
[~, M] = dlqr(AExt, BExt, Q, R);

% simulate
display('Computing MPC response');

%start with yawRate = 0 and yaw = -yawRef 
xHatSimMPC1 = [    0;
                -yawRef;
                rudderBeforeTack];
            
xHatSimMPC = zeros(nx, N);

%init
xHatSimMPC(:, 1) = xHatSimMPC1;
rudHatMPC = zeros(1, N);

%initial value of covariance prediction error
P_k1_k1 = guessP1_1;
xHatEst_k1_k1 = guessX1Hat;

xHatEstMPC = zeros(nx, N-1);
% xHatEstMPC(:, 1) = guessX1Hat;

%rudder value before tacking
u_k1 = rudderBeforeTack;

%assume variable ordering zi = [ui; xi+1] for i=1...N
problem.z1 = zeros(nx+nu,1);%TODO: CHECK THIS!
for k = 1 : N-1
    
   %we are starting now the step k, prediction phase
   [K_k, P_k_k] = kfPrediction(model, convarianceStr, P_k1_k1);
   
   %now we read the corrupted measurements
   meas_k = xHatSimMPC(:, k) + measNoise(:, k);
   
   %update step to predict the real state at step k with info up to k
   [xEst_k_k] = kfUpdate(model, K_k, meas_k, u_k1, xHatEst_k1_k1);

   %save predicted step for later plots
   if(k == 1)
       xHatEstMPC(:, 1) = guessX1Hat;
   else
       xHatEstMPC(:, k) = xEst_k_k;
   end
   
   %compute MPC control input using meas 
   problem.minusAExt_times_x0 = -AExt * xHatEstMPC(:, k);
   problem.Hessians = [sChattering; qYawRate; qYaw; rU]; %H
   problem.HessiansFinal = blkdiag(zeros(nu), M); %H at final step
   problem.lowerBound = [-rudderVelSim; -rudderMax];
   problem.upperBound = [rudderVelSim; rudderMax];
   
   [solverout, exitflag, info] = mpc_boatTack(problem);
   
   if(exitflag == 1)
       rudHatMPC(:, k) = solverout.u0;
   else
       display(info);
       error('Some problem in solver');
   end
   
   %update system dynamic
   xHatSimMPC(:, k+1) = AExt * xHatSimMPC(:, k) + BExt * rudHatMPC(k);
   
   %save kalman filter variables at the end of step k
   P_k1_k1 = P_k_k;
   u_k1 = rudHatMPC(k);
   xHatEst_k1_k1 = xEst_k_k;
end

% translate system response
xHatSimMPC(yawIndex, :) = xHatSimMPC(yawIndex, :) + yawRef;
yawMPC = xHatSimMPC(yawIndex, :);
xHatEstMPC(yawIndex, :) = xHatEstMPC(yawIndex, :) + yawRef;

%from uHatMPC compute rudder sequence for the normal system (not the
%extended one)
rudMPC = cumsum([rudderBeforeTack, rudHatMPC]);

%% plot
% ---- plot comparison ----

%compute simulation time, starting from time 0
time = (0:N) .* meanTsSec; 

%compute yaw overshoot limit
if(strcmp(tack, 'p2s'))
    limitYaw = yawRef - absDeltaYaw;
else
    limitYaw = yawRef + absDeltaYaw;
end

figure;

set(gcf,'name', ...
    ['MPC (forces), type of model: ' typeOfModel], ...
    'numbertitle', 'off');

lW0 = 1.3;
lW1 = 1.9;

leg = {{'\psi real', 'limit', '\psi by KF'}, ...
       'rud'};

yawReal = yawMPC;

%rudder for the real system, NOT for the extended one    
rudder = rudMPC;
    
%state estimated by the time vayring Kalman filter
yawEstVector = xHatEstMPC(2,:);


subplot(2, 1, 1);
title('state');

plot(time(2:end), yawReal .* 180 / pi, ...
    'LineWidth', lW1, 'Color', [88 25 225] ./ 255);
hold on;
plot([time(2) time(end)], [limitYaw limitYaw] .* 180 / pi, 'r-.', 'LineWidth', lW0);
plot(time(2:end-1), yawEstVector .* 180 / pi, 'c-.', ...
    'LineWidth', lW1, 'Color', [245 86 1] ./ 255);
grid on;
ylabel('[deg]');
xlabel('Time [sec]');
xlim([time(1) time(end)]);
legend(leg{1});

%add yawRef as tick in Y axis
tmp = get(gca, 'YTick');
tmp = [tmp yawRef * 180 / pi];
tmp = sort(tmp);
set(gca, 'YTick', tmp);


subplot(2, 1, 2);
title('input');

plot(time, rudder, 'm', 'LineWidth', lW1);
hold on;
plot([time(1) time(end)], [rudderMax rudderMax], 'r-.', 'LineWidth', lW0);
plot([time(1) time(end)], [-rudderMax -rudderMax], 'r-.', 'LineWidth', lW0);
grid on;
ylabel('Rudder [cmd]');
xlabel('Time [sec]');
xlim([time(1) time(end)]);
legend(leg{2}, 'Location', 'southwest');


