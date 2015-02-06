clc;
clear;
close all;

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


% tunable parameters
%take the linear model to use in MPC
eval(['model = linModel.' nameModel ';']);

%prediction horizon, in number of simulation steps
predHor = 10;
    
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

%use extended state space model
[nx, nu] = size(BExt);

%take the sample time of the selected model, in seconds
meanTsSec = model.Dt;
display(['Horizon MPC: ' num2str(predHor * meanTsSec) ' [sec].']);

% %constraints
% rudderMax = 0.9;%0.9;
% rudderVelocity = 1.8 / 0.5; %command/sec: rudder can go full right to full lest in 0.5 sec
% absDeltaYaw = 10 * pi / 180;

%convert rudderVelocity from command/sec to command/simulationStep

%every simulation step lasts meanTsSec seconds.
%rudderVelSim = rudderVelocity * meanTsSec;

%% FORCES multistage form
%assume variable ordering zi = [uHat_{i}; xHat_{i+1}] for i=1...N
%useful index with this ordering
indUHat = 1; %uHat = u_{k} - u_{k-1}
indW = 2;
indY = 3;
indU_k_minus_1 = 4;

stages = MultistageProblem(predHor+1);

for i = 1 : predHor+1
    
        % dimension
        stages(i).dims.n = nx + nu; % number of stage variables
        stages(i).dims.r = nx;    % number of equality constraints        
        stages(i).dims.l = 2; % number of lower bounds
        stages(i).dims.u = 2; % number of upper bounds
        
        stages(i).cost.f = zeros(nx+nu,1); % linear cost terms
        
        % lower bounds on rudder saturation (using indU_k_minus_1) and
        % bounds on rudder velocity (using indUHat)
        stages(i).ineq.b.lbidx = [indUHat; indU_k_minus_1]; % lower bound acts on these indices
        
        % upper bounds on rudder saturation (using indU_k_minus_1) and
        % bounds on rudder velocity (using indUHat)
        stages(i).ineq.b.ubidx = [indUHat; indU_k_minus_1]; % upper bound acts on these indices
        
        % equality constraints
        if(i < predHor+1)
            stages(i).eq.C = [zeros(nx, nu), AExt];
        end
        if(i > 1)
            stages(i).eq.c = zeros(nx, 1);
        end
        stages(i).eq.D = [BExt, -eye(nx)];
        
end

params(1) = newParam('minusAExt_times_x0', 1, 'eq.c'); % RHS of first eq. constr. is a parameter   
params(2) = newParam('Hessians', 1:predHor, 'cost.H', 'diag');
params(3) = newParam('HessiansFinal', predHor+1, 'cost.H');
params(4) = newParam('lowerBound', 1:predHor+1, 'ineq.b.lb');
params(5) = newParam('upperBound', 1:predHor+1, 'ineq.b.ub');

%% Define outputs of the solver
outputs(1) = newOutput('u0', 1, 1:nu);

%% Solver settings
codeoptions = getOptions('mpc_boatTack');

%% Generate code
generateCode(stages, params, codeoptions, outputs);
