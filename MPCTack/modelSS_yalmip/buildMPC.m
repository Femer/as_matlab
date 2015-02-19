function mpcController = buildMPC(model, predHor, Q, R, typeTack, rudderMax, rudderVel_cmd_sec)

addpath(genpath('yalmip/'));

yalmip('clear');

display('Building MPC');

%compute rudder velocity based on the sample time of the model.
%rudderVel_cmd_sec MUST be expressed in cmd/sec
rudderVelMPC = rudderVel_cmd_sec * model.Dt;

%usefull index
yawRateIndex = 1;
yawIndex = 2;
%another useful index for the extended state
lastRudderIndex = 3;

%model matrix
A = model.A;
B = model.B;

%yalmip options
options = sdpsettings('solver', 'mosek', 'verbose', 1);

% Number of states and inputs
[nx, nu] = size(B); 

uHatMPC = sdpvar(repmat(nu,1,predHor), repmat(1,1,predHor));
xHatMPC = sdpvar(repmat(nx,1,predHor+1), repmat(1,1,predHor+1));

constraints = [];
objective = 0;

for k = 1 : predHor
    %remember that the extended model has the form:
    % xMPC{k} = [yawRate_k, yaw_k, rudder_{k-1}],
    % uMPC{k} = [rudder_{k} - rudder{k-1}];
    
    %update cost function
    objective = objective + ...
                norm(Q * xHatMPC{k}, 2) + ...
                norm(R * uHatMPC{k}, 2);
               
    %add system dynamic to constraints
    constraints = [constraints, xHatMPC{k+1} == A * xHatMPC{k} + B * uHatMPC{k}];
    
    %limit input action to be within feasible set
    %rudder to real system = uHatMPC{k} + xHatMPC{k}(lastRudderIndex)
    constraints = [constraints, ...
                  -rudderMax <= uHatMPC{k} + xHatMPC{k}(lastRudderIndex) <= rudderMax];
    
    %limit input velocity
    constraints = [constraints, abs(uHatMPC{k}) <= rudderVelMPC];
 
%     %limit maximum yaw overshoot based on which haul you want to have at
%     %the end
%     if(strcmp(typeTack, 'p2s'))
%         constraints = [constraints, xHatMPC{k}(yawIndex)  >= -absDeltaYaw];
%     else
%         constraints = [constraints, xHatMPC{k}(yawIndex) <= absDeltaYaw];
%     end

 
end

%compute Riccati solution and use it as final cost
[~, M, ~] = dlqr(A, B, Q, R);

%add final cost 
objective = objective + norm(M * xHatMPC{predHor + 1}, 2);

parameters_in = xHatMPC{1};
solutions_out = uHatMPC{1};

mpcController = optimizer(constraints, objective, options, parameters_in, solutions_out);

end

