function buildMPC(predHoriz_steps)
%% FORCES multistage form

%assume variable ordering zi = [uHat_{i}; xHat_{i+1}] for i=1...N
%useful index with this ordering
indUHat = 1; %uHat = u_{k} - u_{k-1}
indW = 2;
indY = 3;
indU_k_minus_1 = 4;

%problem with 3 states and 1 input
nx = 3;
nu = 1;

stages = MultistageProblem(predHoriz_steps+1);

for i = 1 : predHoriz_steps+1
    
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
%         if(i < predHoriz_steps+1)
%             stages(i).eq.C = [zeros(nx, nu), AExt];
%         end
        if(i > 1)
            stages(i).eq.c = zeros(nx, 1);
        end
        %stages(i).eq.D = [BExt, -eye(nx)];
        
end

params(1) = newParam('minusAExt_times_x0', 1, 'eq.c'); % RHS of first eq. constr. is a parameter   
params(2) = newParam('Hessians', 1:predHoriz_steps, 'cost.H', 'diag');
params(3) = newParam('HessiansFinal', predHoriz_steps+1, 'cost.H');
params(4) = newParam('lowerBound', 1:predHoriz_steps+1, 'ineq.b.lb');
params(5) = newParam('upperBound', 1:predHoriz_steps+1, 'ineq.b.ub');
%equality constraints
params(6) = newParam('C', 1:predHoriz_steps, 'eq.C');
params(7) = newParam('D', 1:predHoriz_steps+1, 'eq.D');

%% Define outputs of the solver
outputs(1) = newOutput('u0', 1, 1:nu);

%% Solver settings
name = ['mpc_boatTack_h' num2str(predHoriz_steps)];
codeoptions = getOptions(name);
codeoptions.platform = 'ARM Cortex-M4 (FPU)';
codeoptions.accuracy.ineq = 1e-4;  % infinity norm of residual for inequalities
codeoptions.accuracy.eq = 1e-4;    % infinity norm of residual for equalities
codeoptions.accuracy.mu = 1e-4;    % absolute duality gap 
codeoptions.accuracy.rdgap = 1e-4; % relative duality gap := (pobj-dobj)/pobj
codeoptions.floattype = 'float';
codeoptions.printlevel = 0;
%test suggested by Alex L.
codeoptions.mu0 = 1;

%% Generate code
generateCode(stages, params, codeoptions, outputs);

end

