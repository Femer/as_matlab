clc;
clear;
close all;

%choose kp anc cp values
kp = 0.73661977;%0.35 %P_tack kp = 0.73661977
cp = 0.1;%0.35;       %P_tack cp = 0.1

%final time simulation
tF = 11;

%rudder constraints
rudMax = 1;
rudVel_cmd_s = 4;

%load numeric models
models = load('defaultIdentifiedModels');
%select which model you want to use
m = models.model1;

%steps
N = round(tF / m.Dt);

%rudder velocity in simulation time
rudVel = rudVel_cmd_s * m.Dt;

%state space sampling M_x0 = [x0_1, x0_2, etc]
yawRateSamples = deg2rad(-30 : 5 :30);
yawSamples = deg2rad([-110 : 10 : -70, 70 : 10 : 110]);

count = 1;
for yR = yawRateSamples
    for y = yawSamples
        M_x0(:,count) = [yR; y];
        count = count + 1;
    end
end

%for every initial conditions compute the system response
for k = 1 : size(M_x0,2)
    %state evolution
    x = zeros(2, N+1);
    x(:, 1) = M_x0(:,k);
    rudder = zeros(1,N+1);
    %rudder before tacking = 0
    rudder(1) = 0;
    
    for i = 2 : N+1
        rudder(i) = (kp / (1 + cp * abs(x(2,i-1)))) * x(2,i-1);
        %saturation
        if(abs(rudder(i)) > rudMax)
            rudder(i) = sign(rudder(i)) * rudMax;
        end
        %velocity limit
        if(abs(rudder(i) - rudder(i-1)) > rudVel)
            rudder(i) = rudder(i-1) + sign(rudder(i)) * rudVel;
        end
        %system evolution
        x(:, i) = m.A * x(:,i-1) + m.B * rudder(i-1);
    end
    
    %plots
    time = (0:N) .* m.Dt;
    if(k == 1)
        figure;
    else
        hold on;
    end
    
    lW = 1.0;
    
    plot(rad2deg(x(1,:)), rad2deg(x(2,:)), 'LineWidth', lW);
end

xlabel('X1');
ylabel('X2');
grid on;
%xlim([time(1) time(end)]);
