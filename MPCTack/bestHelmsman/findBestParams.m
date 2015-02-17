clc;
clear;
close all;

%load samped trajectories
load('sampledTraj');

%best params minimize cost funtion J = sum((x - xRef)' * Q * (x - xRef)) + 
% + r * rud ^ 2.
%where x = [w; y]

%xRef
xRef = [ 0; 
         deg2rad(-45)];

%weights
qW = 0.01; % 0.01 short peak    | 0.5 very short peak
qY = 10; % 10                   |  10
r = 0.5;  % 5                   |  0.1

Q = blkdiag(qW, qY);

%% J function for Helmsman0 and Helmsman1 tack
for helm = [0, 1]
    display(['Computing J' num2str(helm)]);
    eval(['h = h' num2str(helm) ';']);
    J = zeros(length(h.w), 1);
    for i = 1 : length(h.w)

        x = [h.w{i};
             h.y{i}];
        rud = h.rud{i};
        %compute J
        J(i) = 0;
        for k = 1 : length(x)
            J(i) = J(i) + (x(:,k) - xRef)' * Q * (x(:,k) - xRef);
        end
        
        for k = 1 : length(rud)
            J(i) = J(i) + r * rud(k)^2;
        end

    end
    
    eval(['J' num2str(helm) ' = J;']);
end

%% find best params
indexMinJ0 = find(J0 == min(J0));
indexMinJ1 = find(J1 == min(J1));

indexMinJ0 = indexMinJ0(1);
indexMinJ10 = indexMinJ1(1);

display('---------- Best Helmsman 0 ----------');
param0 = h0.param{indexMinJ0};
display(['[X1, X2] star :' num2str(rad2deg(param0.xVector)) '   [rad].']);
display(['Y1star :' num2str(param0.yValue) ' [cmd].']);

display('---------- Best Helmsman 1 ----------');
param1 = h1.param{indexMinJ1};
display(['[X1, X2] star :' num2str(rad2deg(param1.xVector)) '   [rad].']);
display(['Y1star :' num2str(param1.yValue) ' [cmd].']);

figure;
lW = 1.9;
%Remeber: alpha_dumas = - yaw_sensor

time_state_h0 = (0 : length(h0.y{indexMinJ0}) - 1) .* h0.Dt;
time_rud_h0 = (0 : length(h0.rud{indexMinJ0}) - 1) .* h0.Dt;

time_state_h1 = (0 : length(h1.y{indexMinJ1}) - 1) .* h1.Dt;
time_rud_h1 = (0 : length(h1.rud{indexMinJ1}) - 1) .* h1.Dt;

subplot(221);
plot(time_state_h0, -rad2deg(h0.y{indexMinJ0}), 'LineWidth', lW);
hold on;
plot([time_state_h0(1) time_state_h0(end)], ...
    -rad2deg([xRef(2) xRef(2)]), ...
    'r--', 'LineWidth', lW);

legend('\alpha H0', '\alpha star');
xlabel('Time [sec]');
ylabel('deg');
grid on;

subplot(222);
plot(time_state_h1, -rad2deg(h1.y{indexMinJ1}), 'LineWidth', lW);
hold on;
plot([time_state_h0(1) time_state_h0(end)], ...
    -rad2deg([xRef(2) xRef(2)]), ...
    'r--', 'LineWidth', lW);


legend('\alpha H1', '\alpha star');
xlabel('Time [sec]');
ylabel('deg');
grid on;

subplot(223);
plot(time_rud_h0,h0.rud{indexMinJ0}, 'LineWidth', lW);
legend('rud H0');
xlabel('Time [sec]');
ylabel('cmd');
grid on;

subplot(224);
plot(time_rud_h1,h1.rud{indexMinJ1}, 'LineWidth', lW);
legend('rud H1');
xlabel('Time [sec]');
ylabel('cmd');
grid on;



