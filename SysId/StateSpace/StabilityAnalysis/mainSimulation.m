clc;
clear;
%close all;

%choose kp anc cp values
kp = 0.6366;%0.35 %P_tack kp = 0.73661977
cp = 0;%0.35;       %P_tack cp = 0.1  

%initial state x0 = [yawRate0; yaw0], in deg!
x0 = [0; -90];
x0 = deg2rad(x0);

%final time simulation
tF = 9;

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

%state evolution
x = zeros(2, N+1);
x(:, 1) = x0;
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
figure;

set(gcf,'name', ...
    ['kp = ' num2str(kp) '; cp = ' num2str(cp)], ...
    'numbertitle', 'off');

lW = 1.8;

subplot(311);
plot(time, rad2deg((x(1,:))), 'LineWidth', lW);
xlabel('Time [s]');
ylabel('w [deg/s]');
grid on;
xlim([time(1) time(end)]);

subplot(312);
plot(time, rad2deg((x(2,:))), 'LineWidth', lW);
xlabel('Time [s]');
ylabel('\psi [deg]');
grid on;
xlim([time(1) time(end)]);

subplot(313);
plot(time, rudder, 'LineWidth', lW, 'Color', 'm');
xlabel('Time [s]');
ylabel('rud [cmd]');
grid on;
xlim([time(1) time(end)]);

