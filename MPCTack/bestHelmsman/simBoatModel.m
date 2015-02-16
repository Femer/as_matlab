function [w, y, rud] = simBoatModel(model, x0, param, tf_sec)

simSteps = round(tf_sec / model.Dt);

x = zeros(2, simSteps);
rud = zeros(simSteps - 1, 1);

%init x0
x(1, 1) = x0(1);
x(2, 1) = x0(2);

%sim model
for k = 1 : simSteps - 1
    %now tack_p2s_rudder wants the yaw angle as first argument
    rud(k) = tack_p2s_rudder(x(2,k), param);
    
    x(:, k + 1) = model.A * x(:, k) + model.B  * rud(k);

end

%save yawRate and yaw angle
w = x(1,:);
y = x(2,:);

