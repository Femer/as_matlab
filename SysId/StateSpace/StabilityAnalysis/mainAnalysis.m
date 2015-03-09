clc;
clear;
close all;

%choose kp anc cp values
kp = 0.35;
cp = 0.35;
%x2 sampling steps (in rad)
x2_step = deg2rad(0.5);
%load numeric models
models = load('defaultIdentifiedModels');
%select which model you want to use
m = models.model1;

%stability if x2 == 0
display(['If x2 == 0: eig = [' num2str(eig(m.A)') ']']);

%stability if x2 > 0
x2 = deg2rad(0.1);

while(1)
    A_cl = m.A + [zeros(2,1), m.B * (kp / (1 + cp * x2))];
    
    if(x2 >= deg2rad(90) || any(abs(eig(A_cl)) > 1))
        break;
    end
    x2 = x2 + x2_step;
end

display(['C.L. system stable if x2 > 0 untill x2 < ' num2str(rad2deg(x2))]);

%stability if x2 < 0
x2 = deg2rad(-0.1);

while(1)
    A_cl = m.A + [zeros(2,1), m.B * (kp / (1 - cp * x2))];
    
    if(x2 <= deg2rad(-90) || any(abs(eig(A_cl)) > 1))
        break;
    end
    x2 = x2 - x2_step;
end

display(['C.L. system stable if x2 < 0 untill x2 < ' num2str(rad2deg(x2))]);
