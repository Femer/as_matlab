function [xEst_k_k] = kfUpdate(model, K_k, measurements_k, u_k1, xEst_k1_k1)
A = model.A;
B = model.B;
C = model.C;

%predict state at step k with information up to k-1
xEst_k_k1 = A * xEst_k1_k1 + B * u_k1;

%update state prediction
xEst_k_k = xEst_k_k1 + K_k * (measurements_k - C * xEst_k_k1);
end

