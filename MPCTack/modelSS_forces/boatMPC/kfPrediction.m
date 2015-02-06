function [K_k, P_k_k] = kfPrediction(model, convarianceStr, P_k1_k1)

A = model.A;
C = model.C;
R1 = convarianceStr.R1;
R2 = convarianceStr.R2;

%predict error covariance matrix at step k with information up to k-1
P_k_k1 = A * P_k1_k1 * A' + R1;

%compute Kalman gain at step k
K_k = P_k_k1 * C' * inv(C * P_k_k1 * C' + R2);

%compute error covariance at step k
P_k_k = (eye(size(C)) - K_k * C) * P_k_k1;

end

