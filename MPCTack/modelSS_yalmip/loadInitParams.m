%guess on the initial state of the KF
guessX1Hat = [  0 * pi / 180;
                -yawRef + (0 * pi / 180);
                rudderBeforeTack];
            
%guess on the variance matrix of the KF
guessP1_1 = blkdiag(1 * eye(2), 0);

%start of the realModel used to simulate the boat in the C.L. MPC
xHatSimMPC1 = [    3.103 * pi / 180;
                -yawRef - sign(yawRef) * deg2rad(7);
                rudderBeforeTack];
            
            
%start of the realModel used to simulate the boat in the C.L. LQR
xHatSimLQR1 = [    3.103 * pi / 180;
                -yawRef - sign(yawRef) * deg2rad(7);
                rudderBeforeTack];