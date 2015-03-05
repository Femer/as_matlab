function printCDebugcode(mpcModel, mpcData, deltas, ...
    mpcParams, predHor_steps)

fileID = fopen('debugC.txt', 'w');

%print on screen yawRateSimMPC
yawRateMPC = mpcData.xHatEst(1,2:end);

fprintf(fileID, 'float yawRateMPC[] = {%ff', yawRateMPC(1));
for yR = yawRateMPC(2:end)
    fprintf(fileID, ',%ff', yR);
end
fprintf(fileID, '};\n');


%print on screen yawSimMPC and yawLength
yawMPC = mpcData.xHatEst(2,2:end);

fprintf(fileID, 'float yawMPC[] = {%ff', yawMPC(1));
for y = yawMPC(2:end)
    fprintf(fileID, ',%ff', y);
end
fprintf(fileID, '};\n');
fprintf(fileID, 'int yawLength = %d;\n', length(yawMPC));

fprintf(fileID, 'static int index_data_mpc = 0;\n');

%print h
h = mpcParams.Hessians;
fprintf(fileID, 'float h[] = {%ff,%ff,%ff,%ff};\n', h(1), h(2), h(3), h(4));


%print HFinal matrix
hf = mpcParams.HessiansFinal(2:end,2:end);
fprintf(fileID, 'float h_final[][3] = {{%ff,%ff,%ff}, \n {%ff,%ff,%ff}, \n {%ff,%ff,%ff}};\n',...
                                      hf(1,1), hf(1,2), hf(1,3), ...
                                      hf(2,1), hf(2,2), hf(2,3), ...
                                      hf(3,1), hf(3,2), hf(3,3));
%print lB and uB
low = mpcParams.lowerBound;
up = mpcParams.upperBound;

fprintf(fileID, 'float lb[] = {%ff,%ff};\n', low(1), low(2));
fprintf(fileID, 'float ub[] = {%ff,%ff};\n', up(1), up(2));

%print real A and B
A = mpcModel.A;

fprintf(fileID, 'float A[][2] = {{%ff,%ff}, {%ff,%ff}};\n', ...
                                A(1,1), A(1,2), A(2,1), A(2,2));
                            
B = mpcModel.B;

fprintf(fileID, 'float B[] = {%ff,%ff};\n', B(1), B(2));

%print sampling time in uSec
fprintf(fileID, 'int32_t mpc_sampling_time_us = %d;\n', round(mpcModel.Dt * 1e6));

%print prediction horizon in steps
fprintf(fileID, 'int32_t pred_horz_steps = %d;\n', predHor_steps);

%print delta values and 2 other params
fprintf(fileID, 'float delta[] = {%ff, %ff};\n', deg2rad(deltas(1)), deg2rad(deltas(2)));
fprintf(fileID, 'float min_time_s = %ff;\n', 2);
fprintf(fileID, 'float safety_time_stop_tack_s = %ff;\n', 320);


%print real Rud
correctRud = mpcData.rudCmd(2:end);

fprintf(fileID, 'float correctRud[%d] = {%ff', length(correctRud), correctRud(1));

for r = correctRud(2:end)
    fprintf(fileID, ',%ff', r);
end
fprintf(fileID, '};\n');


fclose(fileID);
end

