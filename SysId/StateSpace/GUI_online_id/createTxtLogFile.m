nameFile = 'data20_10_02_2015.txt';

tack = stepTacks.tack20;

y = tack.dataYaw.OutputData;
w = tack.dataYawRate.OutputData;
rud = tack.dataYaw.InputData;
t = tack.time_sysId;

t = (t - t(1)) ./ 1e3;

fileID = fopen(nameFile,'w');

fprintf(fileID, 'TIMESTAMPms	M1ATTITUDEyawspeed	M1ATTITUDEyaw	Rud\n');
for i = 1 : length(t)
   fprintf(fileID, '%i\t%f\t%f\t%f\n', [round(t(i)) w(i) y(i) rud(i)]); 
end
fclose(fileID);
