function [realModelExt, uncertModelExt, uncertMPCModelExt] = extentendModels(realModel, ...
                                                                             uncertModel, ...
                                                                             uncertMPCModel)

% extended model xHat = [yawRate_k, yaw_k, rudder_{k-1}],
% uHat = [rudder_{k} - rudder{k-1}];
% we use the extended model to achieve the same cost function in the
% LQR and in the MPC.
olMod = {realModel, uncertModel, uncertMPCModel};

for i = 1 : length(olMod)
   extMod.A = [olMod{i}.A,                      olMod{i}.B;
               zeros(1, length(olMod{i}.A)),    1]; 
           
   extMod.B = [olMod{i}.B;
               1]; 
   extMod.C = eye(size(extMod.A));
   extMod.Dt = olMod{i}.Dt;
   
   if i == 1
       realModelExt = extMod;
   elseif i == 2
       uncertModelExt = extMod;
   else
       uncertMPCModelExt = extMod;
   end
   
end

end

