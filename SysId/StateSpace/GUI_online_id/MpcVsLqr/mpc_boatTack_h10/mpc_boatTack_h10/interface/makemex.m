% mpc_boatTack_h10 : A fast customized optimization solver.
% 
% Copyright (C) 2013-2015 EMBOTECH GMBH [info@embotech.com]. All rights reserved.
% 
% 
% This software is intended for simulation and testing purposes only. 
% Use of this software for any commercial purpose is prohibited.
% 
% This program is distributed in the hope that it will be useful.
% EMBOTECH makes NO WARRANTIES with respect to the use of the software 
% without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
% PARTICULAR PURPOSE. 
% 
% EMBOTECH shall not have any liability for any damage arising from the use
% of the software.
% 
% This Agreement shall exclusively be governed by and interpreted in 
% accordance with the laws of Switzerland, excluding its principles
% of conflict of laws. The Courts of Zurich-City shall have exclusive 
% jurisdiction in case of any dispute.
% 

mex -c -O -DUSEMEXPRINTS ../src/mpc_boatTack_h10.c 
mex -c -O -DMEXARGMUENTCHECKS mpc_boatTack_h10_mex.c
if( ispc )
    mex mpc_boatTack_h10.obj mpc_boatTack_h10_mex.obj -output "mpc_boatTack_h10" 
    delete('*.obj');
elseif( ismac )
    mex mpc_boatTack_h10.o mpc_boatTack_h10_mex.o -output "mpc_boatTack_h10"
    delete('*.o');
else % we're on a linux system
    mex mpc_boatTack_h10.o mpc_boatTack_h10_mex.o -output "mpc_boatTack_h10" -lrt
    delete('*.o');
end
copyfile(['mpc_boatTack_h10.',mexext], ['../../mpc_boatTack_h10.',mexext], 'f');
copyfile( 'mpc_boatTack_h10.m', '../../mpc_boatTack_h10.m','f');
