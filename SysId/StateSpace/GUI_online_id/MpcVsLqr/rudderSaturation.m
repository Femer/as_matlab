function newRudHat = rudderSaturation( 	rudHat, lastRudCmd,...
                                        rudderMax, rudderVelocity_cmd_s,...
                                        simDt)


%rudder command to the real system (not the extended one)
rud = rudHat + lastRudCmd;

%rudder velocity in simulation time
rudderVelSim = rudderVelocity_cmd_s * simDt;

%velocity constrain
if(rudHat >= rudderVelSim)
    %velocity constrain violated
    if(rudHat >= 0)
        rud = lastRudCmd + rudderVelSim;
    else
        rud = lastRudCmd - rudderVelSim;
    end
end

%saturation constrin
if(rud > rudderMax)
    rud = rudderMax;
elseif(rud < -rudderMax)
    rud = -rudderMax;
end

%compute new rudHat
newRudHat = rud - lastRudCmd;
   
end

