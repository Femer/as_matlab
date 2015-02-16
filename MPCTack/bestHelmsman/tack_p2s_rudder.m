function dr = tack_p2s_rudder(yaw, param)

%dummy init
dr = 0;

%alpha in dums convention = - yaw
alpha = -yaw;

%extract x1 and x2 param
x1 = param.xVector(1);
x2 = param.xVector(2);
% constant of the function

% max rudder angle when steering from right to left
y1 = param.yValue;

%helmens type 0
if(param.typeHelmsman == 0)
   
    if(alpha <= deg2rad(-30))
        dr = y1;
    elseif(alpha <= x1)
        dr = (-y1 / (x1 - deg2rad(-30))) * (alpha - x1);
    elseif(alpha <= x2)
        dr = (y1 / (x2 - x1)) * (alpha - x1);
    elseif(alpha <= deg2rad(40))
        dr = (-y1 / (deg2rad(40) - x2)) * (alpha - deg2rad(40));
    else
        dr = 0;
    end
else
    
    if(alpha <= x1)
        dr = y1;
    elseif(alpha <= x2)
        dr = (-y1 / (x2 - x1)) * (alpha - x2);
    else
        dr = 0;
    end
    
end
end

