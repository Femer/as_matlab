function [u, v, w, time] = tool_extractVelocites(tStart, tStop, est0Struct)

    est_index = est0Struct.time_est0 >= tStart & est0Struct.time_est0 <= tStop;
    
    time = est0Struct.time_est0(est_index);

    est0_vel_n = est0Struct.EST0_s4(est_index);
    est0_vel_e = est0Struct.EST0_s5(est_index);
    est0_vel_d = est0Struct.EST0_s6(est_index);
    
    q0 = est0Struct.EST0_s0(est_index);
    q1 = est0Struct.EST0_s1(est_index);
    q2 = est0Struct.EST0_s2(est_index);
    q3 = est0Struct.EST0_s3(est_index);
    u = zeros(length(q0), 1);
    v = zeros(length(q0), 1);
    w = zeros(length(q0), 1);
    
    for j = 1 : length(q0)
        rotation_matrix = tool_R_body_ned(q0(j),  q1(j), q2(j), q3(j));
        app = rotation_matrix * [est0_vel_n(j); est0_vel_e(j); est0_vel_d(j)];
        u(j) = app(1); 
        v(j) = app(2); 
        w(j) = app(3); 
    end


end

