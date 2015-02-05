function [tack] = tool_selectData(start, stop, data)



tack.start = start;
tack.stop = stop;

att_index = data.att.time_att >= tack.start & data.att.time_att <= tack.stop;
bgud_index = data.bgud.time_bgud >= tack.start & data.bgud.time_bgud <= tack.stop;
rc_index = data.rc.time_rc >= tack.start & data.rc.time_rc <= tack.stop;
gps_index = data.gps.time_gps >= tack.start & data.gps.time_gps <= tack.stop;
wsai_index = data.wsai.time_wsai >= tack.start & data.wsai.time_wsai <= tack.stop;
est_index = data.est.time_est0 >= tack.start & data.est.time_est0 <= tack.stop;
gpos_index = data.gpos.time_gpos >= tack.start & data.gpos.time_gpos <= tack.stop;
%stupid work around, take last value for type of tack
qgc2_index = data.qgc2.time_qgc2 <= tack.start;
%REMEMBER that data.qgc2.time_qgc2(1) is Nan, so to have the correct index
%you have to add 1!!!!
index_last_past_changed = sum(qgc2_index) + 1;
valTack = data.qgc2.QGC2_TypTck(index_last_past_changed);
TypTck = valTack * ones(length(data.bgud.time_bgud(bgud_index)), 1);

tack.bgud_time = data.bgud.time_bgud(bgud_index);
tack.tack = data.bgud.BGUD_ShldTck(bgud_index);
tack.alpha = data.bgud.BGUD_Alpha(bgud_index);
tack.alpha_star = data.bgud.BGUD_AlphaStar(bgud_index);
tack.rudder = data.bgud.BGUD_Rudder(bgud_index);
tack.sails = data.bgud.BGUD_Sail(bgud_index);
tack.twd_avg = data.bgud.BGUD_TwdAvg(bgud_index);
tack.x = data.bgud.BGUD_XRace(bgud_index);
tack.y = data.bgud.BGUD_YRace(bgud_index);
tack.next_grid = data.bgud.BGUD_NextGrid(bgud_index);


tack.rc_time = data.rc.time_rc(rc_index);
tack.mode = data.rc.RC_Ch4(rc_index);

tack.att_time = data.att.time_att(att_index);
tack.roll = data.att.ATT_Roll(att_index);
tack.yaw = data.att.ATT_Yaw(att_index);
tack.yawRate = data.att.ATT_YawRate(att_index);

tack.cog = data.gps.GPS_Cog(gps_index);
tack.gps_time = data.gps.time_gps(gps_index);

tack.qgc2_time =  data.bgud.time_bgud(bgud_index);
tack.type_tack = TypTck;

tack.wsai_time = data.wsai.time_wsai(wsai_index);
tack.twd_raw = data.wsai.WSAI_AngleTrue(wsai_index);


tack.est_time = data.est.time_est0(est_index);
q0 = data.est.EST0_s0(est_index);
q1 = data.est.EST0_s1(est_index);
q2 = data.est.EST0_s2(est_index);
q3 = data.est.EST0_s3(est_index);

tack.velU = zeros(length(q0), 1);
tack.velV = zeros(length(q0), 1);
interp_vel_n = spline(data.gpos.time_gpos(gpos_index), ...
                      data.gpos.GPOS_VelN(gpos_index), ...
                      data.est.time_est0(est_index));
interp_vel_e = spline(data.gpos.time_gpos(gpos_index), ...
                      data.gpos.GPOS_VelE(gpos_index), ...
                      data.est.time_est0(est_index));
interp_vel_d = spline(data.gpos.time_gpos(gpos_index), ...
                      data.gpos.GPOS_VelD(gpos_index), ...
                      data.est.time_est0(est_index));
for j = 1 : length(q0)
        rotation_matrix = tool_R_body_ned(q0(j),  q1(j), q2(j), q3(j));
        app = rotation_matrix * [interp_vel_n(j); interp_vel_e(j); interp_vel_d(j)];
        tack.velU(j) = app(1); 
        tack.velV(j) = app(2); 
end

%compute alpha_yaw by interpolating
yaw_interp = spline(tack.att_time, tack.yaw, tack.bgud_time);
tack.yaw_interp = yaw_interp;
tack.alpha_yaw = tack.twd_avg - yaw_interp;

yawRate_interp = spline(tack.att_time, tack.yawRate, tack.bgud_time);
tack.yawRate_interp = yawRate_interp;

end

