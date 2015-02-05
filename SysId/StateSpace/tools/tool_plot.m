function [] = tool_plot(tack, name)
bgud_rel_time = (tack.bgud_time - tack.bgud_time(1)) ./ 1e6;
rc_rel_time = (tack.rc_time - tack.rc_time(1)) ./ 1e6;
att_rel_time = (tack.att_time - tack.att_time(1)) ./ 1e6;
gps_rel_time = (tack.gps_time - tack.gps_time(1)) ./ 1e6;
wsai_rel_time = (tack.wsai_time - tack.wsai_time(1)) ./ 1e6;
qgc2_rel_time = (tack.qgc2_time - tack.qgc2_time(1)) ./ 1e6;
est_rel_time = (tack.est_time - tack.est_time(1)) ./ 1e6;


rowSub = 4;
colSub = 2;
figure;
set(gcf,'name', name,'numbertitle','off');

subplot(rowSub,colSub,1);
plot(bgud_rel_time, tack.tack, 'g');
hold on;
plot(rc_rel_time, tack.mode, 'r');
plot(qgc2_rel_time, tack.type_tack .* 0.8, 'k');


ylim([-1.1 1.1]);
legend('tack', 'mode', 'tackType');
xlabel('Time [s]');

subplot(rowSub,colSub,2);
plot(att_rel_time, tack.roll .* 180 / pi, 'r');
hold on;
plot(att_rel_time, tack.yaw .* 180 / pi, 'b');
plot(bgud_rel_time, tack.tack .* 100, 'g');
legend('roll', 'yaw', 'tack');
xlabel('Time [s]');
ylabel('Angle [deg]');


subplot(rowSub,colSub,3);
plot(bgud_rel_time, tack.alpha .* 180 / pi, 'Color', 'r');
hold on;
plot(bgud_rel_time, tack.alpha_star .* 180 / pi, 'b');
plot(bgud_rel_time, tack.alpha_yaw .* 180 / pi, 'c--');
legend('\alpha', '\alpha Star', '\alpha Yaw');
xlabel('Time [s]');
ylabel('Angle [deg]');

subplot(rowSub,colSub,4);
plot(bgud_rel_time, tack.rudder, 'r');
hold on;
plot(bgud_rel_time, tack.sails, 'b');
plot(bgud_rel_time, tack.tack, 'g');
plot(bgud_rel_time, tack.alpha_yaw .* 180 / pi / 100, 'c--');
legend('rud', 'sails', 'tack', '\alpha Yaw');
xlabel('Time [s]');
ylabel('Command');

subplot(rowSub,colSub,5);
plot(att_rel_time, tack.yaw .* 180 / pi, 'b');
hold on;
plot(gps_rel_time, tack.cog .* 180 / pi, 'r');
plot(bgud_rel_time, tack.tack .* 100, 'g');
legend('yaw', 'cog', 'tack');
xlabel('Time [s]');
ylabel('Angle[deg]');

subplot(rowSub,colSub,6);
plot(wsai_rel_time, tack.twd_raw .* 180 / pi, 'b');
hold on;
plot(bgud_rel_time, tack.twd_avg .* 180 / pi, 'r');
legend('twd', 'twd avg');
xlabel('Time [s]');
ylabel('Angle[deg]');

subplot(rowSub,colSub,7);
plot(est_rel_time, tack.velU);
hold on;
plot(bgud_rel_time, tack.tack, 'g');
legend('u', 'tack');
xlabel('Time [s]');
ylabel('Vel[m/s]');

subplot(rowSub,colSub,8);
plot(est_rel_time, tack.velV);
hold on;
plot(bgud_rel_time, tack.tack, 'g');
legend('v', 'tack');
xlabel('Time [s]');
ylabel('Vel[m/s]');

%subplot(rowSub,colSub,7);
% figure;
% a = 25;
% c = linspace(1,10,length(tack.x));
% scatter(tack.y, tack.x, a, c, 'filled');
% hold on;
% xlabel('Y [m]');
% ylabel('X[m]');
% axis equal;
% grid on;
% %take grid lines
% grids = unique(tack.next_grid);
% for i = 1 : length(grids)
%    val = grids(i);
%    plot([min(tack.y) max(tack.y)], [val val], 'c--', 'LineWidth', 2.5); 
% end

end

