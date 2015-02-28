function tool_plotPSD(logStr)

%compute mean sampling frequency
avgFs = tool_computeAvgSamplingF(logStr.time);

%plot yawRate nad yaw PSD
figure;

h1 = subplot(211);
[psdYR, f] = pwelch(logStr.yawRate, [], [], [], avgFs);
plot(f, psdYR, 'LineWidth', 1.7);
grid on;
xlabel('Frequency [Hz]');
ylabel('Magnitude [dB]');
title('YawRate');

h2 = subplot(212);
[psdY, f] = pwelch(logStr.yaw, [], [], [], avgFs);
plot(f, psdY, 'LineWidth', 1.7);
grid on;
xlabel('Frequency [Hz]');
ylabel('Magnitude [dB]');
title('YawRate');

%link axes when zoom or move
linkaxes([h1, h2], 'x');
%xlim(h1, [time_sec(1) time_sec(end)]);
end

