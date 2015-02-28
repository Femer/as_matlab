function tool_plotPSD(logStr)

%compute mean sampling frequency
avgFs = tool_computeAvgSamplingF(logStr.time);

%plot yawRate nad yaw PSD
figure;

subplot(211);
[psdYR, f] = pwelch(logStr.yawRate, [], [], [], avgFs);
plot(f, psdYR, 'LineWidth', 1.7);
grid on;
xlabel('Frequency [Hz]');
ylabel('Magnitude [dB]');
title('YawRate');

subplot(212);
[psdY, f] = pwelch(logStr.yaw, [], [], [], avgFs);
plot(f, psdY, 'LineWidth', 1.7);
grid on;
xlabel('Frequency [Hz]');
ylabel('Magnitude [dB]');
title('YawRate');


end

