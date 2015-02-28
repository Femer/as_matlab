function avgFs = tool_computeAvgSamplingF(time_ms)

%time in milliseconds
avgFs = 1 / mean(diff(time_ms) / 1e3);

end

