function [errorInHeader, logStr, logName] = tool_loadTxtLog(file_name, path_name)

full_path = fullfile(path_name, file_name);

if ~isequal(file_name, 0) %if valid files has been selected
    %read file
    fileID = fopen(full_path);
    header = textscan(fileID, '%s %s %s %s', 1);
    rows = textscan(fileID, '%f %f %f %f');
    fclose(fileID);
    
    %make sure in the txt file loaded there are timeStamp, yawRate, yaw and
    %rudder command.
    errorInHeader = 0;
    
    %index of each needed field in the header struct
    timeIndex = [];
    yawRateIndex = [];
    yawIndex = [];
    rudderIndex = [];
    
    if(length(header) ~= 4)
        errorInHeader = 1;
    else
        %use reg exp
        neededField = {'TIME', 'yawspeed$', 'yaw$', 'asrud$'};
        for i = 1 : length(header)
            for j = 1 : length(neededField)
                resRegExp = regexp(header{i}, neededField{j});
                if(~isempty(resRegExp{1}))
                    %found a needed field in header{i}
                    if(j == 1)
                        timeIndex = j;
                    elseif(j == 2)
                        yawRateIndex = j;
                    elseif(j == 3)
                        yawIndex = j;
                    elseif(j == 4)
                        rudderIndex = j;
                    end
                end
            end
        end
        %check if every needed field was found
        if(isempty(timeIndex) || ...
                isempty(yawRateIndex) || ...
                isempty(yawIndex) || ...
                isempty(rudderIndex))
            
            errorInHeader = 1;
        else
            %everything ok
            errorInHeader = 0;
            %remove rows with nan value, if any
            rows = tool_removeNan(rows);
            logStr.time = rows{timeIndex};
            logStr.yawRate = rows{yawRateIndex};
            %UNWRAP yaw angle!
            logStr.yaw = unwrap(rows{yawIndex});
            logStr.rudder = rows{rudderIndex};
            %call the log as the name of the loaded file
            pointIndex = regexp(file_name, '\.');
            logName = file_name(1 : pointIndex - 1);
        end
        
    end
else
    errorInHeader = 1;%no file selected
end

end

