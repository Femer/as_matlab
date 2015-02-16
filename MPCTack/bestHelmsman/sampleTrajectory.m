clc;
clear;
close all;

%% load identified model

%simulation time, in seconds
tf_sec = 4.5;

%initila state in the simulation, tack port 2 starboard
x0 = [0;
      deg2rad(45)];

%estimated model with a and B or with A and B
typeOfModel = 'capital'; %little or capital

if(strcmp(typeOfModel, 'little'))
    %load identified model with scalar value for a and b
    load('linModelScalar-21-01');
    linModels = linModelScalar;
    display('Model estimated with a and b');
    %select the sequence from which you want to load the identified model
    nameSeqId = 'tack8';
    dateSysId = '21-01-2015';
else
    load('linModelFull-21-01');
    linModels = linModelFull;
    display('Model estimated with A and B');
    %select the sequence from which you want to load the identified model
    nameSeqId = 'tack6';
    dateSysId = '21-01-2015';
end

display(['Model estimated with data from ' nameSeqId '; estimated using ' typeOfModel ' model.']);
eval(['model = linModels.' nameSeqId ';']);

%% sample parameters x1, x2, y1 for helmsman0 and helmsman1
stepX1 = deg2rad(2.5);
stepX2 = deg2rad(2.5);

stepY1 = -0.025;
startY1 = 0.9;
stopY1 = 0;


startX1_h0 = deg2rad(-25);
stopX1_h0 = deg2rad(35);


stopX2_h0 = deg2rad(40);

startX1_h1 = deg2rad(-25);
startX2_h1 = deg2rad(-30);
stopX2_h1 = deg2rad(45);

%sample param for helmsman0
count = 1;
for x1 = startX1_h0 : stepX1 : stopX1_h0
    for x2 = (x1 + stepX2) : stepX2 : stopX2_h0
        for y1 = startY1 : stepY1 : stopY1
            param.xVector = [x1, x2];
            param.yValue = y1;
            param.typeHelmsman = 0;
            paramH0{count} = param;
            count = count + 1;
        end
    end
end

%sample param for helmsman1
count = 1;
for x2 = startX2_h1 : stepX2 : stopX2_h1
    for x1 = startX1_h1 : stepX1 : (x2 - stepX1)
        for y1 = startY1 : stepY1 : stopY1
            param.xVector = [x1, x2];
            param.yValue = y1;
            param.typeHelmsman = 1;
            paramH1{count} = param;
            count = count + 1;
        end
    end
end

%% sim tack helmsman and helmsman1 using sampled parameters
display('Sim Helm0');
for i = 1 : length(paramH0)
    [w, y, rud] = simBoatModel(model, x0, paramH0{i}, tf_sec);
    
    h0.w{i} = w;
    h0.y{i} = y;
    h0.rud{i} = rud;
    h0.param{i} = paramH0{i};
    
    if(mod(i, 100) == 0)
        display(['#it: ' num2str(i)]);
    end
end
h0.Dt = model.Dt;

display('Sim Helm1');
for i = 1 : length(paramH1)
    [w, y, rud] = simBoatModel(model, x0, paramH1{i}, tf_sec);
    
    h1.w{i} = w;
    h1.y{i} = y;
    h1.rud{i} = rud;
    h1.param{i} = paramH1{i};
    
    if(mod(i, 100) == 0)
        display(['#it: ' num2str(i)]);
    end
end
h1.Dt = model.Dt;

%save
save('sampledTraj', 'h0', 'h1');
