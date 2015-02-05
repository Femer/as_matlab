%take data when PI controllore was in charge to control the course of the
%boat. Use function selectTackData only to select data based on time
%% init
clc;
clear;

data.att = load('test1Att');
data.rc = load('test1RC');
data.bgud = load('test1BGUD');
data.gps = load('test1GPS');
data.qgc2 = load('test1QGC2');
data.wsai = load('test1WSAI');

data.est = load('test1EST0');
data.gpos = load('test1GPOS');

%add path for extra helper functions
addpath('../tools/');

% %test unifor mean sample time for every sequence, in uSec
% meanTs = 9855.6988;

%% course control 1

start = 3703813234; %usec
stop =  3842653235;%sec

track = tool_selectData(start, stop, data);
stepTacks.track1 = tool_createIdData(track);

%plotTackChar(track1);

%% course control 2

start = 4072144807 + 60 * 1e6; %usec
stop =  4310345930;%sec

track = tool_selectData(start, stop, data);
stepTacks.track2 = tool_createIdData(track);

%plotTackChar(track2);
%pay attention: here in some parts uVel is < 0

%% course control 3

start =  5184461609; %usec
stop =  5320354095;%sec

track = tool_selectData(start, stop, data);
stepTacks.track3 = tool_createIdData(track);

%plotTackChar(track3);
%pay attention: here in some parts uVel is < 0

%% course control 4

start =  4908548305; %usec
stop =  4995163518;%sec

track = tool_selectData(start, stop, data);
stepTacks.track4 = tool_createIdData(track);

%plotTackChar(track4);
%pay attention: here in some parts uVel is < 0

%% course control 5

start =  4723915637; %usec
stop =  4807945357;%sec

track = tool_selectData(start, stop, data);
stepTacks.track5 = tool_createIdData(track);

%plotTackChar(track5);
%pay attention: here in some parts uVel is < 0

%% course control 6

start =  4535888709; %usec
stop =  4592436372;%sec

track = tool_selectData(start, stop, data);
stepTacks.track6 = tool_createIdData(track);

