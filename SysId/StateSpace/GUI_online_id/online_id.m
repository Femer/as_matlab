function varargout = online_id(varargin)
% ONLINE_ID MATLAB code for online_id.fig
%      ONLINE_ID, by itself, creates a new ONLINE_ID or raises the existing
%      singleton*.
%
%      H = ONLINE_ID returns the handle to a new ONLINE_ID or the handle to
%      the existing singleton*.
%
%      ONLINE_ID('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ONLINE_ID.M with the given input arguments.
%
%      ONLINE_ID('Property','Value',...) creates a new ONLINE_ID or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before online_id_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to online_id_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help online_id

% Last Modified by GUIDE v2.5 21-Feb-2015 15:25:00

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @online_id_OpeningFcn, ...
                   'gui_OutputFcn',  @online_id_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before online_id is made visible.
function online_id_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to online_id (see VARARGIN)

% Choose default command line output for online_id
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes online_id wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = online_id_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)



% --- Executes on button press in b_load_logs.
function b_load_logs_Callback(hObject, eventdata, handles)
[file_name, path_name] = uigetfile('*.txt');

%load('datiProva');

if ~isequal(file_name, 0) %if valid files has been selected
    
    
    full_path = fullfile(path_name, file_name);
    
    fileID = fopen(full_path);
    
    header = textscan(fileID, '%s %s %s %s', 1);
    celldisp(header)
    
    rows = textscan(fileID, '%f %f %f %f');
    celldisp(rows)
    
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
        neededField = {'TIME', 'yawspeed$', 'yaw$', 'Rud$'};
        for i = 1 : length(header)
           for j = 1 : length(neededField)
               if(~isempty(regexp(header{i}, neededField{j})))
                   %found a needed field in header{i}
                   if(j == 1)
                       timeIndex = j;
                   elseif(j == 2)
                       yawRateIndex = j;
                   elseif(j == 3)
                       yawIndex = j;
                   else
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
        end
    end

    %safety check and alert box
    if(errorInHeader == 1)
        msgbox('Log yawspeed, yaw and rudder in QGC', 'Error','error');
    else
        %no errors, save data in handles
        log.time = rows{timeIndex};
        log.yawSpeed = rows{yawRateIndex};
        log.yaw = rows{yawIndex};
        log.rudder = rows{rudderIndex};
        %call the log as the name of the loaded file
        pointIndex = regexp(file_name, '\.');
        logName = file_name(1 : pointIndex - 1);
        eval(['handles.logs.' logName ' = log;']);
        guidata(hObject, handles);
        
        %debug
        assignin('base', 'h', handles);
        assignin('base', 'file_name', file_name);
    end
    
end


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
