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

% Last Modified by GUIDE v2.5 21-Feb-2015 17:08:36

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
%[file_name, path_name] = uigetfile('*.txt');


%if ~isequal(file_name, 0) %if valid files has been selected
if 1
file_name = 'test.txt';    
    
%     full_path = fullfile(path_name, file_name);
%     %read file
%     fileID = fopen(full_path); 
    fileID = fopen('test.txt'); 
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
        neededField = {'TIME', 'yawspeed$', 'yaw$', 'Rud$'};
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
        end
    end

    %safety check and alert box
    if(errorInHeader == 1)
        msgbox('Please log yawspeed, yaw and rudder in QGC', 'Error','error');
    else
        %no errors, save data in handles
        log.time = rows{timeIndex};
        log.yawRate = rows{yawRateIndex};
        log.yaw = rows{yawIndex};
        log.rudder = rows{rudderIndex};
        %call the log as the name of the loaded file
        pointIndex = regexp(file_name, '\.');
        logName = file_name(1 : pointIndex - 1);
        eval(['handles.logs.' logName ' = log;']);
        guidata(hObject, handles);
        %update log list
        tool_updateLogList(handles);
        %debug
        assignin('base', 'h', handles);
    end
    
end


% --- Executes on selection change in p_logList.
function p_logList_Callback(hObject, eventdata, handles)

contents = cellstr(get(hObject,'String')); %returns p_logList contents as cell array
selectedLog = contents{get(hObject,'Value')}; %returns selected item from p_logList

if(strcmp(selectedLog, 'log list') ~= 1)
   %plot selected log
   
   eval(['logStr = handles.logs.' selectedLog ';']);
   tool_plotLog(handles, logStr);
end


% --- Executes during object creation, after setting all properties.
function p_logList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to p_logList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in b_identify.
function b_identify_Callback(hObject, eventdata, handles)
%see which model has been selected
selectedModel = get(handles.p_typeModel,'Value'); 

%take the log selected in p_logList, it must be ~= first string
contents = cellstr(get(handles.p_logList,'String')); %returns p_logList contents as cell array
selectedLog = contents{get(handles.p_logList,'Value')}; %returns selected item from p_logList

if(strcmp(selectedLog, 'log list') ~= 1)
    %call ginput to allow the user selecting starting and ending of Id data
    [timeSelected, ~] = ginput(2);
    
    %identify model, if possibile
    eval(['logStr = handles.logs.' selectedLog ';']);    
    [retVal, model] = tool_idModel(selectedModel, timeSelected, logStr);
    
    %add model to the list of identified model
    modelType = {'_black', '_grey'};
    eval(['handles.idModels.' selectedLog modelType{selectedModel} ' = model;']);
    guidata(hObject, handles);
    %debug
    assignin('base', 'h', handles);
    %update idModel list
    tool_updateModelList(handles);
else
    %error, no valid log selected
    msgbox('Please select a valid log', 'Error','error');
end



% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in p_typeModel.
function p_typeModel_Callback(hObject, eventdata, handles)
% hObject    handle to p_typeModel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns p_typeModel contents as cell array
%        contents{get(hObject,'Value')} returns selected item from p_typeModel


% --- Executes during object creation, after setting all properties.
function p_typeModel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to p_typeModel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in p_idModels.
function p_idModels_Callback(hObject, eventdata, handles)

contents = cellstr(get(hObject,'String')); 
nameModel = contents{get(hObject,'Value')};

%if nameModel ~= identified models
if(strcmp(nameModel, 'identified models') ~= 1)
    %take the model
    eval(['modelSelected = handles.idModels.' nameModel ';']);
    %print the model selected in the text box
    tool_printIdModel(handles, modelSelected);
end


% --- Executes during object creation, after setting all properties.
function p_idModels_CreateFcn(hObject, eventdata, handles)
% hObject    handle to p_idModels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in b_validate.
function b_validate_Callback(hObject, eventdata, handles)
%take the selected log and use it as validation data
contents = cellstr(get(handles.p_logList,'String')); %returns p_logList contents as cell array
selectedLog = contents{get(handles.p_logList,'Value')}; %returns selected item from p_logList
%see which model has been selected
contents = cellstr(get(handles.p_idModels,'String')); 
nameModel = contents{get(handles.p_idModels,'Value')};

%take resaple factor
resampleFactor = str2double(get(handles.e_resampleFactor, 'String'));

%make sure selectedLog ~= log list AND nameModel ~= identified model
if(strcmp(selectedLog, 'log list') ~= 1)
    if(strcmp(nameModel, 'identified model') ~= 1)
        
        eval(['validationLog = handles.logs.' selectedLog ';']);   
        eval(['modelSelected = handles.idModels.' nameModel ';']);
        tool_validateModel(handles, modelSelected, validationLog, resampleFactor);
    else
        msgbox('Please select an identified model', 'Error','error');
    end
else
    msgbox('Please select a log to use as validation data', 'Error','error');
end




function e_resampleFactor_Callback(hObject, eventdata, handles)

sampleFactor = str2double(get(hObject, 'String'));

%make sure sampleFactor is a number
if(isnan(sampleFactor))
    msgbox('Please insert a number', 'Error','error');
    set(hObject, 'String', '-1');
end


% --- Executes during object creation, after setting all properties.
function e_resampleFactor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to e_resampleFactor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
