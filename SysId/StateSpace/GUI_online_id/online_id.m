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

% Last Modified by GUIDE v2.5 04-Mar-2015 09:31:31

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

%add tool folder to the path
addpath('tool');

%add MpcVsLqr folder to the path
addpath('MpcVsLqr');

%default values for weights, deltas and constraints
weights = [1, 3, 1, 35]; 
deltas = [12, 0.16]; 
constraints = [1, 3.6];
set(handles.t_weights, 'Data', weights);
set(handles.t_deltas, 'Data', deltas);
set(handles.t_constraints, 'Data', constraints);

%load default identified models
handles.idModels = tool_loadDefaultModels();
tool_updateModelList(handles);

% Update handles structure
guidata(hObject, handles);

%clean
clc;

% %debug
% path_name = ...
% 'C:\Users\IFA_sailing\Documents\FirmawareMarco\Matlab\as_matlab\SysId\StateSpace\GUI_online_id\exampleLogs';
% % path_name = ...
% %     'C:\Users\femer_000\Desktop\Tesi\mioCodice\Matlab\as_matlab\SysId\StateSpace\GUI_online_id';
% [~, logStr, logName] = tool_loadTxtLog('tack8_21_01_2015.txt', path_name);
% eval(['handles.logs.' logName ' = logStr;']);
% guidata(hObject, handles);
% 
% [~, logStr, logName] = tool_loadTxtLog('tack5B_21_01_2015.txt', path_name);
% eval(['handles.logs.' logName ' = logStr;']);
% guidata(hObject, handles);
% 
% [~, logStr, logName] = tool_loadTxtLog('cross6_18_02_2015.txt', path_name);
% eval(['handles.logs.' logName ' = logStr;']);
% guidata(hObject, handles);
% 
% tool_updateLogList(handles);

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
%open load window
[file_name, path_name] = uigetfile('*.txt', 'MultiSelect', 'on');

%if only a file has been selected, convert file name to call type
if(iscell(file_name) == 0)
    file_name = {file_name};
end

errorInHeader = 0;
%try to load each selected txt file
for i = 1 : length(file_name)
    [errorInHeader, logStr, logName] = tool_loadTxtLog(file_name{i}, path_name);
    
    %safety check and alert box
    if(errorInHeader == 1)
        msgbox('Please log yawspeed, yaw and rudder in QGC', 'Error','error');
        break;
    else
        %no errors, save data in handles
        eval(['handles.logs.' logName ' = logStr;']);
        %update log list
        tool_updateLogList(handles);
    end
end

%if no errors, save handles
if(errorInHeader == 0)
    guidata(hObject, handles);
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

%take sampling time to use in d2d command
resamplingTime = str2double(get(handles.e_newSampling, 'String'));

if(strcmp(selectedLog, 'log list') ~= 1)
    %call ginput to allow the user selecting starting and ending of Id data
    %[timeSelected, ~] = ginput(2);
    [timeSelected, ~] = ginputc(2, 'Color', 'r', 'LineWidth', 2);
    
    %identify model, if possibile
    eval(['logStr = handles.logs.' selectedLog ';']);    
    [retVal, model] = tool_idModel(selectedModel, timeSelected, logStr, resamplingTime);
    
    %add model to the list of identified models
    modelType = {'_black', '_grey'};
    eval(['handles.idModels.' ...
         selectedLog modelType{selectedModel} '_sampFact' num2str(resamplingTime) ...
         ' = model;']);
    guidata(hObject, handles);
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

%get selected id model
[error, modelSelected, indexModel] = tool_getSelectedIdModel(handles);

%check error
if(error == 0)
    %print the model selected in the text box and put the same model
    % in the MPC and LQR model
    tool_printIdModel(handles, modelSelected, indexModel);
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
[error, modelSelected, ~] = tool_getSelectedIdModel(handles);

%check if any error
if(error == 0)
    
    updateState = get(handles.c_updateState, 'Value');
    if(updateState == 1)
        %after how much time the model state has to be updated with the validation
        %state?
        stateUpdateTime = tool_getPredHorSec(handles, modelSelected);
    else
        %do not update the model state with the real one
        stateUpdateTime = -1;
    end
    
    %make sure selectedLog ~= log list
    if(strcmp(selectedLog, 'log list') ~= 1)
        eval(['validationLog = handles.logs.' selectedLog ';']);
        tool_validateModel(handles, modelSelected, validationLog, stateUpdateTime);
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



function e_newSampling_Callback(hObject, eventdata, handles)
% hObject    handle to e_newSampling (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
samplingFactor = str2double(get(hObject, 'String'));

%make sure samplingFactor is an integer number >= 1 
if(isnan(samplingFactor))
    msgbox('Please insert a number', 'Error','error');
    set(hObject, 'String', '1');
elseif(samplingFactor < 1)
    msgbox('Please insert a number >= 1', 'Error','error');
    set(hObject, 'String', '1');
else
    %force samplingFactor to be an integer
    set(hObject, 'String', num2str(round(samplingFactor)));
end



% --- Executes during object creation, after setting all properties.
function e_newSampling_CreateFcn(hObject, eventdata, handles)
% hObject    handle to e_newSampling (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in b_exportModel.
function b_exportModel_Callback(hObject, eventdata, handles)
% hObject    handle to b_exportModel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in b_exportData.
function b_exportData_Callback(hObject, eventdata, handles)

%see which models have been selected to design the LQR and MPC.
[error, lqrModel, nameLqr, mpcModel, nameMpc] = tool_getLqrMpcModels(handles);

%error check
if(error == 0) 
    %export to txt file new weights, deltas and constraints values
    weights = get(handles.t_weights, 'Data');
    deltas = get(handles.t_deltas, 'Data');
    constraints = get(handles.t_constraints, 'Data');
    
    %take the prediction horizon, in steps
    contents = cellstr(get(handles.p_predHorizon,'String')); 
    strPredHor = contents{get(handles.p_predHorizon,'Value')};
    predHor_steps = str2double(strPredHor);
    
    tool_exportData(lqrModel, nameLqr, mpcModel, nameMpc, ...
                    weights, deltas, constraints, predHor_steps);

end


% --- Executes on button press in b_MpcVsLqr.
function b_MpcVsLqr_Callback(hObject, eventdata, handles)
%see which models have been selected to design the LQR and MPC.
[errorDesignModel, lqrModel, ~, mpcModel, ~] = tool_getLqrMpcModels(handles);

%take the model to be considered as "real" in the simulation
[errorRealModel, realModel] = tool_getRealModel(handles);

%check errors
if(errorDesignModel == 0 && errorRealModel == 0)
    %take the prediction horizon, in steps
    contents = cellstr(get(handles.p_predHorizon,'String')); 
    strPredHor = contents{get(handles.p_predHorizon,'Value')};
    predHor_steps = str2double(strPredHor);
    
    %take weigths, deltas and constrints value
    weights = get(handles.t_weights, 'Data');
    deltas = get(handles.t_deltas, 'Data');
    constraints = get(handles.t_constraints, 'Data');
    
    %which type of tack should we make?
    typeTack = get(handles.p_typeTack, 'Value');
    
    %do we have to add noise in the simulation?
    addNoise = get(handles.c_addNoise, 'Value');
    %debug
%     assignin('base', 'realModel', realModel);
%     assignin('base', 'lqrModel', lqrModel);
%     assignin('base', 'mpcModel', mpcModel);
%     assignin('base', 'predHor_steps', predHor_steps);
%     assignin('base', 'weights', weights);
%     assignin('base', 'deltas', deltas);
%     assignin('base', 'constraints', constraints);
%     assignin('base', 'typeTack', typeTack);
    %simulate MPC and LQR response
    sim_MpcVsLqr(realModel, lqrModel, mpcModel, predHor_steps, ...
                 weights, deltas, constraints, typeTack, addNoise);
end


% --- Executes on selection change in p_predHorizon.
function p_predHorizon_Callback(hObject, eventdata, handles)
%take the selected identified model
[error, modelSelected, ~] = tool_getSelectedIdModel(handles);

%check error
if(error == 0)
    %based on the prediction horizon, in steps, compute the prediction horizon
    %in seconds
    predHor_s = tool_getPredHorSec(handles, modelSelected);

    %show it on screen
    set(handles.t_predHorizon, 'String', ...
        ['prediction horizon: ' num2str(predHor_s) ' [sec].']);
end


% --- Executes during object creation, after setting all properties.
function p_predHorizon_CreateFcn(hObject, eventdata, handles)
% hObject    handle to p_predHorizon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in c_updateState.
function c_updateState_Callback(hObject, eventdata, handles)
% hObject    handle to c_updateState (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of c_updateState


% --- Executes on selection change in p_realModel.
function p_realModel_Callback(hObject, eventdata, handles)
% hObject    handle to p_realModel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns p_realModel contents as cell array
%        contents{get(hObject,'Value')} returns selected item from p_realModel


% --- Executes during object creation, after setting all properties.
function p_realModel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to p_realModel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in p_typeTack.
function p_typeTack_Callback(hObject, eventdata, handles)
% hObject    handle to p_typeTack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns p_typeTack contents as cell array
%        contents{get(hObject,'Value')} returns selected item from p_typeTack


% --- Executes during object creation, after setting all properties.
function p_typeTack_CreateFcn(hObject, eventdata, handles)
% hObject    handle to p_typeTack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in p_lqrModel.
function p_lqrModel_Callback(hObject, eventdata, handles)
% hObject    handle to p_lqrModel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns p_lqrModel contents as cell array
%        contents{get(hObject,'Value')} returns selected item from p_lqrModel


% --- Executes during object creation, after setting all properties.
function p_lqrModel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to p_lqrModel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in p_mpcModel.
function p_mpcModel_Callback(hObject, eventdata, handles)
% hObject    handle to p_mpcModel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns p_mpcModel contents as cell array
%        contents{get(hObject,'Value')} returns selected item from p_mpcModel


% --- Executes during object creation, after setting all properties.
function p_mpcModel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to p_mpcModel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in c_addNoise.
function c_addNoise_Callback(hObject, eventdata, handles)
% hObject    handle to c_addNoise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of c_addNoise


% --- Executes on button press in b_plotPSD.
function b_plotPSD_Callback(hObject, eventdata, handles)
[error, logStr] = tool_getSelectedLog(handles);

if(error == 0)
    tool_plotPSD(logStr);
end



function e_cutOffFreq_Callback(hObject, eventdata, handles)
% hObject    handle to e_cutOffFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cutOffFreq = str2double(get(hObject, 'String'));

%make sure cutOffFreq is > 0
if(isnan(cutOffFreq))
    msgbox('Please insert a number', 'Error','error');
    set(hObject, 'String', '0.5');
elseif(cutOffFreq < 0)
    msgbox('Please insert a number > 0', 'Error','error');
    set(hObject, 'String', '0.5');
end


% --- Executes during object creation, after setting all properties.
function e_cutOffFreq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to e_cutOffFreq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function e_orderFilt_Callback(hObject, eventdata, handles)
% hObject    handle to e_orderFilt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
filterOrder = str2double(get(hObject, 'String'));

%make sure filterOrder is an integer number >= 1 
if(isnan(filterOrder))
    msgbox('Please insert a number', 'Error','error');
    set(hObject, 'String', '10');
elseif(filterOrder < 1)
    msgbox('Please insert a number >= 1', 'Error','error');
    set(hObject, 'String', '10');
else
    %force filterOrder to be an integer
    set(hObject, 'String', num2str(round(filterOrder)));
end


% --- Executes during object creation, after setting all properties.
function e_orderFilt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to e_orderFilt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in b_filt.
function b_filt_Callback(hObject, eventdata, handles)
%filtfilt yawRate and yaw
[error, filtFiltLog, ~] = tool_filtFiltLogStr(handles);

%check error
if(error == 0)
    %plot original log Vs filtfilt one
    tool_plotFiltFiltLog(handles, filtFiltLog);
end



% --- Executes on button press in b_addToLog.
function b_addToLog_Callback(hObject, eventdata, handles)
%filtfilt yawRate and yaw
[error, filtFiltLog, nameFiltFilt] = tool_filtFiltLogStr(handles);

%check error
if(error == 0)
    %add filt filt log to log list
    eval(['handles.logs.' nameFiltFilt ' = filtFiltLog;']);
    guidata(hObject, handles);
    %update log list
    tool_updateLogList(handles);
end


% --- Executes on button press in b_exportToWork.
function b_exportToWork_Callback(hObject, eventdata, handles)
[error, modelSelected, ~, modelName] = tool_getSelectedIdModel(handles);

%check error
if(error == 0)
    %save into workspace
    eval(['assignin(''base'', ''' modelName ''', modelSelected);']);
end
