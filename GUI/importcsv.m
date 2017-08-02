function varargout = importcsv(varargin)
% IMPORTCSV MATLAB code for importcsv.fig
%      IMPORTCSV, by itself, creates a new IMPORTCSV or raises the existing
%      singleton*.
%
%      H = IMPORTCSV returns the handle to a new IMPORTCSV or the handle to
%      the existing singleton*.
%
%      IMPORTCSV('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IMPORTCSV.M with the given input arguments.
%
%      IMPORTCSV('Property','Value',...) creates a new IMPORTCSV or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before importcsv_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to importcsv_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help importcsv

% Last Modified by GUIDE v2.5 28-Jul-2017 12:25:50

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @importcsv_OpeningFcn, ...
                   'gui_OutputFcn',  @importcsv_OutputFcn, ...
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


% --- Executes just before importcsv is made visible.
function importcsv_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to importcsv (see VARARGIN)

% Choose default command line output for importcsv
h = findobj('Tag','figureMain');
% if exists (not empty)
if ~isempty(h)
    % get handles and other user-defined data associated
    figmain = guidata(h);
    handles.Qsyn = figmain.Qsyn;
    handles.dateformat = figmain.dateformat;
    handles.na = figmain.na;
    handles.axes = figmain.axes;
end

handles.output = hObject;
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes importcsv wait for user response (see UIRESUME)
% uiwait(handles.figureImportCSV);


% --- Outputs from this function are returned to the command line.
function varargout = importcsv_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function editFile_Callback(hObject, eventdata, handles)
handles.file = get(hObject,'String');
guidata(hObject, handles);

function editFile_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function pushbuttonBrowse_Callback(hObject, eventdata, handles)
[FileName,PathName] = uigetfile('*.csv');
filepath = [PathName FileName];
set(handles.editFile,'String',filepath);
handles.file = get(handles.editFile,'String');
guidata(hObject, handles);

function pushbuttonOK_Callback(hObject, eventdata, handles)
handles.Qsyn.importts(handles.file, handles.dateformat, handles.na);
handles.Qsyn.perennial();
handles.Qsyn.checkobs();
handles.Qsyn.gaps();
axes(handles.axes);
plot(handles.axes,handles.Qsyn.tt_obs.Date, handles.Qsyn.tt_obs.Q,'blue','visible','on');
title('Daily Streamflow');
grid;
xlabel('Date','FontWeight','bold');
ylabel('Q [m^3/s]','FontWeight','bold');
guidata(hObject, handles);
delete(handles.figureImportCSV);
msgbox('Loaded .csv successfully');

function pushbuttonCancel_Callback(hObject, eventdata, handles)
delete(handles.figureImportCSV);
