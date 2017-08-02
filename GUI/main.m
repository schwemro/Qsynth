function varargout = main(varargin)
% MAIN MATLAB code for main.fig
%      MAIN, by itself, creates a new MAIN or raises the existing
%      singleton*.
%
%      H = MAIN returns the handle to a new MAIN or the handle to
%      the existing singleton*.
%
%      MAIN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MAIN.M with the given input arguments.
%
%      MAIN('Property','Value',...) creates a new MAIN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before main_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to main_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
%      Requires MY_XTICKLABELS.
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help main

% Last Modified by GUIDE v2.5 27-Jul-2017 15:23:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @main_OpeningFcn, ...
    'gui_OutputFcn',  @main_OutputFcn, ...
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


% --- Executes just before main is made visible.
function main_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to main (see VARARGIN)

% Choose default command line output for main
handles.output = hObject;
handles.Qsyn = Qsynth(); % intilialize class
handles.approach = 1; % 1 = 'AR' or 2 = 'ARMA'
handles.l = 10; % number of iterations
handles.dateformat = [];
handles.na = [];
handles.startdate = [];
handles.enddate = [];
cla(handles.axes,'reset') % fresh empty plot
cla(handles.axes1,'reset')
cla(handles.axes2,'reset')
cla(handles.axes3,'reset')
cla(handles.axes4,'reset')
cla(handles.axes5,'reset')
set(handles.axes,'visible','on');
set(handles.axes1,'visible','off');
set(handles.axes2,'visible','off');
set(handles.axes3,'visible','off');
set(handles.axes4,'visible','off');
set(handles.axes5,'visible','off');
guidata(hObject, handles); % Update handles structure

% UIWAIT makes main wait for user response (see UIRESUME)
% uiwait(handles.figureMain);


% --- Outputs from this function are returned to the command line.
function varargout = main_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function editDate_Callback(hObject, eventdata, handles)
handles.dateformat = get(hObject,'String');
guidata(hObject,handles);

function editDate_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editNA_Callback(hObject, eventdata, handles)
handles.na = get(hObject,'String');
guidata(hObject,handles);

function editNA_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function pushbuttonLoadCSV_Callback(hObject, eventdata, handles)
if isempty(handles.dateformat)==1
    errordlg('Date format is not defined!','Error');
    return
elseif isempty(handles.na)==1
    errordlg('NA is not defined!','Error');
    return
end
importcsv(handles);


function popupmenuDailyStreamflow_Callback(hObject, eventdata, handles)
h = findobj('Tag','figureImportCSV');
if ~isempty(h)
    % get handles and other user-defined data associated
    figimport = guidata(h);
    handles.Qsyn = figimport.Qsyn;
end
str = get(hObject, 'String');
val = get(hObject,'Value');
switch str{val}
    case ''
        cla(handles.axes,'reset')
        cla(handles.axes1,'reset')
        cla(handles.axes2,'reset')
        cla(handles.axes3,'reset')
        cla(handles.axes4,'reset')
        cla(handles.axes5,'reset')
        set(handles.axes,'visible','on');
        set(handles.axes1,'visible','off');
        set(handles.axes2,'visible','off');
        set(handles.axes3,'visible','off');
        set(handles.axes4,'visible','off');
        set(handles.axes5,'visible','off');
    case 'Daily Streamflow'
        cla(handles.axes,'reset')
        cla(handles.axes1,'reset')
        cla(handles.axes2,'reset')
        cla(handles.axes3,'reset')
        cla(handles.axes4,'reset')
        cla(handles.axes5,'reset')
        set(handles.axes,'visible','on');
        set(handles.axes1,'visible','off');
        set(handles.axes2,'visible','off');
        set(handles.axes3,'visible','off');
        set(handles.axes4,'visible','off');
        set(handles.axes5,'visible','off');
        axes(handles.axes);
        plot(handles.axes,handles.Qsyn.tt_obs.Date, handles.Qsyn.tt_obs.Q,'blue');
        title('Daily Streamflow');
        grid;
        xlabel('Date [Days]','FontWeight','bold');
        ylabel('Q [m^3/s]','FontWeight','bold');
end

function popupmenuDailyStreamflow_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popupmenuSeas_Callback(hObject, eventdata, handles)
h = findobj('Tag','figureImportCSV');
if ~isempty(h)
    % get handles and other user-defined data associated
    figimport = guidata(h);
    handles.Qsyn = figimport.Qsyn;
end
str = get(hObject, 'String');
val = get(hObject,'Value');
switch str{val}
    case ''
        cla(handles.axes,'reset')
        cla(handles.axes1,'reset')
        cla(handles.axes2,'reset')
        cla(handles.axes3,'reset')
        cla(handles.axes4,'reset')
        cla(handles.axes5,'reset')
        set(handles.axes,'visible','on');
        set(handles.axes1,'visible','off');
        set(handles.axes2,'visible','off');
        set(handles.axes3,'visible','off');
        set(handles.axes4,'visible','off');
        set(handles.axes5,'visible','off');
    case 'Runoff Regime (daily)'
        cla(handles.axes,'reset')
        cla(handles.axes1,'reset')
        cla(handles.axes2,'reset')
        cla(handles.axes3,'reset')
        cla(handles.axes4,'reset')
        cla(handles.axes5,'reset')
        set(handles.axes,'visible','on');
        set(handles.axes1,'visible','off');
        set(handles.axes2,'visible','off');
        set(handles.axes3,'visible','off');
        set(handles.axes4,'visible','off');
        set(handles.axes5,'visible','off');
        axes(handles.axes);
        plot(handles.axes,handles.Qsyn.Q_regime_daily.DD, handles.Qsyn.Q_regime_daily.nanmean_Q, '-b','visible','on');
        title('Runoff Regime (daily)');
        grid;
        xlim([1 365]);
        md = nanmean(handles.Qsyn.Q_regime_daily.nanmean_Q);
        hline = refline([0 md]);
        hline.Color = 'blue';
        hline.LineStyle = '-.';
        ylim_up = ceil(max(handles.Qsyn.Q_regime_daily.nanmean_Q))+1;
        ylim([0,ylim_up]);
        xlabel('Day of Year','FontWeight','bold');
        ylabel('Q [m^3/s]','FontWeight','bold');
    case 'Runoff Regime (monthly)'
        cla(handles.axes,'reset')
        cla(handles.axes1,'reset')
        cla(handles.axes2,'reset')
        cla(handles.axes3,'reset')
        cla(handles.axes4,'reset')
        cla(handles.axes5,'reset')
        set(handles.axes,'visible','on');
        set(handles.axes1,'visible','off');
        set(handles.axes2,'visible','off');
        set(handles.axes3,'visible','off');
        set(handles.axes4,'visible','off');
        set(handles.axes5,'visible','off');
        axes(handles.axes);
        plot(handles.axes,handles.Qsyn.Q_regime_monthly.MM, handles.Qsyn.Q_regime_monthly.nanmean_Q, '-bo','visible','on');
        title('Runoff Regime (monthly)');
        grid;
        mm = nanmean(handles.Qsyn.Q_regime_monthly.nanmean_Q);
        hline = refline([0 mm]);
        hline.Color = 'blue';
        hline.LineStyle = '-.';
        xlim([1 12]);
        ylim_up = ceil(max(handles.Qsyn.Q_regime_monthly.nanmean_Q));
        ylim([0,ylim_up]);
        xlabel('Month','FontWeight','bold');
        ylabel('Q [m^3/s]','FontWeight','bold');
    case 'Parde Coefficient'
        cla(handles.axes,'reset')
        cla(handles.axes1,'reset')
        cla(handles.axes2,'reset')
        cla(handles.axes3,'reset')
        cla(handles.axes4,'reset')
        cla(handles.axes5,'reset')
        set(handles.axes,'visible','on');
        set(handles.axes1,'visible','off');
        set(handles.axes2,'visible','off');
        set(handles.axes3,'visible','off');
        set(handles.axes4,'visible','off');
        set(handles.axes5,'visible','off');
        axes(handles.axes);
        plot(handles.axes,handles.Qsyn.Q_regime_monthly.MM, handles.Qsyn.Q_regime_monthly.parde, '-bo','visible','on');
        title('Parde Coefficient');
        grid;
        hline = refline([0 1]);
        hline.Color = 'blue';
        hline.LineStyle = '-.';
        xlim([1 12]);
        ylim_up = ceil(max(handles.Qsyn.Q_regime_monthly.parde));
        ylim([0,ylim_up]);
        xlabel('Month','FontWeight','bold');
        ylabel('[-]','FontWeight','bold');
end

guidata(hObject,handles)

function popupmenuSeas_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editStartDate_Callback(hObject, eventdata, handles)
handles.startdate = get(hObject,'String');
guidata(hObject,handles);

function editStartDate_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editEndDate_Callback(hObject, eventdata, handles)
handles.enddate = get(hObject,'String');
guidata(hObject,handles);

function editEndDate_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function pushbuttonRun_Callback(hObject, eventdata, handles)
if isempty(handles.startdate)==1
    errordlg('Start date is not defined!','Error');
    return
elseif isempty(handles.enddate)==1
    errordlg('End date is not defined!','Error');
    return
end
startdate = datetime(handles.startdate,'InputFormat',handles.dateformat);
enddate = datetime(handles.enddate,'InputFormat',handles.dateformat);
startDD = day(startdate);
startMM = month(startdate);
endDD = day(enddate);
endMM = month(enddate);
if startDD~=1 && startMM~=1
    errordlg({'Start date is not a 1 January.' ...
        'In order to run the streamflow' ...
        'Generator complete years are necessary.'},'Error');
    return
elseif endDD~=31 && endMM~=12
    errordlg({'End date is not a 31 December.' ...
        'In order to run the streamflow' ...
        'Generator complete years are necessary.'},'Error');
    return
end
h = findobj('Tag','figureImportCSV');
if ~isempty(h)
    % get handles and other user-defined data associated
    figimport = guidata(h);
    handles.Qsyn = figimport.Qsyn;
end
h = waitbar(0,'Please wait...','Name','Generating streamflow...', ...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(h,'canceling',0);
opt_val = NaN(10,1);
handles.Qsyn.settimeperiod(handles.startdate, handles.enddate, handles.dateformat);
handles.Qsyn.rmseas(handles.Qsyn.Q_regime_daily, handles.Qsyn.Q_std_daily);
handles.Qsyn.determineorder();
handles.Qsyn.selectmodel(handles.Qsyn.tt_obs.Q_trans_stand_d);
Q_max = max(handles.Qsyn.tt_obs.Q)*1.1;
clusters = {1,handles.l}; % initializing clusters for iteration
for i = 1:handles.l
    if getappdata(h,'canceling')
        delete(h)
        delete(clusters)
        break
    end
    clusters{1,i} = Qsynth();
    clusters{1,i}.settimeperiod(handles.startdate, handles.enddate, handles.dateformat);
    clusters{1,i}.tt_obs = handles.Qsyn.tt_obs;
    clusters{1,i}.Q_regime_daily = handles.Qsyn.Q_regime_daily;
    clusters{1,i}.Q_std_daily = handles.Qsyn.Q_std_daily;
    clusters{1,i}.N_obs = handles.Qsyn.N_obs;
    clusters{1,i}.p = handles.Qsyn.p;
    clusters{1,i}.EstMdl = handles.Qsyn.EstMdl;
    clusters{1,i}.generaterunoff(clusters{1,i}.EstMdl,clusters{1,i}.Q_regime_daily, clusters{1,i}.Q_std_daily);
    clusters{1,i}.cutpeaks(Q_max,clusters{1,i}.EstMdl,clusters{1,i}.tt_syn);
    clusters{1,i}.optMAwMWS();
    clusters{1,i}.tt_syn = clusters{1,i}.tt_syn; 
    clusters{1,i}.N_sim = clusters{1,i}.N_sim; 
    clusters{1,i}.w = clusters{1,i}.w;
    clusters{1,i}.Qx = clusters{1,i}.Qx;
    clusters{1,i}.mws = clusters{1,i}.mws;
    opt_val(i,1) = clusters{1,i}.opt_val;
    waitbar(i/handles.l)
end
[val, idx] = min(opt_val);
handles.res = clusters{1,idx}; % assign cluster with best results to handles
handles.res.tt_syn.Q_sim_re = handles.res.MAwMWS(handles.res.tt_syn.Q_sim_re,handles.res.w,handles.res.Qx);
handles.res.testnorm();
handles.res.compareIHA();
handles.res.volume();
handles.res.teststats();
set(handles.textModel, 'String', ['AR(' num2str(handles.res.p) ')'],'FontWeight','bold');
delete(h)
msgbox('Streamflow generation was successful!');
% p = gcp;
% delete(p) % shutting parallel pool down
guidata(hObject,handles);

function popupmenuSyntheticStreamflow_Callback(hObject, eventdata, handles)
str = get(hObject, 'String');
val = get(hObject,'Value');
switch str{val}
    case ''
        cla(handles.axes,'reset')
        cla(handles.axes1,'reset')
        cla(handles.axes2,'reset')
        cla(handles.axes3,'reset')
        cla(handles.axes4,'reset')
        cla(handles.axes5,'reset')
        set(handles.axes,'visible','on');
        set(handles.axes1,'visible','off');
        set(handles.axes2,'visible','off');
        set(handles.axes3,'visible','off');
        set(handles.axes4,'visible','off');
        set(handles.axes5,'visible','off');
    case 'Observed & Synthetic Streamflow'
        cla(handles.axes,'reset')
        cla(handles.axes1,'reset')
        cla(handles.axes2,'reset')
        cla(handles.axes3,'reset')
        cla(handles.axes4,'reset')
        cla(handles.axes5,'reset')
        set(handles.axes,'visible','on');
        set(handles.axes1,'visible','off');
        set(handles.axes2,'visible','off');
        set(handles.axes3,'visible','off');
        set(handles.axes4,'visible','off');
        set(handles.axes5,'visible','off');
        axes(handles.axes);
        plot(handles.axes,handles.res.tt_obs.Date, handles.res.tt_obs.Q, 'b');
        hold on
        plot(handles.axes,handles.res.tt_syn.Date, handles.res.tt_syn.Q_sim_re, 'r');
        hold off
        title('Observed & Synthetic Streamflow');
        grid;
        xlabel('Date [Days]','FontWeight','bold');
        ylabel('Q [m^3/s]','FontWeight','bold');
        xtickformat('dd-MMM-yyyy');
        legend({'Q_{obs}','Q_{syn}'},'Box','off','FontSize',12);
    case 'Synthetic Streamflow'
        cla(handles.axes,'reset')
        cla(handles.axes1,'reset')
        cla(handles.axes2,'reset')
        cla(handles.axes3,'reset')
        cla(handles.axes4,'reset')
        cla(handles.axes5,'reset')
        set(handles.axes,'visible','on');
        set(handles.axes1,'visible','off');
        set(handles.axes2,'visible','off');
        set(handles.axes3,'visible','off');
        set(handles.axes4,'visible','off');
        set(handles.axes5,'visible','off');
        axes(handles.axes);
        plot(handles.axes,handles.res.tt_syn.Date, handles.res.tt_syn.Q_sim_re, 'r');
        title('Synthetic Streamflow');
        grid;
        xlabel('Date [Days]','FontWeight','bold');
        ylabel('Q [m^3/s]','FontWeight','bold');
        xtickformat('dd-MMM-yyyy');
        legend({'Q_{syn}'},'Box','off','FontSize',12);
end
guidata(hObject,handles)

function popupmenuSyntheticStreamflow_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function popupmenuFilter_Callback(hObject, eventdata, handles)
str = get(hObject, 'String');
val = get(hObject,'Value');
switch str{val}
    case ''
        cla(handles.axes,'reset')
        cla(handles.axes1,'reset')
        cla(handles.axes2,'reset')
        cla(handles.axes3,'reset')
        cla(handles.axes4,'reset')
        cla(handles.axes5,'reset')
        set(handles.axes,'visible','on');
        set(handles.axes1,'visible','off');
        set(handles.axes2,'visible','off');
        set(handles.axes3,'visible','off');
        set(handles.axes4,'visible','off');
        set(handles.axes5,'visible','off');
    case 'Window Size vs Streamflow'
        cla(handles.axes,'reset')
        cla(handles.axes1,'reset')
        cla(handles.axes2,'reset')
        cla(handles.axes3,'reset')
        cla(handles.axes4,'reset')
        cla(handles.axes5,'reset')
        set(handles.axes,'visible','on');
        set(handles.axes1,'visible','off');
        set(handles.axes2,'visible','off');
        set(handles.axes3,'visible','off');
        set(handles.axes4,'visible','off');
        set(handles.axes5,'visible','off');
        mws_1 = [0;handles.res.mws(1,2)+2]';
        mws_2 = [max(handles.res.tt_syn.Q_sim_re);0]';
        mws = [mws_1;handles.res.mws;mws_2];
        axes(handles.axes);
        stairs(mws(:,1), mws(:,2), 'k');
        title('Window Size vs Streamflow');
        grid;
        yticks(flipud(mws(:,2)));
        xlabel('Q [m^3/s]','FontWeight','bold');
        ylabel('Window Size [Days]','FontWeight','bold');
end
guidata(hObject,handles)

function popupmenuFilter_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popupmenuACF_Callback(hObject, eventdata, handles)
str = get(hObject, 'String');
val = get(hObject,'Value');
switch str{val}
    case ''
        cla(handles.axes,'reset')
        cla(handles.axes1,'reset')
        cla(handles.axes2,'reset')
        cla(handles.axes3,'reset')
        cla(handles.axes4,'reset')
        cla(handles.axes5,'reset')
        set(handles.axes,'visible','on');
        set(handles.axes1,'visible','off');
        set(handles.axes2,'visible','off');
        set(handles.axes3,'visible','off');
        set(handles.axes4,'visible','off');
        set(handles.axes5,'visible','off');
    case 'Observed & Synthetic Streamflow'
        cla(handles.axes,'reset')
        cla(handles.axes1,'reset')
        cla(handles.axes2,'reset')
        cla(handles.axes3,'reset')
        cla(handles.axes4,'reset')
        cla(handles.axes5,'reset')
        set(handles.axes,'visible','off');
        set(handles.axes1,'visible','on');
        set(handles.axes2,'visible','on');
        set(handles.axes3,'visible','off');
        set(handles.axes4,'visible','off');
        set(handles.axes5,'visible','off');
        axes(handles.axes1);
        autocorr(handles.res.tt_obs.Q, 100);
        title('ACF - Observed');
        xlim([0 100]);
        xlabel('Lag [Days]','FontWeight','bold');
        ylabel('ACF','FontWeight','bold');
        axes(handles.axes2);
        autocorr(handles.res.tt_syn.Q_sim_re, 100);
        title('ACF - Synthetic');
        xlim([0 100]);
        xlabel('Lag [Days]','FontWeight','bold');
        ylabel('ACF','FontWeight','bold');
    case 'January'
        cla(handles.axes,'reset')
        cla(handles.axes1,'reset')
        cla(handles.axes2,'reset')
        cla(handles.axes3,'reset')
        cla(handles.axes4,'reset')
        cla(handles.axes5,'reset')
        set(handles.axes,'visible','off');
        set(handles.axes1,'visible','on');
        set(handles.axes2,'visible','on');
        set(handles.axes3,'visible','off');
        set(handles.axes4,'visible','off');
        set(handles.axes5,'visible','off');
        axes(handles.axes1);
        autocorr(handles.res.tt_obs.Q(handles.res.tt_obs.MM==1), 100);
        title('ACF - Observed');
        xlim([0 100]);
        xlabel('Lag [Days]','FontWeight','bold');
        ylabel('ACF','FontWeight','bold');
        axes(handles.axes2);
        autocorr(handles.res.tt_syn.Q_sim_re(handles.res.tt_syn.MM==1), 100);
        title('ACF - Synthetic');
        xlim([0 100]);
        xlabel('Lag [Days]','FontWeight','bold');
        ylabel('ACF','FontWeight','bold');
    case 'February'
        cla(handles.axes,'reset')
        cla(handles.axes1,'reset')
        cla(handles.axes2,'reset')
        cla(handles.axes3,'reset')
        cla(handles.axes4,'reset')
        cla(handles.axes5,'reset')
        set(handles.axes,'visible','off');
        set(handles.axes1,'visible','on');
        set(handles.axes2,'visible','on');
        set(handles.axes3,'visible','off');
        set(handles.axes4,'visible','off');
        set(handles.axes5,'visible','off');
        axes(handles.axes1);
        autocorr(handles.res.tt_obs.Q(handles.res.tt_obs.MM==2), 100);
        title('ACF - Observed');
        xlim([0 100]);
        xlabel('Lag [Days]','FontWeight','bold');
        ylabel('ACF','FontWeight','bold');
        axes(handles.axes2);
        autocorr(handles.res.tt_syn.Q_sim_re(handles.res.tt_syn.MM==2), 100);
        title('ACF - Synthetic');
        xlim([0 100]);
        xlabel('Lag [Days]','FontWeight','bold');
        ylabel('ACF','FontWeight','bold');
    case 'March'
        cla(handles.axes,'reset')
        cla(handles.axes1,'reset')
        cla(handles.axes2,'reset')
        cla(handles.axes3,'reset')
        cla(handles.axes4,'reset')
        cla(handles.axes5,'reset')
        set(handles.axes,'visible','off');
        set(handles.axes1,'visible','on');
        set(handles.axes2,'visible','on');
        set(handles.axes3,'visible','off');
        set(handles.axes4,'visible','off');
        set(handles.axes5,'visible','off');
        axes(handles.axes1);
        autocorr(handles.res.tt_obs.Q(handles.res.tt_obs.MM==3), 100);
        title('ACF - Observed');
        xlim([0 100]);
        xlabel('Lag [Days]','FontWeight','bold');
        ylabel('ACF','FontWeight','bold');
        axes(handles.axes2);
        autocorr(handles.res.tt_syn.Q_sim_re(handles.res.tt_syn.MM==3), 100);
        title('ACF - Synthetic');
        xlim([0 100]);
        xlabel('Lag [Days]');
        ylabel('ACF');
    case 'April'
        cla(handles.axes,'reset')
        cla(handles.axes1,'reset')
        cla(handles.axes2,'reset')
        cla(handles.axes3,'reset')
        cla(handles.axes4,'reset')
        cla(handles.axes5,'reset')
        set(handles.axes,'visible','off');
        set(handles.axes1,'visible','on');
        set(handles.axes2,'visible','on');
        set(handles.axes3,'visible','off');
        set(handles.axes4,'visible','off');
        set(handles.axes5,'visible','off');
        axes(handles.axes1);
        autocorr(handles.res.tt_obs.Q(handles.res.tt_obs.MM==4), 100);
        title('ACF - Observed');
        xlim([0 100]);
        xlabel('Lag [Days]','FontWeight','bold');
        ylabel('ACF','FontWeight','bold');
        axes(handles.axes2);
        autocorr(handles.res.tt_syn.Q_sim_re(handles.res.tt_syn.MM==4), 100);
        title('ACF - Synthetic');
        xlim([0 100]);
        xlabel('Lag [Days]','FontWeight','bold');
        ylabel('ACF','FontWeight','bold');
    case 'May'
        cla(handles.axes,'reset')
        cla(handles.axes1,'reset')
        cla(handles.axes2,'reset')
        cla(handles.axes3,'reset')
        cla(handles.axes4,'reset')
        cla(handles.axes5,'reset')
        set(handles.axes,'visible','off');
        set(handles.axes1,'visible','on');
        set(handles.axes2,'visible','on');
        set(handles.axes3,'visible','off');
        set(handles.axes4,'visible','off');
        set(handles.axes5,'visible','off');
        axes(handles.axes1);
        autocorr(handles.res.tt_obs.Q(handles.res.tt_obs.MM==5), 100);
        title('ACF - Observed');
        xlim([0 100]);
        xlabel('Lag [Days]','FontWeight','bold');
        ylabel('ACF','FontWeight','bold');
        axes(handles.axes2);
        autocorr(handles.res.tt_syn.Q_sim_re(handles.res.tt_syn.MM==5), 100);
        title('ACF - Synthetic');
        xlim([0 100]);
        xlabel('Lag [Days]','FontWeight','bold');
        ylabel('ACF','FontWeight','bold');
    case 'June'
        cla(handles.axes,'reset')
        cla(handles.axes1,'reset')
        cla(handles.axes2,'reset')
        cla(handles.axes3,'reset')
        cla(handles.axes4,'reset')
        cla(handles.axes5,'reset')
        set(handles.axes,'visible','off');
        set(handles.axes1,'visible','on');
        set(handles.axes2,'visible','on');
        set(handles.axes3,'visible','off');
        set(handles.axes4,'visible','off');
        set(handles.axes5,'visible','off');
        axes(handles.axes1);
        autocorr(handles.res.tt_obs.Q(handles.res.tt_obs.MM==6), 100);
        title('ACF - Observed');
        xlim([0 100]);
        xlabel('Lag [Days]','FontWeight','bold');
        ylabel('ACF','FontWeight','bold');
        axes(handles.axes2);
        autocorr(handles.res.tt_syn.Q_sim_re(handles.res.tt_syn.MM==6), 100);
        title('ACF - Synthetic');
        xlim([0 100]);
        xlabel('Lag [Days]','FontWeight','bold');
        ylabel('ACF','FontWeight','bold');
    case 'July'
        cla(handles.axes,'reset')
        cla(handles.axes1,'reset')
        cla(handles.axes2,'reset')
        cla(handles.axes3,'reset')
        cla(handles.axes4,'reset')
        cla(handles.axes5,'reset')
        set(handles.axes,'visible','off');
        set(handles.axes1,'visible','on');
        set(handles.axes2,'visible','on');
        set(handles.axes3,'visible','off');
        set(handles.axes4,'visible','off');
        set(handles.axes5,'visible','off');
        axes(handles.axes1);
        autocorr(handles.res.tt_obs.Q(handles.res.tt_obs.MM==7), 100);
        title('ACF - Observed');
        xlim([0 100]);
        xlabel('Lag [Days]','FontWeight','bold');
        ylabel('ACF','FontWeight','bold');
        axes(handles.axes2);
        autocorr(handles.res.tt_syn.Q_sim_re(handles.res.tt_syn.MM==7), 100);
        title('ACF - Synthetic');
        xlim([0 100]);
        xlabel('Lag [Days]','FontWeight','bold');
        ylabel('ACF','FontWeight','bold');
    case 'August'
        cla(handles.axes,'reset')
        cla(handles.axes1,'reset')
        cla(handles.axes2,'reset')
        cla(handles.axes3,'reset')
        cla(handles.axes4,'reset')
        cla(handles.axes5,'reset')
        set(handles.axes,'visible','off');
        set(handles.axes1,'visible','on');
        set(handles.axes2,'visible','on');
        set(handles.axes3,'visible','off');
        set(handles.axes4,'visible','off');
        set(handles.axes5,'visible','off');
        axes(handles.axes1);
        autocorr(handles.res.tt_obs.Q(handles.res.tt_obs.MM==8), 100);
        title('ACF - Observed');
        xlim([0 100]);
        xlabel('Lag [Days]','FontWeight','bold');
        ylabel('ACF','FontWeight','bold');
        axes(handles.axes2);
        autocorr(handles.res.tt_syn.Q_sim_re(handles.res.tt_syn.MM==8), 100);
        title('ACF - Synthetic');
        xlim([0 100]);
        xlabel('Lag [Days]','FontWeight','bold');
        ylabel('ACF','FontWeight','bold');
    case 'September'
        cla(handles.axes,'reset')
        cla(handles.axes1,'reset')
        cla(handles.axes2,'reset')
        cla(handles.axes3,'reset')
        cla(handles.axes4,'reset')
        cla(handles.axes5,'reset')
        set(handles.axes,'visible','off');
        set(handles.axes1,'visible','on');
        set(handles.axes2,'visible','on');
        set(handles.axes3,'visible','off');
        set(handles.axes4,'visible','off');
        set(handles.axes5,'visible','off');
        axes(handles.axes1);
        autocorr(handles.res.tt_obs.Q(handles.res.tt_obs.MM==9), 100);
        title('ACF - Observed');
        xlim([0 100]);
        xlabel('Lag [Days]','FontWeight','bold');
        ylabel('ACF','FontWeight','bold');
        axes(handles.axes2);
        autocorr(handles.res.tt_syn.Q_sim_re(handles.res.tt_syn.MM==9), 100);
        title('ACF - Synthetic');
        xlim([0 100]);
        xlabel('Lag [Days]','FontWeight','bold');
        ylabel('ACF','FontWeight','bold');
    case 'October'
        cla(handles.axes,'reset')
        cla(handles.axes1,'reset')
        cla(handles.axes2,'reset')
        cla(handles.axes3,'reset')
        cla(handles.axes4,'reset')
        cla(handles.axes5,'reset')
        set(handles.axes,'visible','off');
        set(handles.axes1,'visible','on');
        set(handles.axes2,'visible','on');
        set(handles.axes3,'visible','off');
        set(handles.axes4,'visible','off');
        set(handles.axes5,'visible','off');
        axes(handles.axes1);
        autocorr(handles.res.tt_obs.Q(handles.res.tt_obs.MM==10), 100);
        title('ACF - Observed');
        xlim([0 100]);
        xlabel('Lag [Days]','FontWeight','bold');
        ylabel('ACF','FontWeight','bold');
        axes(handles.axes2);
        autocorr(handles.res.tt_syn.Q_sim_re(handles.res.tt_syn.MM==10), 100);
        title('ACF - Synthetic');
        xlim([0 100]);
        xlabel('Lag [Days]','FontWeight','bold');
        ylabel('ACF','FontWeight','bold');
    case 'ACF - November'
        cla(handles.axes,'reset')
        cla(handles.axes1,'reset')
        cla(handles.axes2,'reset')
        cla(handles.axes3,'reset')
        cla(handles.axes4,'reset')
        cla(handles.axes5,'reset')
        set(handles.axes,'visible','off');
        set(handles.axes1,'visible','on');
        set(handles.axes2,'visible','on');
        set(handles.axes3,'visible','off');
        set(handles.axes4,'visible','off');
        set(handles.axes5,'visible','off');
        axes(handles.axes1);
        autocorr(handles.res.tt_obs.Q(handles.res.tt_obs.MM==11), 100);
        title('ACF - Observed');
        xlim([0 100]);
        xlabel('Lag [Days]');
        ylabel('ACF');
        axes(handles.axes2);
        autocorr(handles.res.tt_syn.Q_sim_re(handles.res.tt_syn.MM==11), 100);
        title('ACF - Synthetic');
        xlim([0 100]);
        xlabel('Lag [Days]','FontWeight','bold');
        ylabel('ACF','FontWeight','bold');
    case 'December'
        cla(handles.axes,'reset')
        cla(handles.axes1,'reset')
        cla(handles.axes2,'reset')
        cla(handles.axes3,'reset')
        cla(handles.axes4,'reset')
        cla(handles.axes5,'reset')
        set(handles.axes,'visible','off');
        set(handles.axes1,'visible','on');
        set(handles.axes2,'visible','on');
        set(handles.axes3,'visible','off');
        set(handles.axes4,'visible','off');
        set(handles.axes5,'visible','off');
        axes(handles.axes1);
        autocorr(handles.res.tt_obs.Q(handles.res.tt_obs.MM==12), 100);
        title('ACF - Observed');
        xlim([0 100]);
        xlabel('Lag [Days]','FontWeight','bold');
        ylabel('ACF','FontWeight','bold');
        axes(handles.axes2);
        autocorr(handles.res.tt_syn.Q_sim_re(handles.res.tt_syn.MM==12), 100);
        title('ACF - Synthetic');
        xlim([0 100]);
        xlabel('Lag [Days]','FontWeight','bold');
        ylabel('ACF','FontWeight','bold');
    case 'Residuals'
        cla(handles.axes,'reset')
        cla(handles.axes1,'reset')
        cla(handles.axes2,'reset')
        cla(handles.axes3,'reset')
        cla(handles.axes4,'reset')
        cla(handles.axes5,'reset')
        set(handles.axes,'visible','on');
        set(handles.axes1,'visible','off');
        set(handles.axes2,'visible','off');
        set(handles.axes3,'visible','off');
        set(handles.axes4,'visible','off');
        set(handles.axes5,'visible','off');
        axes(handles.axes);
        autocorr(handles.res.tt_syn.E, 100);
        xlim([0 100]);
        title('ACF - Residuals');
        xlabel('Lag [Days]','FontWeight','bold');
        ylabel('ACF','FontWeight','bold');
end
guidata(hObject,handles)

function popupmenuACF_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popupmenuHistogram_Callback(hObject, eventdata, handles)
str = get(hObject, 'String');
val = get(hObject,'Value');
switch str{val}
    case ''
        cla(handles.axes,'reset')
        cla(handles.axes1,'reset')
        cla(handles.axes2,'reset')
        cla(handles.axes3,'reset')
        cla(handles.axes4,'reset')
        cla(handles.axes5,'reset')
        set(handles.axes,'visible','on');
        set(handles.axes1,'visible','off');
        set(handles.axes2,'visible','off');
        set(handles.axes3,'visible','off');
        set(handles.axes4,'visible','off');
        set(handles.axes5,'visible','off');
        plot(handles.axes,handles.nan(:,1), handles.nan(:,2));
    case 'Observed & Synthetic Streamflow'
        cla(handles.axes,'reset')
        cla(handles.axes1,'reset')
        cla(handles.axes2,'reset')
        cla(handles.axes3,'reset')
        cla(handles.axes4,'reset')
        cla(handles.axes5,'reset')
        set(handles.axes,'visible','off');
        set(handles.axes1,'visible','on');
        set(handles.axes2,'visible','on');
        set(handles.axes3,'visible','off');
        set(handles.axes4,'visible','off');
        set(handles.axes5,'visible','off');
        axes(handles.axes1);
        histogram(handles.res.tt_obs.Q,handles.res.nbins_nn,'Normalization','probability','FaceColor','b');
        grid; title('Histogram - Observed');
        xlim([handles.res.x1_nn handles.res.x2_nn]);
        ylim([0 handles.res.y2_nn]);
        xlabel('Q [m^3/s]','FontWeight','bold');
        ylabel('Probability','FontWeight','bold');
        axes(handles.axes2);
        histogram(handles.res.tt_syn.Q_sim_re,handles.res.nbins_nn,'Normalization','probability','FaceColor','r');
        grid; title('Histogram - Synthetic');
        xlim([handles.res.x1_nn handles.res.x2_nn]);
        ylim([0 handles.res.y2_nn]);
        xlabel('Q [m^3/s]','FontWeight','bold');
        ylabel('Probability','FontWeight','bold');
    case 'Transformed-to-normal & Standardized'
        cla(handles.axes,'reset')
        cla(handles.axes1,'reset')
        cla(handles.axes2,'reset')
        cla(handles.axes3,'reset')
        cla(handles.axes4,'reset')
        cla(handles.axes5,'reset')
        set(handles.axes,'visible','off');
        set(handles.axes1,'visible','off');
        set(handles.axes2,'visible','off');
        set(handles.axes3,'visible','on');
        set(handles.axes4,'visible','on');
        set(handles.axes5,'visible','on');
        axes(handles.axes3);
        histogram(handles.res.tt_obs.Q_trans,handles.res.nbins_n,'Normalization','probability','FaceColor','b');
        grid; title('Histogram - Transformed-to-normal (Observed)');
        xlim([handles.res.x1_n handles.res.x2_n]);
        ylim([0 handles.res.y2_n]);
        axes(handles.axes4);
        histogram(handles.res.tt_obs.Q_trans_stand_d,handles.res.nbins_n,'Normalization','probability','FaceColor','b');
        grid; title('Histogram - Transformed-to-normal & Standardized (Observed)');
        xlim([handles.res.x1_n handles.res.x2_n]);
        ylim([0 handles.res.y2_n]);
        ylabel('Probability','FontWeight','bold');
        axes(handles.axes5);
        histogram(handles.res.tt_syn.Q_sim,handles.res.nbins_n,'Normalization','probability','FaceColor','r');
        grid; title('Histogram - Transformed-to-normal & Standardized (Synthetic)');
        xlim([handles.res.x1_n handles.res.x2_n]);
        ylim([0 handles.res.y2_n]);
        xlabel('Q [m^3/s]','FontWeight','bold');
    case 'Residuals'
        cla(handles.axes,'reset')
        cla(handles.axes1,'reset')
        cla(handles.axes2,'reset')
        cla(handles.axes3,'reset')
        cla(handles.axes4,'reset')
        cla(handles.axes5,'reset')
        set(handles.axes,'visible','on');
        set(handles.axes1,'visible','off');
        set(handles.axes2,'visible','off');
        set(handles.axes3,'visible','off');
        set(handles.axes4,'visible','off');
        set(handles.axes5,'visible','off');
        axes(handles.axes);
        histogram(handles.res.tt_syn.E,'FaceColor','k','Normalization','probability');
        grid; title('Histogram - Residuals');
        ylabel('Probability','FontWeight','bold');
end
guidata(hObject,handles)

function popupmenuHistogram_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popupmenuIHA_Callback(hObject, eventdata, handles)
str = get(hObject, 'String');
val = get(hObject,'Value');
switch str{val}
    case ''
        cla(handles.axes,'reset')
        cla(handles.axes1,'reset')
        cla(handles.axes2,'reset')
        cla(handles.axes3,'reset')
        cla(handles.axes4,'reset')
        cla(handles.axes5,'reset')
        set(handles.axes,'visible','on');
        set(handles.axes1,'visible','off');
        set(handles.axes2,'visible','off');
        set(handles.axes3,'visible','off');
        set(handles.axes4,'visible','off');
        set(handles.axes5,'visible','off');
        plot(handles.axes,handles.nan(:,1), handles.nan(:,2));
    case 'Group 1'
        cla(handles.axes,'reset')
        cla(handles.axes1,'reset')
        cla(handles.axes2,'reset')
        cla(handles.axes3,'reset')
        cla(handles.axes4,'reset')
        cla(handles.axes5,'reset')
        set(handles.axes,'visible','on');
        set(handles.axes1,'visible','off');
        set(handles.axes2,'visible','off');
        set(handles.axes3,'visible','off');
        set(handles.axes4,'visible','off');
        set(handles.axes5,'visible','off');
        axes(handles.axes);
        hold on
        b = bar([handles.res.IHA_ind_obs_mean(1:12) handles.res.IHA_ind_syn_mean(1:12)]);
        title('IHA - Group 1');
        grid;
        b(1).FaceColor = 'b';
        b(2).FaceColor = 'r';
        ylabel('Q [m^3/s]','FontWeight','bold');
        legend({'Q_{obs}','Q_{syn}'},'Box','off','FontSize',12)
        xticks(1:1:12);
    case 'Group 2'
        cla(handles.axes,'reset')
        cla(handles.axes1,'reset')
        cla(handles.axes2,'reset')
        cla(handles.axes3,'reset')
        cla(handles.axes4,'reset')
        cla(handles.axes5,'reset')
        set(handles.axes,'visible','on');
        set(handles.axes1,'visible','off');
        set(handles.axes2,'visible','off');
        set(handles.axes3,'visible','off');
        set(handles.axes4,'visible','off');
        set(handles.axes5,'visible','off');
        axes(handles.axes);
        hold on
        b = bar([handles.res.IHA_ind_obs_mean(13:22) handles.res.IHA_ind_syn_mean(13:22)]);
        title('IHA - Group 2');
        grid;
        b(1).FaceColor = 'b';
        b(2).FaceColor = 'r';
        ylabel('Q [m^3/s]','FontWeight','bold');
        legend({'Q_{obs}','Q_{syn}'},'Box','off','FontSize',12)
        xp = 1:1:10;
        xt = {{'Min'; '1-day'} {'Min'; '3-day'} ...
            {'Min'; '7-day'} {'Min'; '30-day'}...
            {'Min'; '90-day'} {'Max'; '1-day'}...
            {'Max'; '3-day'} {'Max'; '7-day'}...
            {'Max'; '30-day'} {'Max'; '90-day'}};
        my_xticklabels(xp, xt);
    case 'Group 4'
        cla(handles.axes,'reset')
        cla(handles.axes1,'reset')
        cla(handles.axes2,'reset')
        cla(handles.axes3,'reset')
        cla(handles.axes4,'reset')
        cla(handles.axes5,'reset')
        set(handles.axes,'visible','on');
        set(handles.axes1,'visible','off');
        set(handles.axes2,'visible','off');
        set(handles.axes3,'visible','off');
        set(handles.axes4,'visible','off');
        set(handles.axes5,'visible','off');
        axes(handles.axes);
        hold on
        b = bar([handles.res.IHA_ind_obs_mean(27:30) handles.res.IHA_ind_syn_mean(27:30)]);
        title('IHA - Group 4');
        grid;
        b(1).FaceColor = 'b';
        b(2).FaceColor = 'r';
        legend({'Q_{obs}','Q_{syn}'},'Box','off','FontSize',12)
        xp = 1:1:4;
        xt = {{'No. of low pulses'} {'No. of high pulses'} {'Mean duration'; 'of low pulses'} {'Mean duration'; 'of high pulses'}};
        my_xticklabels(xp, xt);
    case 'Group 5-1'
        cla(handles.axes,'reset')
        cla(handles.axes1,'reset')
        cla(handles.axes2,'reset')
        cla(handles.axes3,'reset')
        cla(handles.axes4,'reset')
        cla(handles.axes5,'reset')
        set(handles.axes,'visible','on');
        set(handles.axes1,'visible','off');
        set(handles.axes2,'visible','off');
        set(handles.axes3,'visible','off');
        set(handles.axes4,'visible','off');
        set(handles.axes5,'visible','off');
        axes(handles.axes);
        hold on
        b = bar([handles.res.IHA_ind_obs_mean(31:32) handles.res.IHA_ind_syn_mean(31:32)]);
        title('IHA - Group 5-1');
        grid;
        b(1).FaceColor = 'b';
        b(2).FaceColor = 'r';
        legend({'Q_{obs}','Q_{syn}'},'Box','off','FontSize',12)
        xp = [1 2];
        xt = {{'Means of all negative'; 'differences between'; 'consecutive daily values'} ...
            {'Means of all postive'; 'differences between'; 'consecutive daily means'}};
        my_xticklabels(xp, xt);
    case 'Group 5-2'
        cla(handles.axes,'reset')
        cla(handles.axes1,'reset')
        cla(handles.axes2,'reset')
        cla(handles.axes3,'reset')
        cla(handles.axes4,'reset')
        cla(handles.axes5,'reset')
        set(handles.axes,'visible','on');
        set(handles.axes1,'visible','off');
        set(handles.axes2,'visible','off');
        set(handles.axes3,'visible','off');
        set(handles.axes4,'visible','off');
        set(handles.axes5,'visible','off');
        axes(handles.axes);
        hold on
        b = bar([handles.res.IHA_ind_obs_mean(33:34) handles.res.IHA_ind_syn_mean(33:34)]);
        title('IHA - Group 5-2');
        grid;
        b(1).FaceColor = 'b';
        b(2).FaceColor = 'r';
        legend({'Q_{obs}','Q_{syn}'},'Box','off','FontSize',12,'Location','north')
        xticks([1 2])
        xticklabels({'No. of falls', 'No. of rises'})
end
guidata(hObject,handles)

function popupmenuIHA_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popupmenuTestStats_Callback(hObject, eventdata, handles)
str = get(hObject, 'String');
val = get(hObject,'Value');
switch str{val}
    case ''
        cla(handles.axes,'reset')
        cla(handles.axes1,'reset')
        cla(handles.axes2,'reset')
        cla(handles.axes3,'reset')
        cla(handles.axes4,'reset')
        cla(handles.axes5,'reset')
        set(handles.axes,'visible','on');
        set(handles.axes1,'visible','off');
        set(handles.axes2,'visible','off');
        set(handles.axes3,'visible','off');
        set(handles.axes4,'visible','off');
        set(handles.axes5,'visible','off');
    case 'Min, Max, Mean, Std & Skew'
        cla(handles.axes,'reset')
        cla(handles.axes1,'reset')
        cla(handles.axes2,'reset')
        cla(handles.axes3,'reset')
        cla(handles.axes4,'reset')
        cla(handles.axes5,'reset')
        set(handles.axes,'visible','on');
        set(handles.axes1,'visible','off');
        set(handles.axes2,'visible','off');
        set(handles.axes3,'visible','off');
        set(handles.axes4,'visible','off');
        set(handles.axes5,'visible','off');
        axes(handles.axes);
        b = bar([handles.res.test_stats(:,1) handles.res.test_stats(:,2)]);
        title('Simple Test Statistic');
        b(1).FaceColor = 'b';
        b(2).FaceColor = 'r';
        grid;
        legend({'Q_{obs}','Q_{syn}'},'Box','off','Location','northwest','FontSize',12)
        ylabel('Q [m^3/s]','FontWeight','bold');
        xticks([1 2 3 4 5])
        xticklabels({'Min', 'Max', 'Mean', 'Std', 'Skew'})
end
guidata(hObject,handles)

function popupmenuTestStats_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popupmenuVolume_Callback(hObject, eventdata, handles)
str = get(hObject, 'String');
val = get(hObject,'Value');
switch str{val}
    case ''
        cla(handles.axes,'reset')
        cla(handles.axes1,'reset')
        cla(handles.axes2,'reset')
        cla(handles.axes3,'reset')
        cla(handles.axes4,'reset')
        cla(handles.axes5,'reset')
        set(handles.axes,'visible','on');
        set(handles.axes1,'visible','off');
        set(handles.axes2,'visible','off');
        set(handles.axes3,'visible','off');
        set(handles.axes4,'visible','off');
        set(handles.axes5,'visible','off');
    case 'Cumulated Volume (annually)'
        cla(handles.axes,'reset')
        cla(handles.axes1,'reset')
        cla(handles.axes2,'reset')
        cla(handles.axes3,'reset')
        cla(handles.axes4,'reset')
        cla(handles.axes5,'reset')
        set(handles.axes,'visible','on');
        set(handles.axes1,'visible','off');
        set(handles.axes2,'visible','off');
        set(handles.axes3,'visible','off');
        set(handles.axes4,'visible','off');
        set(handles.axes5,'visible','off');
        axes(handles.axes);
        hold on
        plot(handles.res.tt_vol_yy.Time, handles.res.tt_vol_yy.vol_obs_cum, 'b');
        plot(handles.res.tt_vol_yy.Time, handles.res.tt_vol_yy.vol_syn_cum, 'r');
        hold off
        title('Volume - Cumulated Volume (annually)');
        grid;
        xlabel('Date','FontWeight','bold');
        ylabel('Volume [km^3]','FontWeight','bold');
        xtickformat('yyyy');
        legend({'Q_{obs}','Q_{syn}'},'Box','off','Location','northwest','FontSize',12);
    case 'Cumulated Volume (monthly)'
        cla(handles.axes,'reset')
        cla(handles.axes1,'reset')
        cla(handles.axes2,'reset')
        cla(handles.axes3,'reset')
        cla(handles.axes4,'reset')
        cla(handles.axes5,'reset')
        set(handles.axes,'visible','on');
        set(handles.axes1,'visible','off');
        set(handles.axes2,'visible','off');
        set(handles.axes3,'visible','off');
        set(handles.axes4,'visible','off');
        set(handles.axes5,'visible','off');
        axes(handles.axes);
        hold on
        plot(handles.res.tt_vol_mm.Time, handles.res.tt_vol_mm.vol_obs_cum, 'b');
        plot(handles.res.tt_vol_mm.Time, handles.res.tt_vol_mm.vol_syn_cum, 'r');
        hold off
        title('Volume - Cumulated Volume (monthly)');
        grid;
        xlabel('Date','FontWeight','bold');
        ylabel('Volume [km^3]','FontWeight','bold');
        xtickformat('MMM-yyyy');
        legend({'Q_{obs}','Q_{syn}'},'Box','off','Location','northwest','FontSize',12);
end
guidata(hObject,handles)

function popupmenuVolume_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function pushbuttonExport_Callback(hObject, eventdata, handles)
export(handles)
