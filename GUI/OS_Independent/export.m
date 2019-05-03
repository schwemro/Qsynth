function varargout = export(varargin)
% EXPORT MATLAB code for export.fig
%      EXPORT, by itself, creates a new EXPORT or raises the existing
%      singleton*.
%
%      H = EXPORT returns the handle to a new EXPORT or the handle to
%      the existing singleton*.
%
%      EXPORT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EXPORT.M with the given input arguments.
%
%      EXPORT('Property','Value',...) creates a new EXPORT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before export_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to export_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help export

% Last Modified by GUIDE v2.5 28-Jul-2017 15:50:37

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @export_OpeningFcn, ...
                   'gui_OutputFcn',  @export_OutputFcn, ...
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


% --- Executes just before export is made visible.
function export_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to export (see VARARGIN)

% Choose default command line output for export
h = findobj('Tag','figureMain');
% if exists (not empty)
if ~isempty(h)
    % get handles and other user-defined data associated
    figmain = guidata(h);
    handles.Qsyn = figmain.Qsyn;
    handles.res = figmain.res;
end

handles.output = hObject;
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes export wait for user response (see UIRESUME)
% uiwait(handles.figureExport);


% --- Outputs from this function are returned to the command line.
function varargout = export_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function editFolder_Callback(hObject, eventdata, handles)
handles.path = get(hObject,'String');
guidata(hObject, handles);

function editFolder_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function pushbuttonBrowse_Callback(hObject, eventdata, handles)
folder_name = uigetdir;
set(handles.editFolder,'String',folder_name);
handles.path = get(handles.editFolder,'String');
guidata(hObject, handles);

function pushbuttonOK_Callback(hObject, eventdata, handles)
delete(handles.figureExport)
h = waitbar(0,'Please wait...','Name','Exporting results...', ...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(h,'canceling',0);
steps = 33;
f = figure('Name','Daily Streamflow - Daily Streamflow','NumberTitle','off', 'Visible','off');
plot(handles.Qsyn.tt_obs.Date, handles.Qsyn.tt_obs.Q,'blue');
grid;
xlabel('Date [Days]');
ylabel('Q [m^3/s]');
saveas(f,[handles.path '/daily_streamflow.pdf']);
if getappdata(h,'canceling')
    delete(h)
    return
end
step = 1;
waitbar(step / steps)

f = figure('Name','Seasonality - Runoff Regime (daily)','NumberTitle','off', 'Visible','off');
plot(handles.Qsyn.Q_regime_daily.DD, handles.Qsyn.Q_regime_daily.nanmean_Q, '-b');
grid;
xlim([1 365]);
md = nanmean(handles.Qsyn.Q_regime_daily.nanmean_Q);
hline = refline([0 md]);
hline.Color = 'blue';
hline.LineStyle = '-.';
ylim_up = ceil(max(handles.Qsyn.Q_regime_daily.nanmean_Q))+1;
ylim([0,ylim_up]);
xlabel('Day of Year');
ylabel('Q [m^3/s]');
saveas(f,[handles.path '/runoff_regime_daily.pdf']);
if getappdata(h,'canceling')
    delete(h)
    return
end
step = 2;
waitbar(step / steps)

f = figure('Name','Seasonality - Runoff Regime (monthly)','NumberTitle','off', 'Visible','off');
plot(handles.Qsyn.Q_regime_monthly.MM, handles.Qsyn.Q_regime_monthly.nanmean_Q, '-bo');
grid;
mm = nanmean(handles.Qsyn.Q_regime_monthly.nanmean_Q);
hline = refline([0 mm]);
hline.Color = 'blue';
hline.LineStyle = '-.';
xlim([1 12]);
ylim_up = ceil(max(handles.Qsyn.Q_regime_monthly.nanmean_Q));
ylim([0,ylim_up]);
xlabel('Month');
ylabel('Q [m^3/s]');
saveas(f,[handles.path '/runoff_regime_monthly.pdf']);
if getappdata(h,'canceling')
    delete(h)
    return
end
step = 3;
waitbar(step / steps)

f = figure('Name','Seasonality - Parde Coefficient','NumberTitle','off', 'Visible','off');
plot(handles.Qsyn.Q_regime_monthly.MM, handles.Qsyn.Q_regime_monthly.parde, '-bo');
grid;
hline = refline([0 1]);
hline.Color = 'blue';
hline.LineStyle = '-.';
xlim([1 12]);
ylim_up = ceil(max(handles.Qsyn.Q_regime_monthly.parde));
ylim([0,ylim_up]);
xlabel('Month');
ylabel('[-]');
saveas(f,[handles.path '/parde_coefficient.pdf']);
if getappdata(h,'canceling')
    delete(h)
    return
end
step = 4;
waitbar(step / steps)

f = figure('Name','Synthetic Streamflow - Observed & Synthetic Streamflow','NumberTitle','off','Visible','off');
plot(handles.res.tt_obs.Date, handles.res.tt_obs.Q, 'b');
hold on
plot(handles.res.tt_syn.Date, handles.res.tt_syn.Q_sim_re, 'r');
hold off
grid;
xlabel('Date [Days]');
ylabel('Q [m^3/s]');
xtickformat('dd-MMM-yyyy');
legend({'Q_{obs}','Q_{syn}'},'Box','off');
saveas(f,[handles.path '/Q_obs_vs_Q_syn.pdf']);
if getappdata(h,'canceling')
    delete(h)
    return
end
step = 5;
waitbar(step / steps)

f = figure('Name','Synthetic Streamflow - Synthetic Streamflow','NumberTitle','off','Visible','off');
plot(handles.res.tt_syn.Date, handles.res.tt_syn.Q_sim_re, 'r');
grid;
xlabel('Date [Days]');
ylabel('Q [m^3/s]');
xtickformat('dd-MMM-yyyy');
legend({'Q_{syn}'},'Box','off');
saveas(f,[handles.path '/Q_syn.pdf']);
if getappdata(h,'canceling')
    delete(h)
    return
end
step = 6;
waitbar(step / steps)

f = figure('Name','Filter - Window Size vs Streamflow','NumberTitle','off','Visible','off');
mws_1 = [0;handles.res.mws(1,2)+2]';
mws_2 = [max(handles.res.tt_syn.Q_sim_re);0]';
mws = [mws_1;handles.res.mws;mws_2];
stairs(mws(:,1), mws(:,2), 'k');
grid;
yticks(flipud(mws(:,2)));
xlabel('Q [m^3/s]');
ylabel('Window Size [Days]');
saveas(f,[handles.path '/MAwMWS.pdf'])
if getappdata(h,'canceling')
    delete(h)
    return
end
step = 7;
waitbar(step / steps)

f = figure('Name','Autocorrelation - Observed & Synthetic Streamflow','NumberTitle','off', 'Visible','off');
subplot(2,1,1);
autocorr(handles.res.tt_obs.Q, 100);
title('Observed');
xlim([0 100]);
xlabel('Lag [Days]');
ylabel('ACF');
subplot(2,1,2);
autocorr(handles.res.tt_syn.Q_sim_re, 100);
title('Synthetic');
xlim([0 100]);
xlabel('Lag [Days]');
ylabel('ACF');
saveas(f,[handles.path '/acf_of_obs_and_syn.pdf']);
if getappdata(h,'canceling')
    delete(h)
    return
end
step = 8;
waitbar(step / steps)

for i = 1:12
    f = figure('Name',['Autocorrelation - ' num2str(i)],'NumberTitle','off','Visible','off');
    subplot(2,1,1);
    autocorr(handles.res.tt_obs.Q(handles.res.tt_obs.MM==i), 100);
    title('Observed');
    xlim([0 100]);
    xlabel('Lag [Days]');
    ylabel('ACF');
    subplot(2,1,2);
    autocorr(handles.res.tt_syn.Q_sim_re(handles.res.tt_syn.MM==i), 100);
    title('Synthetic');
    xlim([0 100]);
    xlabel('Lag [Days]');
    ylabel('ACF');
    saveas(f,[handles.path '/acf_obs_and_syn_' num2str(i) '.pdf']);
    if getappdata(h,'canceling')
        delete(h)
        return
    end
    step = 8+i;
    waitbar(step / steps)
end

f = figure('Name','Autocorrelation - Residuals','NumberTitle','off', 'Visible','off');
autocorr(handles.res.tt_syn.E, 100);
xlim([0 100]);
title('');
xlabel('Lag [Days]');
ylabel('ACF');
saveas(f,[handles.path '/acf_residuals.pdf']);
if getappdata(h,'canceling')
    delete(h)
    return
end
step = 21;
waitbar(step / steps)

f = figure('Name','Histograms - Observed & Synthetic Streamflow','NumberTitle','off', 'Visible','off');
s1 = subplot(2, 1, 1);
histogram(handles.res.tt_obs.Q,handles.res.nbins_nn,'Normalization','probability','FaceColor','b');
grid; title('Observed');
xlim([handles.res.x1_nn handles.res.x2_nn]);
ylim([0 handles.res.y2_nn]);
s2 = subplot(2, 1, 2);
histogram(handles.res.tt_syn.Q_sim_re,handles.res.nbins_nn,'Normalization','probability','FaceColor','r');
grid; title('Synthetic');
xlim([handles.res.x1_nn handles.res.x2_nn]);
ylim([0 handles.res.y2_nn]);
xlabel('Q [m^3/s]');
p1=get(s1,'position');
p2=get(s2,'position');
height=p1(2)+p1(4)-p2(2);
axes('position',[p2(1) p2(2) p2(3) height],'visible','off');
ylabel('Probability','visible','on');
saveas(f,[handles.path '/histograms_obs_and_syn.pdf'])
if getappdata(h,'canceling')
    delete(h)
    return
end
step = 22;
waitbar(step / steps)

f = figure('Name','Histograms - Transformed-to-normal & Standardized','NumberTitle','off', 'Visible','off');
subplot(3, 1, 1);
histogram(handles.res.tt_obs.Q_trans,handles.res.nbins_n,'Normalization','probability','FaceColor','b');
grid; title('Transformed-to-normal');
xlim([handles.res.x1_n handles.res.x2_n]);
ylim([0 handles.res.y2_n]);
subplot(3, 1, 2);
histogram(handles.res.tt_obs.Q_trans_stand_d,handles.res.nbins_n,'Normalization','probability','FaceColor','b');
grid; title('Transformed-to-normal & Standardized (Observed)');
xlim([handles.res.x1_n handles.res.x2_n]);
ylim([0 handles.res.y2_n]);
ylabel('Probability');
subplot(3, 1, 3);
histogram(handles.res.tt_syn.Q_sim,handles.res.nbins_n,'Normalization','probability','FaceColor','r');
grid; title('Transformed-to-normal & Standardized (Synthetic)');
xlim([handles.res.x1_n handles.res.x2_n]);
ylim([0 handles.res.y2_n]);
xlabel('Q [m^3/s]');
saveas(f,[handles.path '/histograms_trans_to_normal_and_standardized.pdf'])
if getappdata(h,'canceling')
    delete(h)
    return
end
step = 23;
waitbar(step / steps)

f = figure('Name','Histograms - Residuals','NumberTitle','off', 'Visible','off');
histogram(handles.res.tt_syn.E,'FaceColor','k','Normalization','probability');
grid;
ylabel('Probability');
saveas(f,[handles.path '/histogram_residuals.pdf'])
if getappdata(h,'canceling')
    delete(h)
    return
end
step = 24;
waitbar(step / steps)

f = figure('Name','IHA - Group 1','NumberTitle','off', 'Visible','off');
hold on
b = bar([handles.res.IHA_ind_obs_mean(1:12) handles.res.IHA_ind_syn_mean(1:12)]);
grid;
b(1).FaceColor = 'b';
b(2).FaceColor = 'r';
ylabel('Q [m^3/s]');
legend({'Q_{obs}','Q_{syn}'},'Box','off')
xticks(1:1:12);
saveas(f,[handles.path '/IHA_group_1.pdf']);
if getappdata(h,'canceling')
    delete(h)
    return
end
step = 25;
waitbar(step / steps)

f = figure('Name','IHA - Group 2','NumberTitle','off', 'Visible','off');
hold on
b = bar([handles.res.IHA_ind_obs_mean(13:22) handles.res.IHA_ind_syn_mean(13:22)]);
grid;
b(1).FaceColor = 'b';
b(2).FaceColor = 'r';
ylabel('Q [m^3/s]');
legend({'Q_{obs}','Q_{syn}'},'Box','off')
xp = 1:1:10;
xt = {{'Min'; '1-day'} {'Min'; '3-day'} ...
    {'Min'; '7-day'} {'Min'; '30-day'}...
    {'Min'; '90-day'} {'Max'; '1-day'}...
    {'Max'; '3-day'} {'Max'; '7-day'}...
    {'Max'; '30-day'} {'Max'; '90-day'}};
my_xticklabels(xp, xt);
saveas(f,[handles.path '/IHA_group_2.pdf']);
if getappdata(h,'canceling')
    delete(h)
    return
end
step = 26;
waitbar(step / steps)

f = figure('Name','IHA - Group 4','NumberTitle','off', 'Visible','off');
hold on
b = bar([handles.res.IHA_ind_obs_mean(27:30) handles.res.IHA_ind_syn_mean(27:30)]);
grid;
b(1).FaceColor = 'b';
b(2).FaceColor = 'r';
legend({'Q_{obs}','Q_{syn}'},'Box','off')
xp = 1:1:4;
xt = {{'No. of low pulses'} {'No. of high pulses'} {'Mean duration'; 'of low pulses'} {'Mean duration'; 'of high pulses'}};
my_xticklabels(xp, xt);
saveas(f,[handles.path '/IHA_group_4.pdf']);
if getappdata(h,'canceling')
    delete(h)
    return
end
step = 27;
waitbar(step / steps)

f = figure('Name','IHA - Group 5-1','NumberTitle','off', 'Visible','off');
hold on
b = bar([handles.res.IHA_ind_obs_mean(31:32) handles.res.IHA_ind_syn_mean(31:32)]);
grid;
b(1).FaceColor = 'b';
b(2).FaceColor = 'r';
legend({'Q_{obs}','Q_{syn}'},'Box','off')
xp = [1 2];
xt = {{'Means of all negative'; 'differences between'; 'consecutive daily values'} ...
    {'Means of all postive'; 'differences between'; 'consecutive daily means'}};
my_xticklabels(xp, xt);
saveas(f,[handles.path '/IHA_group_5-1.pdf']);
if getappdata(h,'canceling')
    delete(h)
    return
end
step = 28;
waitbar(step / steps)

f = figure('Name','IHA - Group 5-2','NumberTitle','off', 'Visible','off');
hold on
b = bar([handles.res.IHA_ind_obs_mean(33:34) handles.res.IHA_ind_syn_mean(33:34)]);
grid;
b(1).FaceColor = 'b';
b(2).FaceColor = 'r';
legend({'Q_{obs}','Q_{syn}'},'Box','off','Location','north')
xticks([1 2])
xticklabels({'No. of falls', 'No. of rises'})
saveas(f,[handles.path '/IHA_group_5-2.pdf']);
if getappdata(h,'canceling')
    delete(h)
    return
end
step = 29;
waitbar(step / steps)

f = figure('Name','Simple test statistic - Min, Max, Mean, Std & Skew','NumberTitle','off', 'Visible','off');
b = bar([handles.res.test_stats(:,1) handles.res.test_stats(:,2)]);
b(1).FaceColor = 'b';
b(2).FaceColor = 'r';
grid;
legend({'Q_{obs}','Q_{syn}'},'Box','off','Location','northwest')
ylabel('Q [m^3/s]');
xticks([1 2 3 4 5])
xticklabels({'Min', 'Max', 'Mean', 'Std', 'Skew'})
saveas(f, [handles.path '/min_max_mean_std_kew.pdf']);
if getappdata(h,'canceling')
    delete(h)
    return
end
step = 30;
waitbar(step / steps)

f = figure('Name','Volume - Cumulated Volume (annually)','NumberTitle','off', 'Visible','off');
hold on
plot(handles.res.tt_vol_yy.Time, handles.res.tt_vol_yy.vol_obs_cum, 'b');
plot(handles.res.tt_vol_yy.Time, handles.res.tt_vol_yy.vol_syn_cum, 'r');
hold off
grid;
xlabel('Date');
ylabel('Volume [km^3]');
xtickformat('yyyy');
legend({'Q_{obs}','Q_{syn}'},'Box','off','Location','northwest');
saveas(f,[handles.path '/cumulated_volume_annually.pdf']);
if getappdata(h,'canceling')
    delete(h)
    return
end
step = 31;
waitbar(step / steps)

f = figure('Name','Volume - Cumulated Volume (monthly)','NumberTitle','off', 'Visible','off');
hold on
plot(handles.res.tt_vol_mm.Time, handles.res.tt_vol_mm.vol_obs_cum, 'b');
plot(handles.res.tt_vol_mm.Time, handles.res.tt_vol_mm.vol_syn_cum, 'r');
hold off
grid;
xlabel('Date');
ylabel('Volume [km^3]');
xtickformat('MMM-yyyy');
legend({'Q_{obs}','Q_{syn}'},'Box','off','Location','northwest');
saveas(f,[handles.path '/cumulated_volume_monthly.pdf']);
if getappdata(h,'canceling')
    delete(h)
    return
end
step = 32;
waitbar(step / steps)

tt_res = timetable(handles.res.tt_syn.Date,handles.res.tt_syn.Q_sim_re);
Q_sim_res = timetable2table(tt_res);
Q_sim_res.Properties.VariableNames = {'DDMMYYYY' 'Q'};
writetable(Q_sim_res,[handles.path  '/Q_syn.csv']);
if getappdata(h,'canceling')
    delete(h)
    return
end
step = 33;
waitbar(step / steps)
delete(h)
msgbox('Done!');

function pushbuttonCancel_Callback(hObject, eventdata, handles)
delete(handles.figureExport);
