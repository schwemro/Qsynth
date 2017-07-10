% Testing robustness of runoff generator Qsynth
%
% MATLAB R2017a
% (c) Copyright 2017, Robin Schwemmle <rschwemmle@yahoo.de>
tic
l = 10;
set(0,'defaultFigureVisible','off')
approach=1; % 1 = 'AR' or 2 = 'ARMA'
start_date = datetime('01-01-1995','InputFormat','dd-MM-yyyy');
end_date = datetime('31-12-2015','InputFormat','dd-MM-yyyy');
Date = [start_date:end_date]';
acf_mc = zeros(101,l);
Q_syn_mc = zeros(length(Date),l);
IHA_mc = zeros(l,34);
test_stats_mc = zeros(l,5);
Qsyn = Qsynth(approach, '01-01-1995', '31-12-2015');
Qsyn.folderres('pre');
Qsyn.importts('input/89006_95_02.csv', 'dd/MM/yyyy');
Qsyn.perennial();
Qsyn.tt_obs.Q = Qsyn.fillgaps(Qsyn.tt_obs.Q);
Qsyn.tt_obs.Q_trans = Qsyn.logtransform(Qsyn.tt_obs.Q);
Qsyn.testseas();
Qsyn.rmseas(Qsyn.Q_regime_daily, Qsyn.Q_std_daily);
Qsyn.determineorder();
Qsyn.selectmodel();
clusters = {};
close all;
for i = 1:l
   clusters{1,i} = Qsyn;
end
clusters_perc = {};
for i = 1:l
   clusters_perc{1,i};
end
parpool('local');
ppm = ParforProgMon('Progress', l);
parfor i = 1:numel(clusters)
    clusters{1,i}.folderres(num2str(i));
    Q_sim1 = clusters{1,i}.generaterunoff(clusters{1,i}.EstMdl,clusters{1,i}.Q_regime_daily, clusters{1,i}.Q_std_daily);
    clusters{1,i}.tt_syn.Q_sim_re = Q_sim1;
    Q_max = max(clusters{1,i}.tt_obs.Q)*1.1;
    Q_cp = clusters{1,i}.cutpeaksrgm(Q_max,clusters{1,i}.EstMdl);
    clusters{1,i}.tt_syn.Q_sim_re = Q_cp;
    clusters{1,i}.optMAwMWS();
    Q_ma = clusters{1,i}.MAwMWS(Q_cp,clusters{1,i}.w,clusters{1,i}.Qx,'lin');
    clusters{1,i}.tt_syn.Q_sim_re = Q_ma;
    clusters{1,i}.teststats();
    test_stats_mc(i,:) = clusters{1,i}.test_stats(:,2);
    Q_syn_mc(:,i) = clusters{1,i}.tt_syn.Q_sim_re;
    acf_mc(:,i) = autocorr(clusters{1,i}.tt_syn.Q_sim_re, 100);
    yy = clusters{1,i}.tt_syn.YYYY(1);
    mm = clusters{1,i}.tt_syn.MM(1);
    dd = clusters{1,i}.tt_syn.dd(1);
    date = [ dd mm yy ];
    perc.q25th = prctile(clusters{1,i}.tt_syn.Q_sim_re,25);
    perc.q75th = prctile(clusters{1,i}.tt_syn.Q_sim_re,75);
    init_date = date;
    init_year = yy;
    [IHA_ind_syn] = IHA_indicators( clusters{1,i}.tt_syn.Q_sim_re, perc, init_date, init_year );
    IHA_mc(i,:) = mean(IHA_ind_syn,2);
    rmdir(num2str(i));
    ppm.increment();
end

acf_lq = prctile(acf_mc,2.5,2);
acf_uq = prctile(acf_mc,97.5,2);
acf_med = prctile(acf_mc,50,2);
acf_obs = autocorr(Qsyn.tt_obs.Q, 100);

Q_syn_lq = prctile(Q_syn_mc,2.5,2);
Q_syn_uq = prctile(Q_syn_mc,97.5,2);
Q_syn_med = prctile(Q_syn_mc,50,2);

save('MC.mat');
mkdir('MC');

f1 = figure('Name','Confidence interval ACF','NumberTitle','off');
x = [0:100]';
xx = [x;flipud(x)];
yy = [acf_lq;flipud(acf_uq)];
fill(xx,yy,[0.7 0.7 0.7],'EdgeColor','None')
hold on
plot(x,acf_med,'k','LineWidth',1.5);
plot(x,acf_obs,'b','LineWidth',1.5);
hold off
xlim([0 100]);
xlabel('Lag [Days]');
ylabel('ACF');
legend({'Confidence interval 95%','Median','Q_{obs}'},'Box','off','Location','northeast');
saveas(f1, 'MC/confidence_interval_ACF.fig');

f2 = figure('Name','Confidence interval Q_syn','NumberTitle','off');
xx = [Date;flipud(Date)];
yy = [Q_syn_lq;flipud(Q_syn_uq)];
fill(xx,yy,[0.7 0.7 0.7],'EdgeColor','None')
hold on
plot(Date,Q_syn_med,'k','LineWidth',1.5);
plot(Qsyn.tt_obs.Date,Qsyn.tt_obs.Q,'b','LineWidth',1.5);
hold off
xlabel('Date');
ylabel('Q [m^3/s]');
legend({'Confidence interval 95%','Median','Q_{obs}'},'Box','off','Location','northeast');
saveas(f2, 'MC/confidence_interval_Q_syn.fig');

f3 = figure('Name','Test statistic','NumberTitle','off');
boxplot(test_stats_mc, 'Colors', 'k','Symbol', 'k+')
ylabel('Q [m^3/s]');
xticks([1 2 3 4 5])
xticklabels({'Mean', 'Std', 'Skew', 'Min', 'Max'})
saveas(f3, 'MC/boxplot_test_stats.fig');

f4 = figure('Name','IHA - Group 1','NumberTitle','off');
boxplot(IHA_mc(:,1:12), 'Colors', 'k','Symbol', 'k+')
ylabel('Q [m^3/s]');
saveas(f4, 'MC/boxplot_IHA_group_1.fig');

f5 = figure('Name','IHA - Group 2','NumberTitle','off');
boxplot(IHA_mc(:,13:22), 'Colors', 'k','Symbol', 'k+')
ylabel('Q [m^3/s]');
xp = 1:1:10;
xt = {{'Min'; '1-day'} {'Min'; '3-day'} ...
    {'Min'; '7-day'} {'Min'; '30-day'}...
    {'Min'; '90-day'} {'Max'; '1-day'}...
    {'Max'; '3-day'} {'Max'; '7-day'}...
    {'Max'; '30-day'} {'Max'; '90-day'}};
ht = my_xticklabels(gca, xp, xt);
saveas(f5, 'MC/boxplot_IHA_group_2.fig');

f6 = figure('Name','IHA - Group 4','NumberTitle','off');
boxplot(IHA_mc(:,27:30), 'Colors', 'k','Symbol', 'k+')
xp = 1:1:4;
xt = {{'No. of low pulses'} {'No. of high pulses'} {'Mean duration'; 'of low pulses'} {'Mean duration'; 'of high pulses'}};
ht = my_xticklabels(gca, xp, xt);
saveas(f6, 'MC/boxplot_IHA_group_4.fig');

f7 = figure('Name','IHA - Group 5-1','NumberTitle','off');
boxplot(IHA_mc(:,31:32), 'Colors', 'k','Symbol', 'k+')
xp = [1 2];
xt = {{'Means of all negative'; 'differences between'; 'consecutive daily values'} ...
    {'Means of all postive'; 'differences between'; 'consecutive daily means'}};
ht = my_xticklabels(gca, xp, xt);
saveas(f7, 'MC/boxplot_IHA_group_5_1.fig');

f8 = figure('Name','IHA - Group 5-2','NumberTitle','off');
boxplot(IHA_mc(:,33:34), 'Colors', 'k','Symbol', 'k+')
xticks([1 2])
xticklabels({'No. of falls', 'No. of rises'})
saveas(f8, 'MC/boxplot_IHA_group_5_2.fig');

t = round(toc/60,1);
disp(['Total runtime: ' num2str(t) ' min'])
