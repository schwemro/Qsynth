% Main file for generating synthetic runoff time series with Qsynth
%
% MATLAB R2017a
% (c) Copyright 2017, Robin Schwemmle <rschwemmle@yahoo.de>

% clear workspace
clear all;
close all;
clear classes;
clc;

%% generate artificial runoff time series
tic
approach=1; % 1 = 'AR' or 2 = 'ARMA'
opt_val = zeros(10,1);
Qsyn = Qsynth(approach, '01-01-2000', '31-12-2019');
Qsyn.folderres('pre');
Qsyn.importts('input/pfinz_00_06.csv', 'dd.MM.yyyy');
Qsyn.perennial();
Qsyn.tt_obs.Q = Qsyn.fillgaps(Qsyn.tt_obs.Q);
Qsyn.tt_obs.Q_trans = Qsyn.logtransform(Qsyn.tt_obs.Q);
Qsyn.testseas();
Qsyn.rmseas(Qsyn.Q_regime_daily, Qsyn.Q_std_daily);
Qsyn.determineorder();
Qsyn.selectmodel();
Q_max = max(Qsyn.tt_obs.Q)*1.1;
i=1;
for i =1:10
    prog = (i/10)*100;
    disp([num2str(prog) '%'])
    Qsyn.folderres(num2str(i));
    Q_sim1 = Qsyn.generaterunoff(Qsyn.EstMdl,Qsyn.Q_regime_daily, Qsyn.Q_std_daily);
    Qsyn.tt_syn.Q_sim_re = Q_sim1;
    Q_cp = Qsyn.cutpeaksrgm(Q_max,Qsyn.EstMdl);
    Qsyn.tt_syn.Q_sim_re = Q_cp;
    Qsyn.optMAwMWS();
    opt_val(i,1) = Qsyn.opt_val;
    Q_ma = Qsyn.MAwMWS(Q_cp,Qsyn.w,Qsyn.Qx,'lin');
    Qsyn.tt_syn.Q_sim_re = Q_ma;
    Qsyn.testnorm();
    Qsyn.testautocor();
    Qsyn.ACFmonths();
    Qsyn.compareIHA();
    Qsyn.volume();
    Qsyn.teststats();
    Qsyn.teststatstoxls();
    Qsyn.Qsimtocsv(Q_ma, 'syn_fin');
    
    f = figure('Name','Q_obs vs Q_syn','NumberTitle','off','defaultFigureVisible','off');
    hold on
    plot(Qsyn.tt_obs.Date, Qsyn.tt_obs.Q, 'b');
    plot(Qsyn.tt_syn.Date, Q_ma, 'r');
    hold off
%     h1 = [Qsyn.tt_syn.Date(1) Qsyn.tt_syn.Date(end)];
%     h2 = [Q_max Q_max];
%     line(h1,h2,'Color','black','LineStyle','-.')
    grid;
    xlabel('Date');
    ylabel('Q [m^3/s]');
    xlim(datetime(2000,[1 12],[1 31]));
    xtickformat('dd-MMM-yyyy');
    legend({'Q_{obs}','Q_{syn}'},'Box','off');
    saveas(f,[Qsyn.dir_results '/Q_obs_vs_Q_syn.fig']);
    close all;
end

[val, idx] = min(opt_val);
disp(['Folder ' num2str(idx) ' shows the best results!'])
dlmwrite('opt_val.txt',opt_val,'delimiter','\t','precision',3);
t = round(toc/60,1);
disp(['Total runtime: ' num2str(t) ' min'])
