% Main file for generating synthetic runoff time series with Qsynth
%
% MATLAB R2017a
% (c) Copyright 2017, Robin Schwemmle <rschwemmle@yahoo.de>

tic
l = 10;
approach=1; % 1 = 'AR' or 2 = 'ARMA'
opt_val = zeros(10,1);
Qsyn = Qsynth(approach, '01-01-2011', '31-12-2030');
Qsyn.folderres('pre');
Qsyn.importts('input/inverliever_11_16.csv', 'dd/MM/yyyy', 'N/A');
Qsyn.perennial();
Qsyn.tt_obs.Q = Qsyn.fillgaps(Qsyn.tt_obs.Q);
Qsyn.tt_obs.Q_trans = Qsyn.logtransform(Qsyn.tt_obs.Q);
Qsyn.testseas();
Qsyn.rmseas(Qsyn.Q_regime_daily, Qsyn.Q_std_daily);
Qsyn.determineorder();
Qsyn.selectmodel(Qsyn.tt_obs.Q_trans_stand_d);
close all;
Q_max = max(Qsyn.tt_obs.Q)*1.1;
clusters = {1,l}; % building clusters which are necessary for the parallelization
for i = 1:l
   clusters{1,i} = Qsyn;
end
parpool('local'); % Parallelization of the loop
% ppm = ParforProgMon('Progress', l, 1); % display progress, not working yet
parfor i = 1:numel(clusters)
    clusters{1,i}.folderres(num2str(i));
    Q_sim1 = clusters{1,i}.generaterunoff(clusters{1,i}.EstMdl,clusters{1,i}.Q_regime_daily, clusters{1,i}.Q_std_daily);
    clusters{1,i}.tt_syn.Q_sim_re = Q_sim1;
    Q_cp = clusters{1,i}.cutpeaks(Q_max,clusters{1,i}.EstMdl,clusters{1,i}.tt_syn);
    clusters{1,i}.tt_syn.Q_sim_re = Q_cp;
    clusters{1,i}.optMAwMWS();
    opt_val(i,1) = clusters{1,i}.opt_val;
    Q_ma = clusters{1,i}.MAwMWS(Q_cp,clusters{1,i}.w,clusters{1,i}.Qx,'lin');
    clusters{1,i}.tt_syn.Q_sim_re = Q_ma;
    clusters{1,i}.testnorm();
    clusters{1,i}.testautocor();
    clusters{1,i}.ACFmonths();
    clusters{1,i}.compareIHA();
    clusters{1,i}.volume();
    clusters{1,i}.teststats();
    clusters{1,i}.teststatstotxt();
    clusters{1,i}.Qtocsv(clusters{1,i}.tt_syn.Date, Q_ma, 'syn_fin');
%     ppm.increment();

    f1 = figure('Name','Q_obs vs Q_syn','NumberTitle','off','defaultFigureVisible','off');
    hold on
    plot(clusters{1,i}.tt_obs.Date, clusters{1,i}.tt_obs.Q, 'b');
    plot(clusters{1,i}.tt_syn.Date, Q_ma, 'r');
    hold off
    grid;
    xlabel('Date');
    ylabel('Q [m^3/s]');
    xlim(datetime(clusters{1,i}.tt_obs.YYYY(1),[1 12],[1 31]));
    xtickformat('dd-MMM-yyyy');
    legend({'Q_{obs}','Q_{syn}'},'Box','off');
    saveas(f1,[clusters{1,i}.dir_results '/Q_obs_vs_Q_syn_' num2str(clusters{1,i}.tt_obs.YYYY(1)) '.fig']);
    saveas(f1,[clusters{1,i}.dir_results '/Q_obs_vs_Q_syn_' num2str(clusters{1,i}.tt_obs.YYYY(1)) '.pdf']);

    f2 = figure('Name','Q_obs vs Q_syn','NumberTitle','off','defaultFigureVisible','off');
    hold on
    plot(clusters{1,i}.tt_obs.Date, clusters{1,i}.tt_obs.Q, 'b');
    plot(clusters{1,i}.tt_syn.Date, Q_ma, 'r');
    hold off
    grid;
    xlabel('Date');
    ylabel('Q [m^3/s]');
    xtickformat('dd-MMM-yyyy');
    legend({'Q_{obs}','Q_{syn}'},'Box','off');
    saveas(f2,[clusters{1,i}.dir_results '/Q_obs_vs_Q_syn.fig']);
    saveas(f2,[clusters{1,i}.dir_results '/Q_obs_vs_Q_syn.pdf']);
    close all;
end

p = gcp;
delete(p)

[val, idx] = min(opt_val);
disp(['Folder ' num2str(idx) ' shows the best results!'])
dlmwrite('opt_val.txt',opt_val,'delimiter','\t','precision',3);

t = round(toc/60,1);
disp(['Total runtime: ' num2str(t) ' min'])
% save('inverliever_AR.mat');
