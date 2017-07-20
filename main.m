% Main file for generating synthetic runoff time series with Qsynth
%
% MATLAB R2017a
% (c) Copyright 2017, Robin Schwemmle <rschwemmle@yahoo.de>

tic
l = 10;
approach=1; % 1 = 'AR' or 2 = 'ARMA'
opt_val = zeros(10,1);
Qsyn = Qsynth();
Qsyn.selectapproach(approach);
Qsyn.settimeperiod('01/01/2011', '31/12/2030', 'dd/MM/yyyy');
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
for h = 1:l
    clusters{1,h} = Qsynth();
    clusters{1,h}.selectapproach(approach);
    clusters{1,h}.settimeperiod('01/01/2011', '31/12/2030', 'dd/MM/yyyy');
    clusters{1,h}.folderres(num2str(h));
    clusters{1,h}.tt_obs = Qsyn.tt_obs;
    clusters{1,h}.Q_regime_daily = Qsyn.Q_regime_daily;
    clusters{1,h}.Q_std_daily = Qsyn.Q_std_daily;
    clusters{1,h}.N_obs = Qsyn.N_obs;
    clusters{1,h}.p = Qsyn.p;
    clusters{1,h}.q = Qsyn.q;
    clusters{1,h}.EstMdl = Qsyn.EstMdl;
end
parpool('local'); % Parallelization of the loop
% ppm = ParforProgMon('Progress', l, 1); % display progress, not working yet
parfor i = 1:numel(clusters)
    clusters{1,i}.generaterunoff(clusters{1,i}.EstMdl,clusters{1,i}.Q_regime_daily, clusters{1,i}.Q_std_daily);
    clusters{1,i}.cutpeaks(Q_max,clusters{1,i}.EstMdl,clusters{1,i}.tt_syn);
    clusters{1,i}.optMAwMWS();
    opt_val(i,1) = clusters{1,i}.opt_val;
    clusters{1,i}.MAwMWS(Q_cp,clusters{1,i}.w,clusters{1,i}.Qx,'lin');
    clusters{1,i}.testnorm();
    clusters{1,i}.testautocor();
    clusters{1,i}.ACFmonths();
    clusters{1,i}.compareIHA();
    clusters{1,i}.volume();
    clusters{1,i}.plotts();
    clusters{1,i}.teststats();
    clusters{1,i}.teststatstotxt();
    clusters{1,i}.Qtocsv(clusters{1,i}.tt_syn.Date, clusters{1,i}.tt_syn.Q_sim_re, 'syn_fin');
%     ppm.increment();
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
