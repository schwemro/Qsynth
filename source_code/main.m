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
Qsyn.folderres('/Users/robinschwemmle/Desktop/MATLAB/Stream_Scale_Model_for_SHP_Efficiency/Qsynth/results/inverliever/AR');
Qsyn.importts('/Users/robinschwemmle/Desktop/MATLAB/Stream_Scale_Model_for_SHP_Efficiency/Qsynth/input/inverliever_11_16.csv', 'dd/MM/yyyy', 'N/A');
Qsyn.perennial();
Qsyn.fillgaps(Qsyn.tt_obs.Q);
Qsyn.tt_obs.Q_trans = Qsyn.logtransform(Qsyn.tt_obs.Q);
Qsyn.testseas();
Qsyn.rmseas(Qsyn.Q_regime_daily, Qsyn.Q_std_daily);
Qsyn.determineorder();
Qsyn.selectmodel(Qsyn.tt_obs.Q_trans_stand_d);
close all;
Q_max = max(Qsyn.tt_obs.Q)*1.1;
clusters = {1,l}; % building clusters which are necessary for the parallelization
for i = 1:l
    clusters{1,i} = Qsynth();
    clusters{1,i}.folderres(['/Users/robinschwemmle/Desktop/MATLAB/Stream_Scale_Model_for_SHP_Efficiency/Qsynth/results/inverliever/AR/' num2str(i)]);
    clusters{1,i}.selectapproach(approach);
    clusters{1,i}.settimeperiod('01/01/2011', '31/12/2030', 'dd/MM/yyyy');
    clusters{1,i}.tt_obs = Qsyn.tt_obs;
    clusters{1,i}.Q_regime_daily = Qsyn.Q_regime_daily;
    clusters{1,i}.Q_std_daily = Qsyn.Q_std_daily;
    clusters{1,i}.N_obs = Qsyn.N_obs;
    clusters{1,i}.p = Qsyn.p;
    clusters{1,i}.q = Qsyn.q;
    clusters{1,i}.EstMdl = Qsyn.EstMdl;
end
%parpool('local'); % Parallelization of the loop
hbar = parfor_progressbar(l,'Please wait...');
parfor i = 1:l
    clusters{1,i}.generaterunoff(clusters{1,i}.EstMdl,clusters{1,i}.Q_regime_daily, clusters{1,i}.Q_std_daily);
    clusters{1,i}.cutpeaks(Q_max,clusters{1,i}.EstMdl,clusters{1,i}.tt_syn);
    clusters{1,i}.optMAwMWS();
    clusters{1,i}.tt_syn = clusters{1,i}.tt_syn; % parallelization is not returning
    clusters{1,i}.N_sim = clusters{1,i}.N_sim; % the class properties automatically.
    clusters{1,i}.w = clusters{1,i}.w; % Very unhandy!
    clusters{1,i}.Qx = clusters{1,i}.Qx;
    clusters{1,i}.w = clusters{1,i}.w;
    opt_val(i,1) = clusters{1,i}.opt_val;
    hbar.iterate(1);
end

close(hbar)

[val, idx] = min(opt_val);
clusters{1,idx}.dir_results = '/Users/robinschwemmle/Desktop/MATLAB/Stream_Scale_Model_for_SHP_Efficiency/Qsynth/results/inverliever/AR/';
clusters{1,idx}.tt_syn.Q_sim_re = clusters{1,idx}.MAwMWS(clusters{1,idx}.tt_syn.Q_sim_re,clusters{1,idx}.w,clusters{1,idx}.Qx,'lin');
clusters{1,idx}.testnorm();
clusters{1,idx}.testautocor();
clusters{1,idx}.ACFmonths();
clusters{1,idx}.compareIHA();
clusters{1,idx}.volume();
clusters{1,idx}.plotts();
clusters{1,idx}.teststats();
clusters{1,idx}.teststatstotxt();
clusters{1,idx}.Qtocsv(clusters{1,idx}.tt_syn.Date, clusters{1,idx}.tt_syn.Q_sim_re, 'syn_fin');
dlmwrite('/Users/robinschwemmle/Desktop/MATLAB/Stream_Scale_Model_for_SHP_Efficiency/Qsynth/results/inverliever/AR/opt_val.txt',opt_val,'delimiter','\t','precision',3);

t = round(toc/60,1);
disp(['Total runtime: ' num2str(t) ' min'])
% save('inverliever_AR.mat');
p = gcp;
delete(p)
close all;
