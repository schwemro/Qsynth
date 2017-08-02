classdef Qsynth < handle
    % Qsynth Implementing an approach to generate artificial daily streamflow
    % time series according to Salas (1993)
    %
    % Build for GUI usage.
    % Requires: MY_XTICKLABELS
    % MATLAB R2017a,
    % (c) Copyright 2017, Robin Schwemmle <rschwemmle@yahoo.de>

    properties
        start_date % Start date of synthetic runoff time series
        end_date % End date of synthetic runoff time series
        tt_obs % Timetable with observed runoff
        tt_syn % Timetable with synthetic runoff
        Q_regime_monthly % Monthly runoff regime
        Q_regime_daily % Daily runoff regime
        Q_std_daily % Monthly runoff regime
        N_obs % Number of observations
        N_sim % Number of simulated observations
        p % AR order
        w % Optimal window size for MA with MWS
        mws % table for MA with MWS
        Qx % Upper threshold for MA with MWS
        x1_nn % lower x-axis boundary nonnormally distributed histogram
        x2_nn % upper x-axis boundary nonnormally distributed histogram
        x1_n % lower x-axis boundary normally distributed histogram
        x2_n % upper x-axis boundary normally distributed histogram
        y2_nn % upper y-axis boundary nonnormally distributed histogram
        y2_n % upper y-axis boundary lower y-axis boundary of nonnormally distributed histogram
        nbins_nn % amount bins nonnormally distributed histogram
        nbins_n % amount bins normally distributed histogram
        IHA_ind_obs_mean % IHA for observed streamflow
        IHA_ind_syn_mean % IHA for synthetic streamflow
        tt_vol_mm % table with monthly cumulated streamflow volume
        tt_vol_yy % table with yearly cumulated streamflow volume
        opt_val % Geometric mean of IHA and ACF for optimal w and Qx
        EstMdl % Time series model with estimated parameters (see MATLAB Documentation "arima class")
        test_stats % Table with simple test statistic
        dir_results % Folder where to store the results
    end

    methods
        function obj = Qsynth()
            % Initializing object.
            if nargin >= 1
                error('Too many input arguments. No input argument is necessary.')
            end
        end

        function settimeperiod(obj, start_date, end_date, dateformat_str)
            % Indicate start and end date of synthetic streamflow time
            % series. Input arguments are start and end date and the
            % corresponding dateformat.
            if nargin < 4
                error('Three input arguments are necessary.')
            end
            if nargin > 4
                error('Too many input arguments. Only three input argument are necessary.')
            end
            obj.start_date = datetime(start_date,'InputFormat',dateformat_str);
            obj.end_date = datetime(end_date,'InputFormat',dateformat_str);
        end

        function importts(obj, filepath_str, dateformat_str, na_str)
            % Import runoff time series. Input arguments are path to
            % .csv-file, format of date column and string which marks NaN.
            if nargin < 4
                error('Three input arguments are necessary.')
            end
            if nargin > 4
                error('Too many input arguments. Only three input argument are necessary.')
            end
            T = readtable(filepath_str,'TreatAsEmpty',na_str);
            T.Properties.VariableNames = {'Date' 'Q'};
            Date = datetime(T.Date,'InputFormat',dateformat_str);
            obj.tt_obs = timetable(Date,T.Q);
            obj.tt_obs.Properties.VariableNames = {'Q'};
            obj.tt_obs.Q(obj.tt_obs.Q<0) = NaN;
            obj.tt_obs.DD = day(obj.tt_obs.Date, 'dayofyear');
            obj.tt_obs.dd = day(obj.tt_obs.Date);
            obj.tt_obs.MM = month(obj.tt_obs.Date);
            obj.tt_obs.YYYY = year(obj.tt_obs.Date);
            obj.N_obs = length(obj.tt_obs.Q);
        end

        function perennial(obj)
           % Testing if river is perennial. Otherwise rasing ERROR!
           if nargin >= 2
               error('Too many input arguments. No input argument is necessary.')
           end
           if any(obj.tt_obs.Q==0)
               errordlg({'River is not perennial. Approach is not appropriate' 'and will generate wrong results!'},'Error')
               return
           end
        end

        function checkobs(obj)
            % Testing if observed streamflow time series starts to a 1
            % January and ends to a 31 December. Otherwise rasing ERROR!
            % Only entire years can be used for the streamflow generator.
            if nargin >= 2
                error('Too many input arguments. No input argument is necessary.')
            end

            if (obj.tt_obs.DD(1)~=1 && obj.tt_obs.MM(1)~=1)
                errordlg({'Time series of observed streamflow does not' ...
                    'start to a 1 January. In order to run the' 'streamflow Generator complete years are necessary.'},'Error')
                return
            end
            if (obj.tt_obs.DD(end)~=31 && obj.tt_obs.MM(end)~=12)
                errordlg({'Time series of observed streamflow does not' ...
                    'end to a 31 December.In order to run the' 'streamflow Generator complete years are necessary.'},'Error')
                return
            end
        end


        function gaps(obj)
            % Testing if any gaps with missing values are available. Rasing ERROR in case gaps are too wide!
            if nargin >= 2
                error('Too many input arguments. No input argument is necessary.')
            end
            bool_na = isnan(obj.tt_obs.Q);
            na_start = find(diff([0; bool_na; 0])==1);
            na_end = find(diff([0; bool_na; 0])==-1);
            durat_na = na_end - na_start;

            if (max(durat_na)>30)
                errordlg({'Gaps with missing values are too wide.' 'Gaps need to be filled accordingly' ...
                    'to run the streamflow generator.'},'Error')
                return
            else
                obj.fillgaps(obj.tt_obs.Q);
                obj.tt_obs.Q_trans = obj.logtransform(obj.tt_obs.Q);
                obj.testseas();
            end
        end

        function fillgaps(obj, Q)
            % Testing if time series exhibits gaps. If there are any gaps,
            % they are filled by using autoregressive modeling using half the
            % length of the gap as number of samples in the estimation.
            % Input argument is time series.
            if nargin < 2
                error('One input argument is necessary.')
            end
            if nargin > 2
                error('Too many input arguments. Only one input argument is necessary.')
            end
            Q_NA = double(isnan(Q));
            N_NA = sum(Q_NA);
            if (N_NA == 0)
                disp('No Gaps in observed runoff time series!')
            elseif (N_NA >= 1)
                NA_idx = find(Q_NA==1);
                idx = [1:1:length(NA_idx)]';
                t=table(idx,NA_idx);
                t_idx = t.idx(diff(t.NA_idx)>1);
                t_idx(end+1) = t.idx(end);
                t_idx(2:end+1)=t_idx;
                t_idx(1)=0;
                tt_obs_NA_idx = NA_idx(diff(t.NA_idx)>1);
                tt_obs_NA_idx(end+1) = t.NA_idx(end);
                x = zeros(length(tt_obs_NA_idx),3);
                x(:,2) = tt_obs_NA_idx;
                x(:,1) = tt_obs_NA_idx - (diff(t_idx)-1);
                x(:,3) = x(:,2) - x(:,1);

                Q_min = nanmin(Q);

                for i = 1:length(tt_obs_NA_idx)
                    maxlen = ceil(x(i,3)/2);
                    if (maxlen == 1)
                        y = fillgaps(Q,maxlen+2);
                        Q(x(i,1):x(i,2),1) = y(x(i,1):x(i,2),1);
                    elseif (maxlen > 1)
                        y = fillgaps(Q,maxlen+1);
                        Q(x(i,1):x(i,2),1) = y(x(i,1):x(i,2),1);
                    elseif (maxlen == 0)
                        y = fillmissing(Q,'movmedian',3);
                        Q(x(i,1),1) = y(x(i,1),1);
                    end
                end
                Q(Q<0) = Q_min;
                obj.tt_obs.Q = Q;
            end
        end

        function folderres(obj, fn)
            % Create folder where to store the results.
            if nargin < 2
                error('One input argument is necessary.')
            end
            if nargin > 2
                error('Too many input arguments. Only one input argument is necessary.')
            end
            obj.dir_results = fn;
            mkdir(fn);
        end

        function Q_log = logtransform(~, Q)
            % Logarithmic transformation of runoff time series. Input argument
            % is time series. Output argument is log transformed time series.
            if nargin < 2
                error('One input argument is necessary.')
            end
            if nargin > 2
                error('Too many input arguments. Only one input argument is necessary.')
            end
            Q_log = log(Q);
        end

        function testseas(obj)
            % Testing visually for prevailing seasonality by plotting the
            % daily and mothly runoff regime and the parde coefficient.
            % TODO: runoff regime and std regime as output args.
            if nargin >= 2
                error('Too many input arguments. No input argument is necessary.')
            end
            Q_mean = nanmean(obj.tt_obs.Q);
            obj.tt_obs.DD = day(obj.tt_obs.Date, 'dayofyear');
            tt_obs_mon = retime(obj.tt_obs, 'monthly', 'mean');
            tt_obs_mon.MM = month(tt_obs_mon.Date);
            func = @nanmean;
            obj.Q_regime_monthly = varfun(func,tt_obs_mon,'GroupingVariables','MM');
            obj.Q_regime_monthly.parde = obj.Q_regime_monthly.nanmean_Q/Q_mean;
            obj.Q_regime_daily = varfun(func,obj.tt_obs,'GroupingVariables','DD');
            func1 = @nanstd;
            obj.Q_std_daily = varfun(func1,obj.tt_obs,'GroupingVariables','DD');
        end

        function rmseas(obj, Q_regime_daily_fit, Q_std_daily_fit)
            % Removing trends and shifts of the runoff time series by
            % daily standardization. Input argument are daily runoff regime
            % and the daily standard deviation regime.
            % TODO: Adjust transfer of input args and the use inside the
            % function.
            if nargin < 3
                error('Two input arguments are necessary.')
            end
            if nargin > 3
                error('Too many input arguments. Only two input arguments are necessary.')
            end
            obj.tt_obs.Q_rm_mean = zeros(obj.N_obs,1);
            obj.tt_obs.Q_rm_std = zeros(obj.N_obs,1);

            for i = 1:obj.N_obs
                d = obj.tt_obs.DD(i);
                obj.tt_obs.Q_rm_mean(i) = Q_regime_daily_fit.nanmean_Q_trans(Q_regime_daily_fit.DD==d);
                obj.tt_obs.Q_rm_std(i) = Q_std_daily_fit.nanstd_Q_trans(Q_std_daily_fit.DD==d);
            end
            obj.tt_obs.Q_trans_stand_d = (obj.tt_obs.Q_trans-obj.tt_obs.Q_rm_mean)./obj.tt_obs.Q_rm_std;
        end

        function determineorder(obj)
            % Determine model order by using the autocorrelation
            % and the partial autocorrelation function.
            if nargin >= 2
                error('Too many input arguments. No input argument is necessary.')
            end
            [pacf,lags_pacf,bounds_pacf] = parcorr(obj.tt_obs.Q_trans_stand_d, 20);
            k = find(pacf<bounds_pacf(1)&pacf>bounds_pacf(2));
            obj.p = k(1) - 1;
        end

        function selectmodel(obj, Q_trans_stand)
            % Initialize object of time series model and estimate the according
            % model paramters. Input argument is the log-transformed and standardized streamflow time
            % series.
            if nargin < 2
                error('One input argument is necessary.')
            end
            if nargin > 2
                error('Too many input arguments. Only one input argument is necessary.')
            end
            Mdl = arima(obj.p,0,0);

            obj.EstMdl = estimate(Mdl,Q_trans_stand);
            obj.p = length(obj.EstMdl.AR);
        end

        function generaterunoff(obj, EstMdl, Q_regime_daily_fit, Q_std_daily_fit)
            % Generating articial streamflow time series. Readding the trend as well as undoing the
            % standardization and logarithmic transformation. Input arguments are the estimated model,
            % daily runoff regime and the daily standard deviation regime.
            if nargin < 4
                error('Three input arguments are necessary.')
            end
            if nargin > 4
                error('Too many input arguments. Only three input argument are necessary.')
            end
            Date = [obj.start_date:obj.end_date]';
            obj.N_sim = length(Date);
            obj.tt_syn = timetable(Date);
            [obj.tt_syn.Q_sim,obj.tt_syn.E] = simulate(EstMdl,obj.N_sim);
            obj.tt_syn.DD = day(obj.tt_syn.Date, 'dayofyear');
            obj.tt_syn.Q_rm_mean = zeros(obj.N_sim,1);
            obj.tt_syn.Q_rm_std = zeros(obj.N_sim,1);
            obj.tt_syn.dd = day(obj.tt_syn.Date);
            obj.tt_syn.MM = month(obj.tt_syn.Date);
            obj.tt_syn.YYYY = year(obj.tt_syn.Date);

            for i = 1:obj.N_sim
                d = obj.tt_syn.DD(i);
                obj.tt_syn.Q_rm_mean(i) = Q_regime_daily_fit.nanmean_Q_trans(Q_regime_daily_fit.DD==d);
                obj.tt_syn.Q_rm_std(i) = Q_std_daily_fit.nanstd_Q_trans(Q_std_daily_fit.DD==d);
            end

            obj.tt_syn.Q_sim_re = exp(obj.tt_syn.Q_sim.*obj.tt_syn.Q_rm_std+obj.tt_syn.Q_rm_mean);
        end

        function Q_sim = regeneraterunoff(~, EstMdl, tt_sim)
            % Regenerating articial streamflow time series. Input arguments are
            % the estimated model and the timetable of the synthetic time series.
            % Output argument is the simulated streamflow time series.
            if nargin < 3
                error('Two input arguments are necessary.')
            end
            if nargin > 3
                error('Too many input arguments. Only two input arguments are necessary.')
            end
            [tt_sim.Q_sim,tt_sim.E] = simulate(EstMdl,size(tt_sim,1));

            Q_sim = exp(tt_sim.Q_sim.*tt_sim.Q_rm_std+tt_sim.Q_rm_mean);
        end

        function cutpeaks(obj, Q_max, EstMdl, tt_sim)
            % Cut unlikely high peaks (Q_sim > Q_max + .1*Q_max) by regenerating
            % time series with model. Starts at 5 day backward shift. Input arguments
            % are the threshold above which the peaks will be cutted, the estimated model
            % and the timetable of the synthetic time series.
            % Output argument is the simulated streamflow time series.
            if nargin < 4
                error('Three input arguments are necessary.')
            end
            if nargin > 4
                error('Too many input arguments. Only three input argument are necessary.')
            end
            if any(obj.tt_syn.Q_sim_re>Q_max)
                tau = 5;
                Q_sim_cp_tt = tt_sim;
                while 1
                    peak = find(Q_sim_cp_tt.Q_sim_re>Q_max,1);
                    if (peak>tau)
                        temp = Q_sim_cp_tt(peak-tau:end,:);
                        Q_sim_cp_tt.Q_sim_re(peak-tau:end) = obj.regeneraterunoff(EstMdl, temp);
                    elseif (peak<=tau && peak>1)
                        temp = Q_sim_cp_tt(peak-(peak-1):end,:);
                        Q_sim_cp_tt.Q_sim_re(peak-(peak-1):end) = obj.regeneraterunoff(EstMdl, temp);
                    elseif (peak==1)
                        temp = Q_sim_cp_tt(peak:end,:);
                        Q_sim_cp_tt.Q_sim_re(peak:end) = obj.regeneraterunoff(EstMdl, temp);
                    end
                    if (Q_sim_cp_tt.Q_sim_re<Q_max)
                        obj.tt_syn.Q_sim_re = Q_sim_cp_tt.Q_sim_re;
                        obj.tt_syn.Q_sim = Q_sim_cp_tt.Q_sim;
                        obj.tt_syn.Q_sim_re = Q_sim_cp_tt.Q_sim_re;
                        break
                    end
                end
            else
                disp('No occurence of unlikely high peaks!')
            end
        end

        function optMAwMWS(obj)
            % Find optimal window size and upper threshold of Q where to start
            % applying a moving average with moving window size.
            if nargin >= 2
                error('Too many input arguments. No input argument is necessary.')
            end
            Q_obs = obj.tt_obs.Q;
            Q_syn = obj.tt_syn.Q_sim_re(obj.tt_syn.Date(1):obj.tt_syn.Date(obj.N_obs));
            yy = obj.tt_obs.YYYY(1);
            mm = obj.tt_obs.MM(1);
            dd = obj.tt_obs.dd(1);
            date =[ dd mm yy ];

            perc.q25th=prctile(Q_obs,25);
            perc.q75th=prctile(Q_obs,75);

            init_date=date;
            init_year=yy;

            [IHA_ind_obs]= IHA_indicators( Q_obs, perc, init_date, init_year );
            IHA_ind_obs_mean_0=mean(IHA_ind_obs,2);
            IHA_ind_obs_mean_1 = [IHA_ind_obs_mean_0(1:22); IHA_ind_obs_mean_0(27:34)];

            acf_obs = autocorr(Q_obs, 100);

            ws = 3:2:9;
            x2 = 10:10:60;
            res = zeros(length(ws),length(x2));

            for i = 1:length(ws)
                for j = 1:length(x2)
                    Q_syn_mws = obj.MAwMWS(Q_syn,ws(i),x2(j));
                    perc.q25th=prctile(Q_syn_mws,25);
                    perc.q75th=prctile(Q_syn_mws,75);

                    [IHA_ind_syn_mws]= IHA_indicators( Q_syn_mws, perc, init_date, init_year );
                    IHA_ind_syn_mws_mean=mean(IHA_ind_syn_mws,2);
                    IHA_ind_syn_mws_mean_1 = [IHA_ind_syn_mws_mean(1:22); IHA_ind_syn_mws_mean(27:34)];
                    IHA_diff = abs(IHA_ind_obs_mean_1-IHA_ind_syn_mws_mean_1);
                    mean_IHA_diff = mean(IHA_diff); %TODO: check if it's correct!

                    acf_syn_mws = autocorr(Q_syn_mws, 100);
                    acf_diff = abs(acf_obs-acf_syn_mws);
                    mean_acf_diff = mean(acf_diff);

                    res(i,j) = geomean([mean_IHA_diff,mean_acf_diff]);
                end
            end
            [obj.opt_val,idx] = min(res(:));
            [I_row, I_col] = ind2sub(size(res),idx);
            obj.w = ws(I_row);
            obj.Qx = x2(I_col);
        end

        function Q = MAwMWS(obj, Q, N_ws, prc)
            % Filtering time series by moving average with moving window
            % size. Input arguments are time series, the starting window size,
            % streamflow percentile as far as the filtering is applied. 
            % Output argument is a filtered time series.
            if nargin < 4
                error('Three input arguments are necessary.')
            end
            if nargin > 4
                error('Too many input arguments. Only three input argument are necessary.')
            end
            t = table();
            t.Q = Q;
            t.Q_sim_dma = zeros(length(t.Q),1);
            x2 = prctile(Q,prc);
            x1 = nanmin(Q);
            obj.mws = zeros(ceil(N_ws/2),2);

            obj.mws(:,1) = linspace(x1,x2,ceil(N_ws/2));
            obj.mws(:,2) = linspace(N_ws,1,ceil(N_ws/2));
            t.mws = zeros(length(t.Q),1);
            for i = 1:ceil(N_ws/2)
                if (i < ceil(N_ws/2))
                    t.mws(t.Q>=obj.mws(i,1)&t.Q<obj.mws(i+1,1)) = obj.mws(i,2);
                elseif (i == ceil(N_ws/2))
                    t.mws(t.Q>=obj.mws(i,1)) = obj.mws(i,2);
                end
            end

            if (N_ws == 3)
                t.ma_3 = smooth(t.Q, 3, 'moving');

                t.Q_sim_dma = t.Q;

                t.Q_sim_dma(t.mws==3) = t.ma_3(t.mws==3);

            elseif (N_ws == 5)
                t.ma_3 = smooth(t.Q, 3, 'moving');
                t.ma_5 = smooth(t.Q, 5, 'moving');

                t.Q_sim_dma = t.Q;

                t.Q_sim_dma(t.mws==3) = t.ma_3(t.mws==3);
                t.Q_sim_dma(t.mws==5) = t.ma_5(t.mws==5);

            elseif (N_ws == 7)
                t.ma_3 = smooth(t.Q, 3, 'moving');
                t.ma_5 = smooth(t.Q, 5, 'moving');
                t.ma_7 = smooth(t.Q, 7, 'moving');

                t.Q_sim_dma = t.Q;

                t.Q_sim_dma(t.mws==3) = t.ma_3(t.mws==3);
                t.Q_sim_dma(t.mws==5) = t.ma_5(t.mws==5);
                t.Q_sim_dma(t.mws==7) = t.ma_7(t.mws==7);

            elseif (N_ws == 9)
                t.ma_3 = smooth(t.Q, 3, 'moving');
                t.ma_5 = smooth(t.Q, 5, 'moving');
                t.ma_7 = smooth(t.Q, 7, 'moving');
                t.ma_9 = smooth(t.Q, 9, 'moving');

                t.Q_sim_dma = t.Q;

                t.Q_sim_dma(t.mws==3) = t.ma_3(t.mws==3);
                t.Q_sim_dma(t.mws==5) = t.ma_5(t.mws==5);
                t.Q_sim_dma(t.mws==7) = t.ma_7(t.mws==7);
                t.Q_sim_dma(t.mws==9) = t.ma_9(t.mws==9);
            end
            Q = t.Q_sim_dma;
        end

        function testnorm(obj)
            % Testing visually for normal distribution by using histograms.
            if nargin >= 2
                error('Too many input arguments. No input argument is necessary.')
            end
            xlims_nn = zeros(2,2);
            xlims_nn(1,1) = min(obj.tt_obs.Q);
            xlims_nn(2,1) = min(obj.tt_syn.Q_sim_re);
            xlims_nn(1,2) = max(obj.tt_obs.Q);
            xlims_nn(2,2) = max(obj.tt_syn.Q_sim_re);

            obj.x1_nn = min(xlims_nn(:,1));
            obj.x2_nn = max(xlims_nn(:,2));

            xlims_n = zeros(3,2);
            xlims_n(1,1) = min(obj.tt_obs.Q_trans);
            xlims_n(2,1) = min(obj.tt_obs.Q_trans_stand_d);
            xlims_n(3,1) = min(obj.tt_syn.Q_sim);
            xlims_n(1,2) = max(obj.tt_obs.Q_trans);
            xlims_n(2,2) = max(obj.tt_obs.Q_trans_stand_d);
            xlims_n(3,2) = max(obj.tt_syn.Q_sim);

            obj.x1_n = min(xlims_n(:,1));
            obj.x2_n = max(xlims_n(:,2));

            h1 = histogram(obj.tt_obs.Q,'Normalization','probability','FaceColor','b');
            h1v = h1.Values;
            h1nb = h1.NumBins;
            h2 = histogram(obj.tt_obs.Q_trans,'Normalization','probability','FaceColor','b');
            h2v = h2.Values;
            h2nb = h2.NumBins;
            h3 = histogram(obj.tt_obs.Q_trans_stand_d,'Normalization','probability','FaceColor','b');
            h3v = h3.Values;
            h3nb = h3.NumBins;
            h4 = histogram(obj.tt_syn.Q_sim,'Normalization','probability','FaceColor','r');
            h4v = h4.Values;
            h4nb = h4.NumBins;
            h5 = histogram(obj.tt_syn.Q_sim_re,'Normalization','probability','FaceColor','r');
            h5v = h5.Values;
            h5nb = h5.NumBins;

            ylims_nn = zeros(2,1);
            ylims_nn(1,1) = max(h1v);
            ylims_nn(2,1) = max(h5v);
            obj.y2_nn = round(max(ylims_nn(:,1)),1)+.1;

            bins_nn = zeros(2,1);
            bins_nn(1,1) = h1nb;
            bins_nn(2,1) = h5nb;
            obj.nbins_nn = max(bins_nn(:,1));

            ylims_n = zeros(3,1);
            ylims_n(1,1) = max(h2v);
            ylims_n(2,1) = max(h3v);
            ylims_n(3,1) = max(h4v);
            obj.y2_n = round(max(ylims_n(:,1)),1)+.02;

            bins_n = zeros(3,1);
            bins_n(1,1) = h2nb;
            bins_n(2,1) = h3nb;
            bins_n(3,1) = h4nb;
            obj.nbins_n = max(bins_n(:,1));
        end

        function compareIHA(obj)
            % Compare the IHA indicators (Richter et al. 1996) for the
            % observed and synthetic streamflow except Group 3.
            if nargin >= 2
                error('Too many input arguments. No input argument is necessary.')
            end
            year = obj.tt_obs.YYYY(1);
            date =[ 01 01 year ];

            perc.q25th=prctile(obj.tt_obs.Q,25);
            perc.q75th=prctile(obj.tt_obs.Q,75);

            init_date=date;
            init_year=year;

            [IHA_ind_obs]= IHA_indicators(obj.tt_obs.Q, perc, init_date, init_year);
            obj.IHA_ind_obs_mean=mean(IHA_ind_obs,2);

            perc.q25th=prctile(obj.tt_syn.Q_sim_re,25);
            perc.q75th=prctile(obj.tt_syn.Q_sim_re,75);

            init_date=date;
            init_year=year;

            [IHA_ind_syn]= IHA_indicators(obj.tt_syn.Q_sim_re, perc, init_date, init_year);
            obj.IHA_ind_syn_mean=mean(IHA_ind_syn,2);
        end

        function volume(obj)
           % Plotting the cumulated volumes of observed and synthetic streamflow
           % at a monthly and a yearly scale and the total volumes
           % for each month.
           if nargin >= 2
               error('Too many input arguments. No input argument is necessary.')
           end
           Q_obs = obj.tt_obs.Q;
           Q_syn = obj.tt_syn.Q_sim_re(obj.tt_syn.Date(1):obj.tt_syn.Date(obj.N_obs));

           tt_vol = timetable(obj.tt_obs.Date,Q_obs, Q_syn);
           tt_vol.Properties.VariableNames = {'Q_obs' 'Q_syn'};
           tt_vol.vol_obs = (tt_vol.Q_obs.*86400)./10^9;
           tt_vol.vol_syn = (tt_vol.Q_syn.*86400)./10^9;

           obj.tt_vol_mm = retime(tt_vol, 'monthly', 'sum');
           obj.tt_vol_mm.vol_obs_cum = cumsum(obj.tt_vol_mm.vol_obs);
           obj.tt_vol_mm.vol_syn_cum = cumsum(obj.tt_vol_mm.vol_syn);
           obj.tt_vol_mm.MM = month(obj.tt_vol_mm.Time);
           obj.tt_vol_yy = retime(tt_vol, 'yearly', 'sum');
           obj.tt_vol_yy.vol_obs_cum = cumsum(obj.tt_vol_yy.vol_obs);
           obj.tt_vol_yy.vol_syn_cum = cumsum(obj.tt_vol_yy.vol_syn);
        end

        function teststats(obj)
            % Calculate test staistsic containing mean, standard deviation,
            % skewness coeeficient, minimum and maximum of the observed and
            % the synthetic streamflow respectively.
            if nargin >= 2
                error('Too many input arguments. No input argument is necessary.')
            end
            obj.test_stats = zeros(5,2);
            obj.test_stats(1,1) = min(obj.tt_obs.Q);
            obj.test_stats(1,2) = min(obj.tt_syn.Q_sim_re);
            obj.test_stats(2,1) = max(obj.tt_obs.Q);
            obj.test_stats(2,2) = max(obj.tt_syn.Q_sim_re);
            obj.test_stats(3,1) = mean(obj.tt_obs.Q);
            obj.test_stats(3,2) = mean(obj.tt_syn.Q_sim_re);
            obj.test_stats(4,1) = std(obj.tt_obs.Q);
            obj.test_stats(4,2) = std(obj.tt_syn.Q_sim_re);
            obj.test_stats(5,1) = skewness(obj.tt_obs.Q);
            obj.test_stats(5,2) = skewness(obj.tt_syn.Q_sim_re);
        end

        function teststatstotxt(obj)
            % Export test statistic to .txt..
            if nargin >= 2
                error('Too many input arguments. No input argument is necessary.')
            end
            B = array2table(obj.test_stats);
            B.Properties.VariableNames = {'Obs' 'Sim'};
            B.Properties.RowNames = {'min' 'max' 'mean' 'std' 'skew'};
            writetable(B,[obj.dir_results '/test_stats.txt'],'WriteRowNames',true,'Delimiter','\t');
        end

        function Qtocsv(obj, Date, Q, fn)
            % Export streamflow time series to .csv.. Input arguments
            % are Date, streamflow and suffix of file name.
            if nargin < 4
                error('Three input arguments are necessary.')
            end
            if nargin > 4
                error('Too many input arguments. Only three input argument are necessary.')
            end
            tt_res = timetable(Date,Q);
            Q_sim_res = timetable2table(tt_res);
            Q_sim_res.Properties.VariableNames = {'DDMMYYYY' 'Q'};
            writetable(Q_sim_res,[obj.dir_results  '/Q_' fn '.csv']);
        end
    end
end
