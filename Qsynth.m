classdef Qsynth < handle
    % Qsynth Implementing an approach to generate artificial daily streamflow
    % time series according to Salas (1993)
    %
    % Build for command line usage
    % Requires: MY_XTICKLABELS
    % MATLAB R2017a,
    % (c) Copyright 2017, Robin Schwemmle <rschwemmle@yahoo.de>

    properties
        approach % 1 = 'AR' or 2 = 'ARMA'
        start_date % Start date of synthetic runoff time series
        end_date % End date of synthetic runoff time series
        tt_obs % Timetable with observed runoff
        tt_syn % Timetable with synthetic runoff
        Q_regime_daily % Daily runoff regime
        Q_std_daily % Monthly runoff regime
        N_obs % Number of observations
        N_sim % Number of simulated observations
        p % AR order
        q % MA order
        w % Optimal window size for MA with MWS
        Qx % Upper threshold for MA with MWS
        opt_val % Geometric mean of IHA and ACF for optimal w and Qx
        EstMdl % Time series model with estimated parameters (see MATLAB Documentation "arima class")
        test_stats % Table with simple test statistic
        dir_results % Folder where to store the results
    end

    properties (Access = protected)
        approach_str = {'AR' 'ARMA'};
    end

    methods
        function obj = Qsynth()
            % Initializing object.
            if nargin >= 1
              error('Too many input arguments. No input argument is necessary.')
            end
        end

        function selectapproach(obj, approach)
            % Selecting the approach either AR or ARMA. Input arguments is
            % approach which is applied.
            if nargin < 2
                error('One input argument is necessary.')
            end
            if nargin > 2
                error('Too many input arguments. Only one input argument is necessary.')
            end
            obj.approach = approach;
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

            f1 = figure('Name','Daily Streamflow','NumberTitle','off');
            plot(obj.tt_obs.Date, obj.tt_obs.Q,'blue');
            grid;
            xlabel('Date');
            ylabel('Q [m^3/s]');
            saveas(f1,[obj.dir_results '/Daily_Runoff.fig']);
        end

        function perennial(obj)
           % Testing if river is perennial. Otherwise rasing ERROR!
           if nargin >= 2
             error('Too many input arguments. No input argument is necessary.')
           end
           if any(obj.tt_obs.Q==0)
               disp('River is not perennial. Approach is not appropriate!')
           end
        end

        function fillgaps(~, Q)
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
                Q_ts_filled = Q;
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
            if nargin < 1
                error('One input argument is necessary.')
            end
            if nargin > 1
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
            if nargin >= 2
              error('Too many input arguments. No input argument is necessary.')
            end
            Q_mean = nanmean(obj.tt_obs.Q);
            obj.tt_obs.DD = day(obj.tt_obs.Date, 'dayofyear');
            tt_obs_mon = retime(obj.tt_obs, 'monthly', 'mean');
            tt_obs_mon.MM = month(tt_obs_mon.Date);
            func = @nanmean;
            Q_regime_monthly = varfun(func,tt_obs_mon,'GroupingVariables','MM');
            Q_regime_monthly.parde = Q_regime_monthly.nanmean_Q/Q_mean;
            obj.Q_regime_daily = varfun(func,obj.tt_obs,'GroupingVariables','DD');
            func1 = @nanstd;
            obj.Q_std_daily = varfun(func1,obj.tt_obs,'GroupingVariables','DD');

            f1 = figure('Name','Runoff Regime (monthly)','NumberTitle','off');
            plot(Q_regime_monthly.MM, Q_regime_monthly.nanmean_Q, '-bo');
            grid;
            mm = nanmean(Q_regime_monthly.nanmean_Q);
            hline = refline([0 mm]);
            hline.Color = 'blue';
            hline.LineStyle = '-.';
            xlim([1 12]);
            ylim_up = ceil(max(Q_regime_monthly.nanmean_Q));
            ylim([0,ylim_up]);
            xlabel('Month');
            ylabel('Q [m^3/s]');
            saveas(f1,[obj.dir_results '/Runoff_Regime_monthly.fig']);

            f2 = figure('Name','Parde Coefficient','NumberTitle','off');
            plot(Q_regime_monthly.MM, Q_regime_monthly.parde, '-bo');
            grid;
            hline = refline([0 1]);
            hline.Color = 'blue';
            hline.LineStyle = '-.';
            xlim([1 12]);
            ylim_up = ceil(max(Q_regime_monthly.parde));
            ylim([0,ylim_up]);
            xlabel('Month');
            ylabel('[-]');
            saveas(f2,[obj.dir_results '/Parde_coefficient.fig']);

            f3 = figure('Name','Runoff Regime (daily)','NumberTitle','off');
            plot(obj.Q_regime_daily.DD, obj.Q_regime_daily.nanmean_Q, '-b');
            grid;
            xlim([1 365]);
            md = nanmean(obj.Q_regime_daily.nanmean_Q);
            hline = refline([0 md]);
            hline.Color = 'blue';
            hline.LineStyle = '-.';
            ylim_up = ceil(max(obj.Q_regime_daily.nanmean_Q))+1;
            ylim([0,ylim_up]);
            xlabel('Day of Year');
            ylabel('Q [m^3/s]');
            saveas(f3,[obj.dir_results '/Runoff_Regime_daily.fig']);
        end

        function rmseas(obj, Q_regime_daily_fit, Q_std_daily_fit)
            % Removing trends and shifts of the runoff time series by
            % daily standardization. Input argument are daily runoff regime
            % and the daily standard deviation regime.
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

            f1 = figure('Name','Standardized Daily Runoff','NumberTitle','off');
            plot(obj.tt_obs.Date, obj.tt_obs.Q_trans_stand_d, 'blue');
            grid;
            xlabel('Date');
            ylabel('Q [m^3/s]');
            saveas(f1,[obj.dir_results '/Standardized_Daily_Runoff.fig']);
        end

        function determineorder(obj)
            % Determine model order by using the autocorrelation
            % and the partial autocorrelation function.
            if nargin >= 2
              error('Too many input arguments. No input argument is necessary.')
            end
            if (obj.approach == 1)
                obj.q = 0;
            elseif (obj.approach == 2)
                acf = autocorr(obj.tt_obs.Q_trans_stand_d, obj.N_obs-1);
                bounds_acf = [0.4;-0.4];
                l = find(acf<bounds_acf(1)&acf>bounds_acf(2));
                obj.q = l(1) - 1;
            end

            [pacf,lags_pacf,bounds_pacf] = parcorr(obj.tt_obs.Q_trans_stand_d, 20);
            k = find(pacf<bounds_pacf(1)&pacf>bounds_pacf(2));
            obj.p = k(1) - 1;

            f1 = figure('Name','ACF & PACF','NumberTitle','off');
            subplot(2,1,1);
            autocorr(obj.tt_obs.Q_trans_stand_d, 100);
            text(obj.q-.5,bounds_acf(1)+.1, ['\bf q=' int2str(obj.q)], 'FontSize', 13);
            xlim([0 100]);
            hline = refline([0 bounds_acf(1)]);
            hline.Color = 'red';
            title('');
            xlabel('Lag [Days]');
            ylabel('ACF');
            subplot(2,1,2);
            parcorr(obj.tt_obs.Q_trans_stand_d, 20);
            xlim([0 20]);
            text(obj.p-.5,bounds_pacf(1)+.15, ['\bf p=' int2str(obj.p)], 'FontSize', 13);
            title('');
            xlabel('Lag [Days]');
            ylabel('PACF');
            saveas(f1,[obj.dir_results '/ACF_and_PACF.fig']);
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
            Mdl = arima(obj.p,0,obj.q);

            obj.EstMdl = estimate(Mdl,Q_trans_stand);
            obj.p = length(obj.EstMdl.AR);
            obj.q = length(obj.EstMdl.MA);
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
            IHA_ind_obs_mean=mean(IHA_ind_obs,2);
            IHA_ind_obs_mean_1 = [IHA_ind_obs_mean(1:22); IHA_ind_obs_mean(27:34)];

            acf_obs = autocorr(Q_obs, 100);

            ws = 3:2:9;
            x2 = 10:10:60;
            res = zeros(length(ws),length(x2));

            for i = 1:length(ws)
                for j = 1:length(x2)
                    Q_syn_mws = obj.MAwMWS(Q_syn,ws(i),x2(j),'lin');
                    perc.q25th=prctile(Q_syn_mws,25);
                    perc.q75th=prctile(Q_syn_mws,75);

                    [IHA_ind_syn_mws]= IHA_indicators( Q_syn_mws, perc, init_date, init_year );
                    IHA_ind_syn_mws_mean=mean(IHA_ind_syn_mws,2);
                    IHA_ind_syn_mws_mean_1 = [IHA_ind_syn_mws_mean(1:22); IHA_ind_syn_mws_mean(27:34)];
                    IHA_diff = abs(IHA_ind_obs_mean_1-IHA_ind_syn_mws_mean_1);
                    mean_IHA_diff = mean(IHA_diff);

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

        function Q = MAwMWS(obj, Q, N_ws, prc, method)
            % Filtering time series by moving average with moving window
            % size. Input arguments are time series, the starting window size,
            % streamflow percentile as far as the filtering is applied and
            % method used for spacing. Output argument is a filtered time
            % series.
            if nargin < 5
                error('Four input arguments are necessary.')
            end
            if nargin > 5
                error('Too many input arguments. Only four input argument are necessary.')
            end
            t = table();
            t.Q = Q;
            t.Q_sim_dma = zeros(length(t.Q),1);
            x2 = prctile(Q,prc);
            x1 = nanmin(Q);
            ws = zeros(ceil(N_ws/2),2);

            if (strcmp(method,'lin'))
                ws(:,1) = linspace(x1,x2,ceil(N_ws/2));
                ws(:,2) = linspace(N_ws,1,ceil(N_ws/2));
                t.ws = zeros(length(t.Q),1);
                for i = 1:ceil(N_ws/2)
                    if (i < ceil(N_ws/2))
                        t.ws(t.Q>=ws(i,1)&t.Q<ws(i+1,1)) = ws(i,2);
                    elseif (i == ceil(N_ws/2))
                        t.ws(t.Q>=ws(i,1)) = ws(i,2);
                    end
                end

            elseif (strcmp(method,'log'))
                ws(:,1) = logspace(x1,x2,ceil(N_ws/2));
                ws(:,2) = linspace(N_ws,1,ceil(N_ws/2));
                t.ws = zeros(length(t.Q),1);
                for i = 1:ceil(N_ws/2)
                    if (i < ceil(N_ws/2))
                        t.ws(log(t.Q)>=ws(i,1)&log(t.Q)<ws(i+1,1)) = ws(i,2);
                    elseif (i == ceil(N_ws/2))
                        t.ws(log(t.Q)>=ws(i,1)) = ws(i,2);
                    end
                end

            elseif (strcmp(method,'perc'))
                perc = linspace(0,prc,ceil(N_ws/2));
                ws(:,1) = prctile(Q,perc);
                ws(:,2) = linspace(N_ws,1,ceil(N_ws/2));
                t.ws = zeros(length(t.Q),1);
                for i = 1:ceil(N_ws/2)
                    if (i < ceil(N_ws/2))
                        t.ws(t.Q>=ws(i,1)&t.Q<ws(i+1,1)) = ws(i,2);
                    elseif (i == ceil(N_ws/2))
                        t.ws(t.Q>=ws(i,1)) = ws(i,2);
                    end
                end
            end

            if (N_ws == 3)
                t.ma_3 = smooth(t.Q, 3, 'moving');

                t.Q_sim_dma = t.Q;

                t.Q_sim_dma(t.ws==3) = t.ma_3(t.ws==3);

            elseif (N_ws == 5)
                t.ma_3 = smooth(t.Q, 3, 'moving');
                t.ma_5 = smooth(t.Q, 5, 'moving');

                t.Q_sim_dma = t.Q;

                t.Q_sim_dma(t.ws==3) = t.ma_3(t.ws==3);
                t.Q_sim_dma(t.ws==5) = t.ma_5(t.ws==5);

            elseif (N_ws == 7)
                t.ma_3 = smooth(t.Q, 3, 'moving');
                t.ma_5 = smooth(t.Q, 5, 'moving');
                t.ma_7 = smooth(t.Q, 7, 'moving');

                t.Q_sim_dma = t.Q;

                t.Q_sim_dma(t.ws==3) = t.ma_3(t.ws==3);
                t.Q_sim_dma(t.ws==5) = t.ma_5(t.ws==5);
                t.Q_sim_dma(t.ws==7) = t.ma_7(t.ws==7);

            elseif (N_ws == 9)
                t.ma_3 = smooth(t.Q, 3, 'moving');
                t.ma_5 = smooth(t.Q, 5, 'moving');
                t.ma_7 = smooth(t.Q, 7, 'moving');
                t.ma_9 = smooth(t.Q, 9, 'moving');

                t.Q_sim_dma = t.Q;

                t.Q_sim_dma(t.ws==3) = t.ma_3(t.ws==3);
                t.Q_sim_dma(t.ws==5) = t.ma_5(t.ws==5);
                t.Q_sim_dma(t.ws==7) = t.ma_7(t.ws==7);
                t.Q_sim_dma(t.ws==9) = t.ma_9(t.ws==9);
            end
            Q = t.Q_sim_dma;

            mws_1 = [0;ws(1,2)+2]';
            mws_2 = [max(obj.tt_syn.Q_sim_re);0]';
            mws = [mws_1;ws;mws_2];
            f1 = figure('Name','Moving average with moving window size','NumberTitle','off');
            stairs(mws(:,1), mws(:,2), 'k');
            grid;
            yticks(flipud(mws(:,2)));
            xlabel('Q [m^3/s]');
            ylabel('Window size [Days]');
            saveas(f1,[obj.dir_results '/MAwMWS.fig'])
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

            x1_nn = min(xlims_nn(:,1));
            x2_nn = max(xlims_nn(:,2));

            xlims_n = zeros(3,2);
            xlims_n(1,1) = min(obj.tt_obs.Q_trans);
            xlims_n(2,1) = min(obj.tt_obs.Q_trans_stand_d);
            xlims_n(3,1) = min(obj.tt_syn.Q_sim);
            xlims_n(1,2) = max(obj.tt_obs.Q_trans);
            xlims_n(2,2) = max(obj.tt_obs.Q_trans_stand_d);
            xlims_n(3,2) = max(obj.tt_syn.Q_sim);

            x1_n = min(xlims_n(:,1));
            x2_n = max(xlims_n(:,2));

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
            y2_nn = round(max(ylims_nn(:,1)),1)+.1;

            bins_nn = zeros(2,1);
            bins_nn(1,1) = h1nb;
            bins_nn(2,1) = h5nb;
            nbins_nn = max(bins_nn(:,1));

            ylims_n = zeros(3,1);
            ylims_n(1,1) = max(h2v);
            ylims_n(2,1) = max(h3v);
            ylims_n(3,1) = max(h4v);
            y2_n = round(max(ylims_n(:,1)),1)+.02;

            bins_n = zeros(3,1);
            bins_n(1,1) = h2nb;
            bins_n(2,1) = h3nb;
            bins_n(3,1) = h4nb;
            nbins_n = max(bins_n(:,1));

            f1 = figure('Name','Histograms','NumberTitle','off');
            subplot(5, 1, 1);
            histogram(obj.tt_obs.Q,nbins_nn,'Normalization','probability','FaceColor','b');
            grid; title('Nonnormally Distributed');
            xlim([x1_nn x2_nn]);
            ylim([0 y2_nn]);
            subplot(5, 1, 2);
            histogram(obj.tt_obs.Q_trans,nbins_n,'Normalization','probability','FaceColor','b');
            grid; title('Transformed-to-normal');
            xlim([x1_n x2_n]);
            ylim([0 y2_n]);
            subplot(5, 1, 3);
            histogram(obj.tt_obs.Q_trans_stand_d,nbins_n,'Normalization','probability','FaceColor','b');
            grid; title('Transformed-to-normal & Standardized (Observed)');
            xlim([x1_n x2_n]);
            ylim([0 y2_n]);
            ylabel('Probability');
            subplot(5, 1, 4);
            histogram(obj.tt_syn.Q_sim,nbins_n,'Normalization','probability','FaceColor','r');
            grid; title(['Transformed-to-normal & Standardized (' obj.approach_str{obj.approach} ')']);
            xlim([x1_n x2_n]);
            ylim([0 y2_n]);
            subplot(5, 1, 5);
            histogram(obj.tt_syn.Q_sim_re,nbins_nn,'Normalization','probability','FaceColor','r');
            grid; title(obj.approach_str{obj.approach});
            xlim([x1_nn x2_nn]);
            ylim([0 y2_nn]);
            xlabel('Q [m^3/s]');
            saveas(f1,[obj.dir_results '/Histograms.fig'])

            f2 = figure('Name','Histograms','NumberTitle','off');
            s1 = subplot(2, 1, 1);
            histogram(obj.tt_obs.Q,nbins_nn,'Normalization','probability','FaceColor','b');
            grid; title('Nonnormally Distributed');
            xlim([x1_nn x2_nn]);
            ylim([0 y2_nn]);
            s2 = subplot(2, 1, 2);
            histogram(obj.tt_syn.Q_sim_re,nbins_nn,'Normalization','probability','FaceColor','r');
            grid; title(obj.approach_str{obj.approach});
            xlim([x1_nn x2_nn]);
            ylim([0 y2_nn]);
            xlabel('Q [m^3/s]');
            p1=get(s1,'position');
            p2=get(s2,'position');
            height=p1(2)+p1(4)-p2(2);
            axes('position',[p2(1) p2(2) p2(3) height],'visible','off');
            ylabel('Probability','visible','on');
            %             saveas(f1,[obj.dir_results '/Histograms.pdf'])
            saveas(f2,[obj.dir_results '/Hist_obs_vs_syn.fig'])

            f3 = figure('Name','Histogram of Residuals','NumberTitle', 'off');
            histogram(obj.tt_syn.E,'FaceColor','k','Normalization','probability');
            grid; title('Distribution of \epsilon');
            xlabel('Q [m^3/s]');
            ylabel('Probability');
            saveas(f3,[obj.dir_results '/Histogram_of_Residuals.fig']);
        end

        function testautocor(obj)
            % Testing visually for correlation in the Residuals and comapring
            % the ACF of the observed and the synthetic streamflow and
            % the PACF respectively.
            if nargin >= 2
               error('Too many input arguments. No input argument is necessary.')
            end
            f1 = figure('Name','ACF of Residuals','NumberTitle','off');
            autocorr(obj.tt_syn.E, 100);
            xlim([0 100]);
            title('');
            xlabel('Lag [Days]');
            ylabel('ACF');
            saveas(f1,[obj.dir_results '/ACF_of_Residuals.fig']);

            f2 = figure('Name','ACF of Obs & Syn','NumberTitle','off');
            subplot(2,1,1);
            autocorr(obj.tt_obs.Q, 100);
            title('Observed');
            xlim([0 100]);
            xlabel('Lag [Days]');
            ylabel('ACF');
            subplot(2,1,2);
            autocorr(obj.tt_syn.Q_sim_re, 100);
            title(obj.approach_str{obj.approach});
            xlim([0 100]);
            xlabel('Lag [Days]');
            ylabel('ACF');
            saveas(f2,[obj.dir_results '/ACF_of_Obs_and_Syn.fig']);

            f3 = figure('Name','PACF of Obs & Syn','NumberTitle','off');
            subplot(2,1,1);
            parcorr(obj.tt_obs.Q, 20);
            xlim([0 20]);
            title('Observed');
            xlabel('Lag [Days]');
            ylabel('PACF');
            subplot(2,1,2);
            parcorr(obj.tt_syn.Q_sim_re, 20);
            xlim([0 20]);
            title(obj.approach_str{obj.approach});
            xlabel('Lag [Days]');
            ylabel('PACF');
            saveas(f3,[obj.dir_results '/PACF_Obs_vs_Syn.fig']);
        end

        function ACFmonths(obj)
            % Splitting the runoff time series into single months,
            % calculate the ACF and compare between observed and synthetic
            % streamflow.
            if nargin >= 2
               error('Too many input arguments. No input argument is necessary.')
            end
            for i = 1:12
                f = figure('Name',['ACF of Obs & Syn ' num2str(i)],'NumberTitle','off','defaultFigureVisible','off');
                subplot(2,1,1);
                autocorr(obj.tt_obs.Q(obj.tt_obs.MM==i), 100);
                title('Observed');
                xlim([0 100]);
                xlabel('Lag [Days]');
                ylabel('ACF');
                subplot(2,1,2);
                autocorr(obj.tt_syn.Q_sim_re(obj.tt_syn.MM==i), 100);
                title(obj.approach_str{obj.approach});
                xlim([0 100]);
                xlabel('Lag [Days]');
                ylabel('ACF');
                saveas(f,[obj.dir_results '/ACF_of_Obs_and_Syn_' num2str(i) '.fig']);
            end
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
            IHA_ind_obs_mean=mean(IHA_ind_obs,2);

            perc.q25th=prctile(obj.tt_syn.Q_sim_re,25);
            perc.q75th=prctile(obj.tt_syn.Q_sim_re,75);

            init_date=date;
            init_year=year;

            [IHA_ind_syn]= IHA_indicators(obj.tt_syn.Q_sim_re, perc, init_date, init_year);
            IHA_ind_syn_mean=mean(IHA_ind_syn,2);

            f1 = figure('Name','IHA - Group 1','NumberTitle','off');
            hold on
            b = bar([IHA_ind_obs_mean(1:12) IHA_ind_syn_mean(1:12)]);
            b(1).FaceColor = 'b';
            b(2).FaceColor = 'r';
            title('Group 1')
            ylabel('Q [m^3/s]');
            legend({'Q_{obs}','Q_{syn}'},'Box','off')
            xticks(1:1:12);
            saveas(f1,[obj.dir_results '/IHA_group_1.fig']);

            f2 = figure('Name','IHA - Group 2','NumberTitle','off');
            hold on
            b = bar([IHA_ind_obs_mean(13:22) IHA_ind_syn_mean(13:22)]);
            b(1).FaceColor = 'b';
            b(2).FaceColor = 'r';
            title('Group 2')
            ylabel('Q [m^3/s]');
            legend({'Q_{obs}','Q_{syn}'},'Box','off')
            xp = 1:1:10;
            xt = {{'Min'; '1-day'} {'Min'; '3-day'} ...
                {'Min'; '7-day'} {'Min'; '30-day'}...
                {'Min'; '90-day'} {'Max'; '1-day'}...
                {'Max'; '3-day'} {'Max'; '7-day'}...
                {'Max'; '30-day'} {'Max'; '90-day'}};
            my_xticklabels(xp, xt);
            saveas(f2,[obj.dir_results '/IHA_group_2.fig']);

            f3 = figure('Name','IHA - Group 4','NumberTitle','off');
            hold on
            b = bar([IHA_ind_obs_mean(27:30) IHA_ind_syn_mean(27:30)]);
            b(1).FaceColor = 'b';
            b(2).FaceColor = 'r';
            title('Group 4')
            legend({'Q_{obs}','Q_{syn}'},'Box','off')
            xp = 1:1:4;
            xt = {{'No. of low pulses'} {'No. of high pulses'} {'Mean duration'; 'of low pulses'} {'Mean duration'; 'of high pulses'}};
            my_xticklabels(xp, xt);
            saveas(f3,[obj.dir_results '/IHA_group_4.fig']);

            f4 = figure('Name','IHA - Group 5-1','NumberTitle','off');
            hold on
            b = bar([IHA_ind_obs_mean(31:32) IHA_ind_syn_mean(31:32)]);
                        b(1).FaceColor = 'b';
            b(2).FaceColor = 'r';
            title('Group 5 - 1')
            legend({'Q_{obs}','Q_{syn}'},'Box','off')
            xp = [1 2];
            xt = {{'Means of all negative'; 'differences between'; 'consecutive daily values'} ...
                {'Means of all postive'; 'differences between'; 'consecutive daily means'}};
            my_xticklabels(xp, xt);
            saveas(f4,[obj.dir_results '/IHA_group_5_1.fig']);

            f5 = figure('Name','IHA - Group 5-2','NumberTitle','off');
            hold on
            b = bar([IHA_ind_obs_mean(33:34) IHA_ind_syn_mean(33:34)]);
            b(1).FaceColor = 'b';
            b(2).FaceColor = 'r';
            title('Group 5 - 2')
            legend({'Q_{obs}','Q_{syn}'},'Box','off')
            xticks([1 2])
            xticklabels({'No. of falls', 'No. of rises'})
            saveas(f5,[obj.dir_results '/IHA_group_5_2.fig']);
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

           tt_vol_mm = retime(tt_vol, 'monthly', 'sum');
           tt_vol_mm.vol_obs_cum = cumsum(tt_vol_mm.vol_obs);
           tt_vol_mm.vol_syn_cum = cumsum(tt_vol_mm.vol_syn);
           tt_vol_mm.MM = month(tt_vol_mm.Time);
           func = @sum;
           tt_vol_mm_sum = varfun(func,tt_vol_mm,'GroupingVariables','MM');
           tt_vol_yy = retime(tt_vol, 'yearly', 'sum');
           tt_vol_yy.vol_obs_cum = cumsum(tt_vol_yy.vol_obs);
           tt_vol_yy.vol_syn_cum = cumsum(tt_vol_yy.vol_syn);

           f1 = figure('Name','Cumulated Volume (monthly)','NumberTitle','off');
           hold on
           plot(tt_vol_mm.Time, tt_vol_mm.vol_obs_cum, 'b');
           plot(tt_vol_mm.Time, tt_vol_mm.vol_syn_cum, 'r');
           hold off
           grid;
           xlabel('Date');
           ylabel('Volume [km^3]');
           xtickformat('MMM-yyyy');
           legend({'Q_{obs}','Q_{syn}'},'Box','off','Location','northwest');
           saveas(f1,[obj.dir_results '/cum_volume_mm.fig']);

           f2 = figure('Name','Cumulated Volume (yearly)','NumberTitle','off');
           hold on
           plot(tt_vol_yy.Time, tt_vol_yy.vol_obs_cum, 'b');
           plot(tt_vol_yy.Time, tt_vol_yy.vol_syn_cum, 'r');
           hold off
           grid;
           xlabel('Date');
           ylabel('Volume [km^3]');
           xtickformat('yyyy');
           legend({'Q_{obs}','Q_{syn}'},'Box','off','Location','northwest');
           saveas(f2,[obj.dir_results '/cum_volume_yy.fig']);

           f3 = figure('Name','Total Volume (monthly)','NumberTitle','off');
           hold on
           b = bar([tt_vol_mm_sum.sum_vol_obs tt_vol_mm_sum.sum_vol_syn]);
           b(1).FaceColor = 'b';
           b(2).FaceColor = 'r';
           legend({'Q_{obs}','Q_{syn}'},'Box','off')
           xlabel('Month');
           ylabel('Volume [km^3]');
           xticks(1:1:12);
           saveas(f3,[obj.dir_results '/total_volume_mm.fig']);
        end

        function plotts(obj)
           % Plotting the generated time series.
           if nargin >= 2
              error('Too many input arguments. No input argument is necessary.')
           end
           f1 = figure('Name','Q_obs vs Q_syn (first year)','NumberTitle','off','defaultFigureVisible','off');
           hold on
           plot(obj.tt_obs.Date, obj.tt_obs.Q, 'b');
           plot(obj.tt_syn.Date, obj.tt_syn.Q_sim_re, 'r');
           hold off
           grid;
           xlabel('Date');
           ylabel('Q [m^3/s]');
           xlim(datetime(obj.tt_obs.YYYY(1),[1 12],[1 31]));
           xtickformat('dd-MMM-yyyy');
           legend({'Q_{obs}','Q_{syn}'},'Box','off');
           saveas(f1,[obj.dir_results '/Q_obs_vs_Q_syn_' num2str(obj.tt_obs.YYYY(1)) '.fig']);

           f2 = figure('Name','Q_obs vs Q_syn (entire)','NumberTitle','off','defaultFigureVisible','off');
           hold on
           plot(obj.tt_obs.Date, obj.tt_obs.Q, 'b');
           plot(obj.tt_syn.Date, obj.tt_syn.Q_sim_re, 'r');
           hold off
           grid;
           xlabel('Date');
           ylabel('Q [m^3/s]');
           xtickformat('dd-MMM-yyyy');
           legend({'Q_{obs}','Q_{syn}'},'Box','off');
           saveas(f2,[obj.dir_results '/Q_obs_vs_Q_syn.fig']);

           f3 = figure('Name','Q_syn (entire)','NumberTitle','off','defaultFigureVisible','off');
           plot(obj.tt_syn.Date, obj.tt_syn.Q_sim_re, 'r');
           grid;
           xlabel('Date');
           ylabel('Q [m^3/s]');
           xtickformat('dd-MMM-yyyy');
           legend({'Q_{obs}','Q_{syn}'},'Box','off');
           saveas(f3,[obj.dir_results '/Q_syn.fig']);

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

            f1 = figure('Name','Test statistic','NumberTitle','off');
            b = bar([obj.test_stats(:,1) obj.test_stats(:,2)]);
            b(1).FaceColor = 'b';
            b(2).FaceColor = 'r';
            grid;
            legend({'Q_{obs}','Q_{syn}'},'Box','off','Location','northwest')
            ylabel('Q [m^3/s]');
            xticks([1 2 3 4 5])
            xticklabels({'Min', 'Max',' Mean', 'Std', 'Skew'})
            saveas(f1, [obj.dir_results '/test_stats.fig']);
        end

        function teststatstotxt(obj)
            % Export test statistic to .txt..
            if nargin >= 2
               error('Too many input arguments. No input argument is necessary.')
            end
            B = array2table(obj.test_stats);
            B.Properties.VariableNames = {'Obs' 'Sim'};
            B.Properties.RowNames = {'mean' 'std' 'skew' 'min' 'max'};
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
