function [IHA_ind] = IHA_indicators( signal, perc, init_date_sim, init_year )
% :::Help:::
%
% IHA_ind = IHA_indicators( signal, perc, init_date_sim, init_year )
%
% The function calculates the values of all 34 IHA indicators, using the
% functions created for each group.
% The indicators are intra-annual indicators, so the function evaluates
% them only for the years in which there is complete availability of data
% (daily mean streamflow present for everyday of the year). 
% 
% Input arguments:
%       signal: daily mean streamflow vector  [m^3/s]          -vector [Nx1]
%       perc:values of the 25th and 75th percentile,        -double number
%       init_date_sim: vector containing the initial date of the data serie
%                      loaded in "signal"                      -vector [1x3]
%                      Example : 
%                      init_date_sim = [ 01 01 1974 ]     ; % first simulation day in the out file
%       init_year: value of the year chosen for the beginning of the
%                  evaluation of the indicators.               -integer number
%                  Example: 1980
% 
% Output arguments:
%       IHA_ind: matrix containing the values of the 34 intra-annual indicators.
%                It has: -number of rows = 34 (number of indicators);
%                        -number of columns = number of years of the serie.
%                For each column:
%                        -in rows from 1 to 12 there are the values of IHA
%                         indicators group 1;
%                        -in rows from 13 to 22 there are the values of IHA
%                         indicators group 2;
%                        -in rows from 23 to 26 there are the values of IHA
%                         indicators group 3;
%                        -in rows from 27 to 30 there are the values of IHA
%                         indicators group 4;
%                        -in rows from 31 to 34 there are the values of IHA
%                         indicators group 5;
%                                                               -matrix [34xn°years] 
%
% Needed functions:
%   -time_date2JD
%   -monthy_means
%   -extreme_a_c 
%   -JD_extreme_a_c
%   -pulses
%   -water_c_changes
%
% Andrea Cominola, Emanuele Mason 29/08/2010
% Lorenzo Gorla 05/01/2012

if init_date_sim(1,3) > init_year
    disp( 'Error: the 1st year in the evaluation horizon drops before the 1st day in the data' );
    return;
end

if (init_date_sim(1,3) == init_year &&((init_date_sim(1,1) ~= 1)||(init_date_sim(1,2) ~= 1)))
    disp( 'Error: the 1st year in the evaluation horizon is not complete:');
    disp( 'the IHA indicators will be calculated starting from the next year' );
    init_year = init_year +1;
end

% Years included in the data
last_day = time_JD2date(time_date2JD(init_date_sim) + length(signal));
if ((last_day(1) == 31) && (last_day(2) == 12))
    last_year = last_day(3);
else
    % disp( ['Warning: the last year (' last_day(3) ') data serie is not complete, so it will be computed the previous one' ]);
    last_year = last_day(3) - 1;
end

N = last_year-init_year+1;
eval_years = (init_year:last_year)';
idx_in  = time_date2JD([ones(N,1) ones(N,1) eval_years]) - time_date2JD( init_date_sim ) + 1 ;
idx_fin = time_date2JD([ones(N,1)*31 ones(N,1)*12 eval_years])  - time_date2JD( init_date_sim ) + 1 ;

% IHA indicators
IHA_ind = NaN(34,N);
for i=1:N
    %% group 1
    IHA_ind(1:12,i) = monthy_means(signal(idx_in(i):idx_fin(i)));
    %% group 2
    temp = extreme_a_c( signal(idx_in(i):idx_fin(i)));
    IHA_ind(13,i) = temp.day1_min;    
    IHA_ind(14,i) = temp.day3_min;    
    IHA_ind(15,i) = temp.day7_min;    
    IHA_ind(16,i) = temp.day30_min;    
    IHA_ind(17,i) = temp.day90_min;    
    IHA_ind(18,i) = temp.day1_max;
    IHA_ind(19,i) = temp.day3_max;
    IHA_ind(20,i) = temp.day7_max;
    IHA_ind(21,i) = temp.day30_max;
    IHA_ind(22,i) = temp.day90_max;
 

    %% group 3
    temp = JD_extreme_a_c( signal(idx_in(i):idx_fin(i)) );
    IHA_ind(23,i) = temp.JD_max_spr;
    IHA_ind(24,i) = temp.JD_min_sum;
    IHA_ind(25,i) = temp.JD_max_aut;
    IHA_ind(26,i) = temp.JD_min_win;
    
    %% group 4
    [temp, mem] = pulses_mio( signal(idx_in(i):idx_fin(i)), perc.q25th, perc.q75th);
    IHA_ind(27,i) = temp.n_low_p;    
    IHA_ind(28,i) = temp.n_high_p;
    IHA_ind(29,i) = temp.sum_durat_low_p;
    IHA_ind(30,i) = temp.sum_durat_high_p;
    
    % tab_low contains [num low periods; bool low cell 1; durat low event 1; bool last low
    % cell; durat last low event] for every year
    % the same for tab_high
    
    tab_low(:,i)=mem.low(:,1);
    tab_high(:,i)=mem.high(:,1);   
    
    %% group 5
    temp = water_c_changes( signal(idx_in(i):idx_fin(i)) );
    IHA_ind(31,i) = temp.mean_neg_diff;
    IHA_ind(32,i) = temp.mean_pos_diff;
    IHA_ind(33,i) = temp.n_falls;
    IHA_ind(34,i) = temp.n_rises;

end


%% correctons of indicators Gr4 folloving convenctions of IHA manual

%if the year starts in a low/high period, and that period started in the previous year, that period doesn't count for
%that year (it will count in the previous one)

%1th year: we suppose that, if the 1th year starts in a low/high period, that period started in the previous year
if (tab_low(1,1)~=0) && (tab_low(2,1)==1)
    IHA_ind(27,1)=IHA_ind(27,1)-1;
    IHA_ind(29,1)=IHA_ind(29,1)-tab_low(3,1);
end
if (tab_high(1,1)~=0) && (tab_high(2,1)==1)
    IHA_ind(28,1)=IHA_ind(28,1)-1;
    IHA_ind(30,1)=IHA_ind(30,1)-tab_high(3,1);
end

for i=2:N
if (tab_low(1,i)~=0) && (tab_low(2,i)==1) && (tab_low(4,i-1)==1)
    IHA_ind(27,i)=IHA_ind(27,i)-1;
    IHA_ind(29,i)=IHA_ind(29,i)-tab_low(3,i);
end
if (tab_high(1,i)~=0) && (tab_high(2,i)==1) && (tab_high(4,i-1)==1)
    IHA_ind(28,i)=IHA_ind(28,i)-1;
    IHA_ind(30,i)=IHA_ind(30,i)-tab_high(3,i);
end
end

% if the year i finishes in a low/high period and the year i+1 starts in a
% low/high period, then the sum of the two will be counted only in the year
% when it started
for i=1:(N-1)
if (tab_high(1,i)~=0) && (tab_high(4,i)==1) && (tab_high(2,i+1)==1)    
    IHA_ind(30,i) = (IHA_ind(30,i) +  tab_high(3,i+1))/IHA_ind(28,i);
else
    IHA_ind(30,i) = IHA_ind(30,i)/IHA_ind(28,i);
end
    
if (IHA_ind(29,i)==0 || IHA_ind(27,i)==0)
    IHA_ind(29,i) = 0;
elseif (tab_low(1,i)~=0) && (tab_low(4,i)==1) && (tab_low(2,i+1)==1)    
    IHA_ind(29,i) = (IHA_ind(29,i) +  tab_low(3,i+1))/IHA_ind(27,i);
else
    IHA_ind(29,i) = IHA_ind(29,i)/IHA_ind(27,i);
end    
end
IHA_ind(30,end) = IHA_ind(30,end)/IHA_ind(28,end);
if (IHA_ind(29,i)==0 || IHA_ind(27,i)==0)
    IHA_ind(29,end) = 0;
else
    IHA_ind(29,i) = IHA_ind(29,i)/IHA_ind(27,i);
end

end
