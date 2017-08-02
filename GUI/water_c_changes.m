function [IHA_g5] = water_c_changes(q)
% ::: Help :::
%
% [IHA_g5] = water_c_changes(q)
% 
% Evaluates the four intra-annual indicators in the 5th group of IHA
% statistics. The group is called "rate and frequency of water condition
% changes".
% This function:
% 1)    calculates the means of all positive differences between consecutive
% daily values within a year and the means of negative ones.
% 2)    calculates the number of rises and falls in the 1-year serie of
% streamflow values. It evaluates the difference between every
% consecutive daily streamflow values; when a positive (negative)
% difference is found, it is counted as a rise (fall); each consecutive
% positive (negative) or null difference doesn't increase the rises (falls)
% counter.
% 
% Input argument:
%   q = vector of daily streamflow means [m^3/s] - vector [365 or 366x1]
%   
% Output argument:
%   IHA_g5 = structure with four IHA-indicators that belong to group 5.
%   It contains four fields:
%    IHA_g5.mean_pos_diff = the mean of all positive differences between
%       consecutive daily values [m^3/s]         - double number
%    IHA_g5.mean_neg_diff = the mean of all negative differences between
%       consecutive daily values [m^3/s]         - double number
%    IHA_g5.n_rises = the number of rises in the daily streamflow serie 
%                                                - double number
%    IHA_g5.n_falls = the number of falls in the daily streamflow serie 
%                                                - double number
%
% Andrea Cominola, Emanuele Mason, Lorenzo Gorla 06/01-2012

if (size(q,1) < 365) || (size(q,1) > 366);
    disp(' Error: q = vector of daily streamflow means [m3/s]; - vector [365 or 366x1]');
    return;
end

q_d = diff(q);
q_d=q_d(q_d ~= 0);
IHA_g5.mean_pos_diff = mean( q_d( q_d > 0 ));
IHA_g5.mean_neg_diff = mean( q_d( q_d < 0 ));

bool_rise = q_d > 0;
bool_fall = q_d < 0;
% Note: each consecutive negative or null difference doesn't increase
% the falls counter.
null = find(q_d(1:end-1) == 0);
bool_rise(null(bool_rise(null+1) == 1)) = 1;
bool_fall(null(bool_fall(null+1) == 1)) = 1;

% The function is not "IHA_g5.n_rises = sum(diff([0; bool_rise; 0]) ==
% 1);", as Andrea Cominola wrote, because "Reversals are analysed on a
% water year by water year basis, so the first change in flow of the water
% cannot be counted as a reversal, since no rising or falling trend exists
% before then" (Manual, pag 7)
IHA_g5.n_rises = sum(diff([bool_rise; 0]) == 1);
IHA_g5.n_falls = sum(diff([bool_fall; 0]) == 1);