function IHA_JD = mean_timing_var(JDs)
% :::Help:::
%
% Input argument:
%       JDs: vector containing the dates of the seasonal max/min of daily mean streamflow, one for each year.  -vector [1xn°years];
% Output argument:
%       IHA_JD: Julian Date of the mean calculated on the values of the JDs vector, using the circular method, as explained below.
%               The Julian Date value minimum is 1 and the maximum is 365.                                     - double number.
%                                                                                      
% In order to deal with the difficulties of calculating statistics for timing variables, the IHA often
% uses a "circular method" for timing statistics. This method puts all the dates in question into
% quarterly bins, which are 1-91, 92-183, 184-275, and 276-366 (depending on the situation,
% these bins are applied to either Julian dates or dates of the water year). If the second or third
% quarter has the greatest frequency of dates, the statistic is computed in the usual way. If the
% largest number of dates is in the first quarter, 366 is temporarily subtracted from all dates in the
% fourth quarter prior to computing the statistic, and if the result is less than 0 then 366 is added
% back onto it. If the largest number of dates is in the fourth quarter, 366 is temporarily added to
% all dates in the first quarter prior to computing the statistics, and if the result is greater than 366
% then 366 is subtracted from it.
%
%
% Andrea Cominola, Emanuele Mason 28/09/2010
JDmax = sort(JDs);
count(1) = sum(JDmax <= 90);
count(2) = sum((JDmax >= 91)&(JDmax <=182));
count(3) = sum((JDmax >= 183)&(JDmax <=274));
count(4) = sum(JDmax >= 275);
if find(count == max(count)) == 1; JDmax(JDmax >= 275) = JDmax(JDmax >= 275) - 365; end;
if find(count == max(count)) == 4; JDmax(JDmax <= 90) = JDmax(JDmax <= 90) + 365; end;
IHA_JD = mean(JDmax);
if IHA_JD < 0; IHA_JD = IHA_JD + 365; end;
if IHA_JD > 365; IHA_JD = IHA_JD - 365; end;