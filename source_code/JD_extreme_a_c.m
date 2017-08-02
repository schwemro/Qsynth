function IHA_g3=JD_extreme_a_c(q)
% :::HELP:::
%
% IHA_g3= extreme_a_c (q)
%
% The function calculates the Julian Date values of annual maxima (spring and autumn)and
% minima (summer and winter) 1-day streamflow means. These are the parameters of group 3 of IHA
% indicators: "timing of annual extreme water conditions".
% If the maximum/minimum streamflow 1-day mean is found more than one time in each season in the data serie,
% the function gives as output only the first date found.
%
% Input arguments:
%       q: vector of daily streamflow means [m^3/s]        -vector [365 or 366x1]
%
% Output arguments:
%       IHA_g3: structure containing 4 fields:
%               IHA_g3.JD_max_spr: field containing the Julian Date of the spring maximum
%               1_day streamflow mean.                                     -integer number 
%               IHA_g3.JD_max_aut: field containing the Julian Date of the spring maximum
%               1_day streamflow mean.                                     -integer number
%               IHA_g3.JD_min_sum: field containing the Julian Date of the summer minimum
%               1_day streamflow mean.                                     -integer number
%               IHA_g3.JD_min_win: field containing the Julian Date of the winter minimum
%               1_day streamflow mean.                                     -integer number
%
% Julian Dates in the IHA method:
% in order to quantify the timings of maximum and minimum flows, the IHA uses
% the concept of Julian dates. Julian dates represent calendar dates by integer values, which start
% with 1 on January 1 and end with 366 on December 31. Note that although the number of
% calendar days varies slightly in each year, depending on whether or not it is a leap year, the
% starting and ending Julian dates in each year are always 1 and 366. In the IHA the difference
% between leap years and non-leap years is that leap years have a Julian date for Feb. 29, which
% is 60, while non-leap years skip from Julian date 59 (Feb. 28) to Julian date 61 (Mar. 1). This
% ensures that each calendar date is represented by the same Julian date in each year.
%
% Users should be aware that this method of defining Julian dates does differ from definitions
% used elsewhere, which usually assign a Julian date of 1 to January 1, 4713 B.C., and then
% increment the date by 1 for each day thereafter, without restarting the count at the beginning of
% each year.
%
% Julian Dates calculated by the function:
% the dates calculated by the function always start with 1 on January 1 and end with 365 on December 31.
% In case of leap year, the 29th of February is given the same Julian Date
% of the 28th of February (59), so that the 1st of March always assumes the
% Julian Date value of 60.
% Before calculating the dates, the function uses the funtion
% "moving_average" with an interval of 6 days on the serie, to produce a
% new serie of daily streamflow mean.
%
%
% Andrea Cominola, Emanuele Mason 28/09/2010

N = length(q);
f=6;
if (N> 366)||(N<365) 
  disp('Error:the timing_ex_a_c function requires a 365x1 or 366x1 streamflow means vector');
  return;
end
% :::Spring and autumn max:::
q_spr=q(1:275);
q_aut=q(276:end);

% Not leap year case
q_mm_spr = moving_average(q_spr,f);
q_mm_aut = moving_average(q_aut,f);
q_max_spr = max(q_mm_spr);
q_max_aut = max(q_mm_aut);

tot_JD_max_spr=find(q_mm_spr==q_max_spr);
tot_JD_max_aut=find(q_mm_aut==q_max_aut)+276;

IHA_g3.JD_max_spr = tot_JD_max_spr(1)+(f/2)-1;    %positions of the max value in q_mm_spr vector.
IHA_g3.JD_max_aut = tot_JD_max_aut(1)+(f/2)-1;    %positions of the max value in q_mm_aut vector.

% Leap year case
if (N==366) && (IHA_g3.JD_max_spr>59)
    IHA_g3.JD_max_spr=IHA_g3.JD_max_spr-1;
end
if N==366
    IHA_g3.JD_max_aut=IHA_g3.JD_max_aut-1;
end

% :::Summer minimun:::
q_sum=q(IHA_g3.JD_max_spr:IHA_g3.JD_max_aut);
q_mm_sum = moving_average(q_sum,f);
q_min_sum=min(q_mm_sum);
tot_JD_min_sum=find(q_mm_sum==q_min_sum)+IHA_g3.JD_max_spr;
IHA_g3.JD_min_sum = tot_JD_min_sum(1)+(f/2)-1;    %positions of the max value in q_mm_sum vector.

%Leap year case
 if (N==366) && (IHA_g3.JD_min_sum>59)
    IHA_g3.JD_min_sum=IHA_g3.JD_min_sum-1;
 end
 
 
%:::Winter minimun:::
q_win=[q(IHA_g3.JD_max_aut:end);q(1:IHA_g3.JD_max_spr)];
q_mm_win = moving_average(q_win,f);
q_min_win=min(q_mm_win);
tot_JD_min_win=find(q_mm_win==q_min_win);
IHA_g3.JD_min_win = tot_JD_min_win(1)+(f/2)-1;    %positions of the max value in q_mm_win vector.
if IHA_g3.JD_min_win < length(q(IHA_g3.JD_max_aut:end))
    %leap year case
    if N==366
        IHA_g3.JD_min_win = IHA_g3.JD_min_win + IHA_g3.JD_max_aut -2;
    %not leap year case
    else
        IHA_g3.JD_min_win = IHA_g3.JD_min_win + IHA_g3.JD_max_aut -1;
    end
else
    %leap year case
    if (N==366) && (IHA_g3.JD_min_win>59)
        IHA_g3.JD_min_win = IHA_g3.JD_min_win + IHA_g3.JD_max_aut - N -2;
    %not leap year case
    else
        IHA_g3.JD_min_win = IHA_g3.JD_min_win + IHA_g3.JD_max_aut - N -1;
    end
end

end

    
    
    