function IHA_g2= extreme_a_c (q)
% :::HELP:::
% 
% IHA_g2= extreme_a_c (q)
%
% The function calculates the values of annual maxima and minima streamflow means, for
% different time intervals. These are the parameters of group 2 of IHA
% indicators: "magnitude and duration of annual extreme water conditions".
% The function requires the use of the function "moving_average".
%
% Input arguments:
%       q: vector of daily streamflow means [m^3/s]        -vector [365 or 366x1]
%
% Output arguments:
%       IHA_g2: structure containing 10 fields:-5 fields for annual maxima values of streamflow means, calculated with a moving average with
%                                                 different time intervals (1,3,7,30,90 days)

%                                              -5 fields for annual minima values of streamflow means, calculated with a moving average with
%                                                 different time intervals (1,3,7,30,90 days)
%                    
%       IHA_g2.day1_max:annual maxima  1_day means [m^3/s]
%       IHA_g2.day3_max:annual maxima  3_day means [m^3/s]
%       IHA_g2.day7_max:annual maxima  7_day means [m^3/s]
%       IHA_g2.day30_max:annual maxima 30_day means [m^3/s]
%       IHA_g2.day90_max:annual maxima 90_day means [m^3/s]
%
%       IHA_g2.day1_min:annual minima  1_day means [m^3/s]
%       IHA_g2.day3_min:annual minima  3_day means [m^3/s]
%       IHA_g2.day7_min:annual minima  7_day means [m^3/s]
%       IHA_g2.day30_min:annual minima 30_day means [m^3/s]
%       IHA_g2.day90_min:annual minima 90_day means [m^3/s]
%
% Andrea Cominola, Emanuele Mason 28/09/2010

if (length(q)> 366)||(length(q)<365) 
  disp('Error:the extreme_a_c function requires a 365x1 or 366x1 streamflow means vector');
  return;
end
 
%Maxima values
       IHA_g2.day1_max= max(q);
       IHA_g2.day3_max=max(moving_average(q, 3));
       IHA_g2.day7_max=max(moving_average(q, 7));
       IHA_g2.day30_max=max(moving_average(q, 30));
       IHA_g2.day90_max=max(moving_average(q, 90)); 
       
%Minima values

       IHA_g2.day1_min= min(q);
       IHA_g2.day3_min=min(moving_average(q, 3));
       IHA_g2.day7_min=min(moving_average(q, 7));
       IHA_g2.day30_min=min(moving_average(q, 30));
       IHA_g2.day90_min=min(moving_average(q, 90)); 
end
  