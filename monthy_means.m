function [IHA_g1] = monthy_means(q)
% :::HELP:::
% 
% [IHA_g1] = monthy_means(q)
% 
% Creates a vector of monthy streamflow means from a 1-year long vector of
% daily streamflow means. The output vector represents group 1 of IHA indicators:
% "magnitude of monthy water conditions".
% 
% Input arguments
%   q = vector of daily streamflow means [m^3/s] - vector [365 or 366x1]
%   
% Output arguments
%   IHA_g1 = vector of monthy streamflow means [m^3/s] - vector [12x1]
%
% Andrea Cominola, Emanuele Mason 28/09/2010

if (size(q,1) < 365) || (size(q,1) > 366);
    disp(' Error: q = vector of daily streamflow means [m3/s]; - vector [365 or 366x1]');
    return;
end

IHA_g1 = NaN(12,1);

% Length of months in a non leap year:
month_length = [31 28 31 30 31 30 31 31 30 31 30 31];

% Condition for a leap year:
if length(q) == 366;
    month_length(2) = 29; 
end;


IHA_g1(1) = nanmean( q( 1 : month_length(1) ,1 ));
for i = 2 : 12 ;
    IHA_g1(i) = nanmean( q( (nansum( month_length(1:i-1))+1) : nansum( month_length(1:i)) ,1 )); %Gorla added "+1", now it is correct.
end