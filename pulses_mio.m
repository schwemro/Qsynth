function [IHA_g4, mem] = pulses_mio(q, quant_25, quant_75)
% ::: Help :::
%
% [IHA_g4] = pulses(q, quant_25, quant_75)
% 
% Evaluates the four intra-annual indicators in the 4th group of IHA
% statistics. The group is called "frequency and duration of high and low
% pulses". An hydrological pulse is defined in the IHA paper as those
% periods within a year in which the daily mean water condition either
% rises above the 75th percentile (high pulse) or drops below the 25th
% percentile (low pulse) of all daily streamflow values for the pre-impact
% time period. The function evaluates number and mean duration of both
% high and low pulses within a year of records, given the parameters.
% 
% Input argument:
%   q = vector of daily streamflow means [m^3/s]          - vector [365 or 366x1]
%   quant_25 = the 25th percentile of all daily streamflow values for the
%              pre-impact time period.                    - double number
%   quant_75 = the 75th percentile of all daily streamflow values for the
%              pre-impact time period.                    - double number
%   
% Output argument:
%   IHA_g4 = structure with four IHA-indicators that belong to group 4.
%   It contains four fields:
%    IHA_g4.mean_durat_high_p = the mean duration of high pulses [day]
%                                                      - double number
%    IHA_g4.mean_durat_low_p = the mean duration of low pulses [day]
%                                                      - double number
%    IHA_g4.n_high_p = the number of high pulses in the daily streamflow
%       serie [.]                                      - double number
%    IHA_g4.n_low_p = the number of low pulses in the daily streamflow
%       serie [.]                                      - double number
% 
%    mem.low=[num low periods; bool low cell 1; durat low event 1; bool low
%    last cell; durat low last event]
%    mem.high=[num high periods; bool high cell 1; durat high event 1; bool
%    high last cell; durat high last event]
% 
% Andrea Cominola, Emanuele Mason, Lorenzo Gorla 06-01-2012

if( (nargin < 1)||(nargin > 3) ) 
  disp(  'Usage: IHA_g4 = pulses(q, quant_25, quant_75)'  );
  return;
end
if (size(q,1) < 365) || (size(q,1) > 366);
    disp(' Error: q = vector of daily streamflow means [m3/s]; - vector [365 or 366x1]');
    return;
end
if( (max(size(quant_25)) ~= 1 )||( max(size(quant_75)) ~= 1) ) 
  disp(  'Error: ''quant_25'' and '' quant_75'' must be scalars'  );
  return;
end
if quant_25 > quant_75
  disp(  'Error: Value of 25th percentile greater than the 75th one'  );
  return;
end

%Commento le due righe seguenti se voglio calcolare gli indicatori con
%percentili fissi, come da software originale; le calcolo se voglio invece
%calcolare i percentili anno per anno.
% quant_25=prctile(q,25);
% quant_75=prctile(q,75);

bool_high_p = q > quant_75;
bool_low_p = q < quant_25;

% Number of low and high pulses within the serie.
IHA_g4.n_high_p = sum(diff([0; bool_high_p; 0]) == 1);
IHA_g4.n_low_p = sum(diff([0; bool_low_p; 0]) == 1);

% Mean duration of the pulses.
high_p_start = find(diff([0; bool_high_p; 0])==1);
low_p_start = find(diff([0; bool_low_p; 0])==1);
high_p_end = find(diff([0; bool_high_p; 0])==-1);
low_p_end = find(diff([0; bool_low_p; 0])==-1);

if size(high_p_start)==[0 1]; durat_high_events=0; else durat_high_events=high_p_end - high_p_start; end
if size(low_p_start)==[0 1]; durat_low_events=0; else durat_low_events=low_p_end - low_p_start; end
IHA_g4.sum_durat_high_p = sum(durat_high_events);%sum al posto di mean
IHA_g4.sum_durat_low_p = sum(durat_low_events);%sum al posto di mean

if IHA_g4.n_high_p == 0; IHA_g4.sum_durat_high_p = 0; end;%sum al posto di mean
if IHA_g4.n_low_p == 0; IHA_g4.sum_durat_low_p = 0; end;%sum al posto di mean

mem.low=NaN(5,1);
mem.low(1,1)=IHA_g4.n_low_p;
mem.low(2,1)=bool_low_p(1,1);
mem.low(3,1)=durat_low_events(1,1);
mem.low(4,1)=bool_low_p(end,1);

mem.high=NaN(5,1);
mem.high(1,1)=IHA_g4.n_high_p;
mem.high(2,1)=bool_high_p(1,1);
mem.high(3,1)=durat_high_events(1,1);
mem.high(4,1)=bool_high_p(end,1);

end






