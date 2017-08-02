function date_num = time_date2JD( date , month , year)
% Copyright (C) 2005 Daniele de Rigo 
% 
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
% 
% ---------------------------------------------------------------------------
% 
% date_num = time_date2JD( date , month , year)
% 
% if only 'date' was passed, 'date' must have one of these formats :
% 
%   'date' format:  'YYYY-MM-DD'
%   'date' format:  'DD/MM/YYYY'
%   'date' format:  [DD MM YYYY]
% 
% else 'date' , 'month' , 'year' must have the format:
% 
% 'date'  format:  DD
% 'month' format:  MM
% 'year'  format:  YYYY
% 
% This algorithm is based on Algorithm 199, CACM, Aug 1963. The
% algorithm is not valid for dates before Sep 14, 1752.
% I.e. these dates are not guaranted to be contiguous, according
% to the gregorian and julian calendar transitions.
% Therefore this algorithm has been tested to guarantee the
% Julian Date (JD) contiguity for contiguous days in the range
% [ 1-st januar 1000 B.C. , 31-th december 4000 A.C. ]:
% 
%   time_date2JD( next_day('DD/MM/YYYY') ) -
%   time_date2JD(          'DD/MM/YYYY'  ) == 1
% 
% You can prove that by executing the code:
%  
%   y_end=4000;  y_begin=-1000; 
%   n=365*(y_end-y_begin+1)+sum(is_leap_year(y_begin:y_end));
% 
%   dd=zeros(n,1);mm=zeros(n,1);yyyy=zeros(n,1);
%   md=[31 28 31 30 31 30 31 31 30 31 30 31];
%   cd=0; 
%   for y=y_begin:y_end mmy=md; mmy(2)=mmy(2)+is_leap_year(y); 
%     for m=1:12 
%       idx=[1:mmy(m)]+cd; 
%       dd(idx)=[1:mmy(m)]; 
%       mm(idx)=repmat(m,mmy(m),1); 
%       yyyy(idx)=repmat(y,mmy(m),1); 
%       cd=cd+mmy(m); 
%     end; 
%   end;
% 
%   JD = time_date2JD(dd,mm,yyyy);
%   
%   disp('Are all the Julian Dates contiguous?')
%   if( all(diff(data)==1) ) disp(' YES '); else disp(' NO '); end;
% 
% 
% 
% 
% Some other details can be read on:
%  http://en.wikipedia.org/wiki/Julian_calendar
% 
% << The Julian calendar was in general use in Europe from the 
% times of the Roman Empire until 1582, when Pope Gregory XIII 
% promulgated the Gregorian Calendar, which was soon adopted 
% by most Catholic countries. The Protestant countries followed 
% later, and the countries of Eastern Europe even later. Great 
% Britain had Thursday 14 September 1752 follow 
% Wednesday 2 September 1752. Sweden adopted the new style calendar 
% in 1753, but also for a twelve-year period starting in 1700 
% used a modified Julian Calendar. Russia remained on the Julian 
% calendar until after the Russian Revolution (which is thus 
% called the 'October Revolution' but occurred in November 
% according to the Gregorian calendar). >>
% 
% and from: http://en.wikipedia.org/wiki/Gregorian_Calendar
% 
% << When the new calendar was put in use, to correct the error  
% already accumulated in the thirteen centuries since the Council  
% of Nicaea, a deletion of ten days was made in the solar calendar.  
% The last day of the Julian calendar was October 4, 1582 and this  
% was followed by the first day of the Gregorian calendar  
% October 15, 1582. Nevertheless, the dates "5 October 1582" to  
% "14 October 1582" (inclusive) are still valid in virtually all  
% countries because even most Roman Catholic countries did not adopt  
% the new calendar on the date specified by the bull, but months or  
% even years later (the last in 1587). >>
% 



if( (nargin < 1)|(nargin > 3) ) 
  disp(  'Usage: date_num = time_date2JD( date , month , year)'  );
  return;
end;

if( nargin == 1)

  n_date=size(date,1);
  if( ( size(date,2)~=10 ) & (size(date,2)~=3 ) )
    disp(  'Error: if only ''date'' was passed, ''date'' must have one of these formats :'  );
    disp(  '               ''date'' format:  ''YYYY-MM-DD'''  );
    disp(  '               ''date'' format:  ''DD/MM/YYYY'''  );
    disp(  '               ''date'' format:  [DD MM YYYY]'  );
    return;
  end;
  if( size(date,2)==10 )
     if( all( (date(:,3)=='/') & (date(:,6)=='/') ) )
       dd   = str2num(date(:,1:2));
       mm   = str2num(date(:,4:5));
       yyyy = str2num(date(:,7:10));
     elseif( all( (date(:,5)=='-') & (date(:,8)=='-') ) )
       dd   = str2num(date(:,9:10));
       mm   = str2num(date(:,6:7));
       yyyy = str2num(date(:,1:4));
     else
       disp(  'Error: if only ''date'' was passed, ''date'' must have one of these formats :'  );
       disp(  '               ''date'' format:  YYYY-MM-DD'  );
       disp(  '               ''date'' format:  DD/MM/YYYY'  );
       disp(  '               ''date'' format:  [DD MM YYYY]'  );
       return;
     end;
     if(min([size(dd,1),size(mm,1),size(yyyy,1)])~=n_date)
       disp(  'Error: if only ''date'' was passed, ''date'' must have one of these formats :'  );
       disp(  '               ''date'' format:  YYYY-MM-DD'  );
       disp(  '               ''date'' format:  DD/MM/YYYY'  );
       disp(  '               ''date'' format:  [DD MM YYYY]'  );
       return;
     end;
  else
    dd   = date(:,1);
    mm   = date(:,2);
    yyyy = date(:,3);
  end;
  
elseif( nargin == 3)
  if( ( min(size(date))~=1 ) | ( min(size(month))~=1 ) | ( min(size(year))~=1 ) )
    disp(  'Error: ''date'',''month'',''year'' must be scalars or vectors'  );
    return;
  end;
  
  date  = reshape( date  , length( date  ) , 1);
  month = reshape( month , length( month ) , 1);
  year  = reshape( year  , length( year  ) , 1);
   
  n_date= max( [length(date) length(month) length(year)] );
  if( (n_date>1) & (length(date)==1) )
     date  = ones(n_date,1)*date ;
  end;
  if( (n_date>1) & (length(month)==1) )
     month = ones(n_date,1)*month;
  end;
  if( (n_date>1) & (length(year)==1) )
     year  = ones(n_date,1)*year ;
  end;
  
  if( ( size(date,1)~=size(month,1) ) | ( size(date,1)~=size(year,1) ) )
    disp(  'Error: ''date'',''month'',''year'' must be scalars or vectors with the same number of rows'  );
    return;
  end;
  
  
  if( (max(month)>12) | (min(month)<1) )
    disp(  'Error: ''month'' must be in [1,12]'  );
    return;
  end;
  dd   = date;
  mm   = month;
  yyyy = year;
else
  disp(  'Error: ''year'' required'  );
  return;
end;

if( any(round(mm)~=mm ) )
  disp(  'Error: month must be integer'  );
  return;
end;
if( any(round(yyyy)~=yyyy ) )
  disp(  'Error: year must be integer'  );
  return;
end;
month=mm;
%      if (yyyy <= 99)
%          yyyy += 1900;
%      end;
%      if (dd <= 0 || dd > Date::DaysInMonth (s.m, s.y))
%          return -1;

    idx=1:length(mm);
    idx_greater_than2=idx(mm > 2);
    idx_lessequal_to2=idx(mm <= 2);
    
    mm(  idx_greater_than2) = (mm( idx_greater_than2) - 3);
    mm(  idx_lessequal_to2) = (mm( idx_lessequal_to2) + 9);
    yyyy(idx_lessequal_to2) = yyyy(idx_lessequal_to2)-1;
    
    
    c  = floor(yyyy / 100);
    yy = floor(yyyy - 100*c);
    date_num = floor((146097 .* c)./4) + floor((1461*yy)./4) + floor((153*mm + 2)/5) + dd + 1721119;

   % date_num2 = 367*yyyy -floor((7*floor(yyyy+(month+9)/12))/4) +floor(275*month/9) +dd +1721014;
                           %floor((7*floor(1993+  13  /12))/4)
   % date_num = [date_num date_num2];


