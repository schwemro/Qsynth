function [date, month, year] = time_JD2date( JD , format )
% Copyright (C) 2005,2006 Daniele de Rigo 
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
% [date, month, year] = time_JD2date( JD , format )
% 
% Input arguments:
%
%   JD      vector of dates in Julian Date format
%   format  string of the required output format
%           (explained by using the example:
%           JD=time_date2JD(31,12,2000) % JD = 2451910):
% 
%      time_JD2date( JD , 'DD/MM/YYYY' ) --->  ans = '31/12/2000'
%      time_JD2date( JD , 'YYYY-MM-DD' ) --->  ans = '2000-12-31'
%           if omitted, it is returned a numeric matrix:
%      time_JD2date( JD )                --->  ans = [ 31 12 2000 ]
% 
% Output arguments:
%
% if only 'date' was required as output, 'date' will have a format
% depending on the 'format' string, as explained before
% 
% else 'date' , 'month' , 'year' will have the format:
% 
% 'date'  format:  DD
% 'month' format:  MM
% 'year'  format:  YYYY
% 
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
%   [dd2 mm2 yyyy2] = time_JD2date( time_date2JD(dd,mm,yyyy) );
% 
%   disp('Are all the Julian Dates correctly rebuilt?')
%   if( all( all([dd2 mm2 yyyy2] == [dd mm yyyy]) ) ) disp(' YES '); else disp(' NO '); end;
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
% History
% 2006-06-16 version 0.2
%    improved the format rendering
%    extended the format choice
%
% New examples:
%   date = [21 4 -753; 1 1 2000; 25 11 -4713; 8 4 -563];
%   JD   = time_date2JD(date);
%   time_JD2date(JD)                             % old basic feature
%   time_JD2date(JD,  'DD/MM/YYYY'  )            % old formatted style (new:spacing corrected) 
%   time_JD2date(JD,  'YYYY-MM-DD'  )            % old formatted style (new:spacing corrected) 
%   time_JD2date(JD(2),  'YYYY-MM-DD'  )         % old formatted style (new:spacing corrected) 
%   time_JD2date(JD,  'YYYY,DDth MM'  )          % new formatted style  
%   time_JD2date(JD,  'year: YYYY  month:MM'  )  % new formatted style (you can omit some element) 
% 


if(nargin < 1)
  disp(  'Error: too few arguments.'  );
  disp(  'Usage: [date month year] = time_JD2date( JD, format )'  );
  return;
end;

if(nargin > 3)  
  disp(  'Error: too many arguments.'  ); 
  disp(  'Usage: [date month year] = time_JD2date( JD, format )'  );
  return;
end;

if( min(size(JD)) ~= 1)
  disp(  '''JD'' must be scalar or vector'  );
  return;
end;

if (nargin > 1 )
  if( ischar(format)==0 )
    disp(  'Error: the format argument must be a string.'  );
    return;
  end;
  if( size(format,1)~=1 )
    disp([  'Error: a vector of '  num2str(size(format,1))  ' strings was passed as format argument.'  ])
    disp(   '       The format argument must be a single string.'  );
    return;
  end;
  
  date_id    =  'DMY'  ;
  date_style =  {  '%0'   ,  '%0'     ,  '%'     };
  name_id    =  {  'day'  ,  'month'  ,  'year'  };
  
  valid_idx  = find(any(repmat(date_id,length(format),1)==repmat(format',1,length(date_id))));

  if ~length(valid_idx)
     disp([  'Error: the 2nd argument ''template'' ('''  format  ''')'  ]);
     disp([  '       must contain at least one of the char codes: '''  date_id  '''.'  ]);
     return;
  end;


  date_id    = date_id(valid_idx);
   
  tmp_replace = {  'date_style'  ,  'name_id'  };
  for i=1:length(tmp_replace);
     eval([  'tmp_out  = sprintf(  '','  tmp_replace{i}  '{%d}'' ,1:length(valid_idx));'  ]);
     tmp_out  = tmp_out(2:end);
     eval([  'tmp_in   = sprintf(  '',''''%s'''''' , '  tmp_replace{i}  '{valid_idx});'  ]);
     tmp_in   = tmp_in(2:end);
     eval(sprintf(  'clear %s; [%s]=deal(%s);'  , tmp_replace{i}, tmp_out , tmp_in ));
  end;
%  date_id
%  date_style
%  name_id
  
  %format(find(format=='d'))='D';
  %format(find(format=='m'))='M';
  %format(find(format=='y'))='Y';
  
  % Explanatory examples about 'format' patterns:
  %
  %  'DD/MM/YYYY'
  %   |  |  |
  %   |  |  +----------------------+
  %   |  +----------------------+  |
  %   +----------------------+  |  |
  %                          |  |  |
  %      token_start_pos = [ 1  4  7 ]
  %      token_length    = [ 2  2  4 ] 
  %      separators      = { '' , '/' , '/' , '' }
  %      sort_idx        = [ 1  2  3 ]
  % 
  %  'inflow_YYYY_MM_DD.txt'
  %          |    |  |
  %          +----|--|-------------+
  %               +--|----------+  |
  %                  +-------+  |  |
  %                          |  |  |
  %      token_start_pos = [16 13  8 ]
  %      token_length    = [ 2  2  4 ] 
  %      separators      = { 'inflow_' , '_' , '_' , '.txt' }
  %      sort_idx        = [ 3  2  1 ]
  % 
  
  token_start_pos = zeros(size(date_id));
  token_length    = zeros(size(date_id));
  separators{length(date_id)+1}  =  ''  ;   
  
  date_id = upper(date_id);  % prevents errors in date_id initialization
  for i = 1:length(date_id)
      token_pos = find(format==date_id(i));
      if length(token_pos)
        if all(diff(token_pos)==1)  % if all the i-th charcodes are contiguous 
          token_start_pos(i) = token_pos(1);
          token_length(i)    = length(token_pos);
        else
          disp(   'Error: wrong format argument.'  );
          disp([  '       The '''  date_id(i)  ''' codes ('  name_id{i}  ') in the ''format'' argument: '''  format  ''' must be contiguous.'  ]);
          return;
        end;
      end;
  end;
  
%  token_start_pos
%  token_length
%  date_id
%  valid_id=find(token_length)
  
  [sort_start_pos, sort_idx ]=sort(token_start_pos);
  bordered_start_pos = [1 sort_start_pos+token_length(sort_idx)];
  bordered_end_pos = [sort_start_pos-1 length(format)];
  
  
  [ day_value , month_value , year_value ] = get_date_from_JD( JD );
  %         rough_date   = [day_value   month_value  year_value ];
  eval( [  'rough_date   = ['  sprintf(  '%s_value '  ,name_id{:})  '];'  ]);
      
  %  max( floor(log10(abs(rough_date)))+1  + (sign(rough_date)<0) ,[],1)
  %  \    \_____________________________/    \_________________/      /
  %   \      number of  decimal digits         plus 1 only if        /
  %    \     for each day, month, year         negative value       /
  %     \__________________________________________________________/
  %        max number of decimal digits of day, month, year
  valid_token_length = max(  token_length , max(floor(log10(abs(rough_date)))+1+(sign(rough_date)<0),[],1)  );
  
  pattern =  ''  ;
%  for i = 1:(length(date_id)+1)
%      separators{i} = format(bordered_start_pos(i):bordered_end_pos(i));
%      pattern = [  pattern separators{i}  ];
%      if i<=length(date_id)
%         pattern = [ pattern  '%0'  num2str(token_length(sort_idx(i)))  'd'   ];
%      end;
%  end;
  for i = 1:(length(date_id)+1)
      separators{i} = format(bordered_start_pos(i):bordered_end_pos(i));
      pattern = [  pattern separators{i}  ];
      if i<=length(date_id)
         pattern = [ pattern  date_style{sort_idx(i)}  num2str(valid_token_length(sort_idx(i)))  'd'   ];
      end;
  end;
%  pattern
    
else
  format = '';
  [ day_value , month_value , year_value ] = get_date_from_JD( JD );
end;




if(nargout <= 1)
  if(nargin >1)
    str_date=sprintf(pattern,rough_date(:,sort_idx)');
    date = reshape( str_date , length(str_date)/length(JD) , length(JD) )';
  else
    date = [   day_value   month_value   year_value  ];
  end;
%  switch format
%    case      ''  
%       date = [   day_value   month_value   year_value  ];
%    case      'DD/MM/YYYY'  
%       date= reshape(sprintf(  '%02d/%02d/%04d'  ,[ day_value   month_value  year_value ]'),10,length(JD))';
%    case      'YYYY-MM-DD'  
%       date= reshape(sprintf(  '%04d-%02d-%02d'  ,[ year_value  month_value  day_value  ]'),10,length(JD))';
%    otherwise
%       date = [   day_value   month_value   year_value  ];
%  end;  
else
  date  = day_value  ;
  month = month_value;
  year  = year_value ;
end;


function [ day_value , month_value , year_value ] = get_date_from_JD( JD )
    j = JD - 1721119;

    year_value  = floor(((j*4) - 1) / 146097);
    j           = (j*4) - 1 - 146097*year_value;
    d           = floor(j/4);
    j           = floor( ((d*4) + 3) / 1461 );
    d           = (d*4) + 3 - 1461*j;
    d           = floor((d + 4)/4);
    month_value = floor( (5*d - 3)/153 );
    d           = 5*d - 3 - 153*month_value;
    day_value   = floor((d + 5)/5);
    year_value  = (100*year_value + j);
    
    idx=1:length(month_value);
    idx_greaterequal_to10=idx(month_value >= 10);
    idx_less_than10=idx(month_value < 10);
    month_value(  idx_less_than10      ) = (month_value( idx_less_than10      ) + 3);
    month_value(  idx_greaterequal_to10) = (month_value( idx_greaterequal_to10) - 9);
    year_value(idx_greaterequal_to10) = year_value(idx_greaterequal_to10) + 1;


