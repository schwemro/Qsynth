function mm = moving_average( y ,f )
% :::HELP:::
%
% mm = media_mobile( y ,f )
%
% The function calculates the mobile average for the vector y, with the
% interval f.
%
% Input arguments:
%       y: vector of values         -vector [Nx1]
%       f: length of the interval to use for calculating the mobile average
%                                   -integer number
%      
% Output arguments:
%       mm: vector of the calculated values of mobile average for the vector y    -vector [N-f+1x1] 
%
% The element in position (i) in vector mm is the mobile average calculated
% taking the elements in position (i:(i+f-1)) in vector y.     
%
% Andrea Cominola, Emanuele Mason 28/09/2010

mm = zeros( (length(y)-f+1), 1 );
for k = 1 : (length(y)-f+1)
    mm( k ) =  mean( y( k : (k+ f-1) ) )  ;
end
end