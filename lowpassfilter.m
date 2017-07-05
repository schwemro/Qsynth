function Q = lowpassfilter(Q, cf)
% Applying a low-pass filter to runoff time series.
% Requires time series as vector and cutoff frequency as input arguments.
% Returns filtered time series.
N = length(Q);
fts = fft(Q);
maxfreq = 1/2;
freq = (1:N/2)/(N/2)*maxfreq;
cf = sum(freq<cf);
Q_sim_fts = zeros(size(fts));
Q_sim_fts(1:cf+1) = 1;
Q_sim_fts(end-cf+1:end) = 1;
Q = ifft(fts.*Q_sim_fts);
end

