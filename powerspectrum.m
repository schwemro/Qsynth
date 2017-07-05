function powerspectrum(Q)
% Plotting the power spectrum as a function of frequency.
% Requires runoff time series with daily values.
N = length(Q);
fts = fft(Q);
fs = 1/86400;                                   % sample frequency
power = abs(fts(1:floor(N/2))).^2;          % power of first half of transform data
maxfreq = 1/2;                                  % maximum frequency
freq = (1:N/2)/(N/2)*maxfreq;           % equally spaced frequency grid

f1 = figure('Name','Power Spectrum','NumberTitle','off');
plot(freq(2:end), power(2:end));
grid;
xlabel('1/Day');
ylabel('Power');
saveas(f1,'Power_spectrum.fig')
end
