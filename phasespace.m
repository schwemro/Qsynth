function phasespace(tiseanPath, Q)
% Phase space reconstruction of runoff time series Q.
acf = autocorr(Q, 365);
tau0 = find(acf<0,1)-1;
tau01 = find(acf<0.1,1)-1;
tau05 = find(acf<0.5,1)-1;
save tisean_in.txt Q -ascii
system([tiseanPath,'mutual -D365 -o mutual_out.txt tisean_in.txt']);
mutual_out = dlmread('mutual_out.txt', ' ', 1, 0);
[val, taumut] = min(mutual_out(:,2));
system([tiseanPath,'delay -d' num2str(tau0) ' -o delay_tau0.txt tisean_in_obs.txt']);
system([tiseanPath,'delay -d' num2str(tau01) ' -o delay_tau01.txt tisean_in_obs.txt']);
system([tiseanPath,'delay -d' num2str(tau05) ' -o delay_tau05.txt tisean_in_obs.txt']);
system([tiseanPath,'delay -d' num2str(taumut) ' -o delay_taumut.txt tisean_in_obs.txt']);
load delay_tau0.txt;load delay_tau01.txt;load delay_tau05.txt;load delay_taumut.txt;

f1 = figure('Name','Phase Space','NumberTitle','off');
subplot(2,2,1)
scatt_obser(delay_tau0(:,1),delay_tau0(:,2),3,'b','filled');
xlabel('Q(t-{\tau})');
ylabel('Q(t)');
title(['ACF: {\tau} = ' num2str(tau0)])
subplot(2,2,2)
scatt_obser(delay_tau01(:,1),delay_tau01(:,2),3,'b','filled');
xlabel('Q(t-{\tau})');
ylabel('Q(t)');
title(['ACF: {\tau} = ' num2str(tau01)])
subplot(2,2,3)
scatt_obser(delay_tau05(:,1),delay_tau05(:,2),3,'b','filled');
xlabel('Q(t-{\tau})');
ylabel('Q(t)');
title(['ACF: {\tau} = ' num2str(tau05)])
subplot(2,2,4)
scatt_obser(delay_taumut(:,1),delay_taumut(:,2),3,'b','filled');
xlabel('Q(t-{\tau})');
ylabel('Q(t)');
title(['MI: {\tau} = ' num2str(taumut)])
saveas(f1,'phasespace_dotty.fig');

f2 = figure('Name','Phase Space Q_obs','NumberTitle','off');
subplot(2,2,1)
plot(delay_tau0(:,1),delay_tau0(:,2),'b');
xlabel('Q(t-{\tau})');
ylabel('Q(t)');
title(['ACF: {\tau} = ' num2str(tau0)])
subplot(2,2,2)
plot(delay_tau01(:,1),delay_tau01(:,2),'b');
xlabel('Q(t-{\tau})');
ylabel('Q(t)');
title(['ACF: {\tau} = ' num2str(tau01)])
subplot(2,2,3)
plot(delay_tau05(:,1),delay_tau05(:,2),'b');
xlabel('Q(t-{\tau})');
ylabel('Q(t)');
title(['ACF: {\tau} = ' num2str(tau05)])
subplot(2,2,4)
plot(delay_taumut(:,1),delay_taumut(:,2),'b');
xlabel('Q(t-{\tau})');
ylabel('Q(t)');
title(['MI: {\tau} = ' num2str(taumut)])
saveas(f2,'phasespace.fig');

delete delay_tau0.txt;delete delay_tau01.txt;delete delay_tau05.txt;delete delay_taumut.txt;
delete tisean_in.txt;delete mutual_out.txt;
end