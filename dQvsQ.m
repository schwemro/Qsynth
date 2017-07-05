function dQvsQ(Date, Q)
% Plot dQ/dt against Q. Requires Date and streamflow both as vector.
TT = timetable(Date, Q);
TT.Properties.VariableNames = {'Q'};
TT1 = lag(TT);
TT2 = synchronize(TT, TT1);
TT2.Properties.VariableNames = {'Q' 'Q_1'};
TT2.dQ = TT2.Q-TT2.Q_1;
TT3 = sortrows(TT2,'Q','ascend');

f1 = figure('Name','dQ vs Q','NumberTitle','off');
scatt_obser(TT3.Q,TT3.dQ,3,'k','filled');
xlabel('Q');
ylabel('dQ/dt');
saveas(f1,'dQ_vs_Q_dotty.fig');

f2 = figure('Name','dQ vs Q','NumberTitle','off');
plot(TT3.Q,TT3.dQ,'k');
xlabel('Q');
ylabel('dQ/dt');
saveas(f2,'dQ_vs_Q.fig');
end

