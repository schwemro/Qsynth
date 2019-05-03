function dQ = dQdt(Date, Q)
% Computes dQ/dt. Input arguments are date and streamflow both. Output
% argument is derivation of streamflow.
TT = timetable(Date, Q);
TT.Properties.VariableNames = {'Q'};
TT1 = lag(TT);
TT2 = synchronize(TT, TT1);
TT2.Properties.VariableNames = {'Q' 'Q_1'};
dQ = TT2.Q-TT2.Q_1;
end

