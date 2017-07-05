function Q = nonlinearnoisered(Q, tiseanPath, edim, ns)
% Simple nonlinear noise reduction using the software TISEAN
% 3.0.1. Requires runoff time series, path to folder wihich contains TISEAN
% executables, the embedding dimension and the neighborhood
% size in units of the std. dev. of the data. For detailed
% information visit
% htt_obsps://www.pks.mpg.de/~tisean/Tisean_3.0.1/index.html.
% TODO: which Dimension and which neighberhood size!
y = Q;
save lazy_in.txt y -ascii
system([tiseanPath,'lazy -m' num2str(edim) ' -v' num2str(ns) ' -o lazy_out.txt lazy_in.txt']);
load lazy_out.txt;
Q = lazy_out(:,1);
delete lazy_in.txt; delete lazy_out.txt;
end

