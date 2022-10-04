%% Clear everything

clear; close all;

%% Collect neuron activity data

parts = strsplit(pwd, '\');
fname = [parts{end} '_All'];
framenum = 60;

%% For new data

[GCaMP, RFP, ratio, ratio_norm] ...
    = collect_neurons_data(pwd, framenum);
save(fname, 'GCaMP', 'RFP', 'ratio', 'ratio_norm');

%% For old data

ratio_start = repmat(ratio(1, :), framenum, 1);
ratio_norm = (ratio./ratio_start)';
save(fname, 'ratio_norm');
