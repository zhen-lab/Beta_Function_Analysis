%% Determine name for the file

parts = strsplit(pwd, '\');
partlast = parts{1,end};
fname = [partlast '_Dorsal_Ventral_All'];

%% Collect raw data

framenum = 60;
[dorsal_ratio, ventral_ratio, dorsal, dorsal_gfp, ventral, ventral_gfp] ...
    = collect_DV_data(pwd, framenum);
save(fname, 'dorsal_ratio', 'ventral_ratio', ...
    'dorsal', 'dorsal_gfp', 'ventral', 'ventral_gfp');

%% Collect mean data

[dorsal_mean, ventral_mean] = collect_DV_data_mean(dorsal, ventral);
dv_mean = [dorsal_mean ventral_mean];
save([fname '_mean'], 'dv_mean');

%% Clear everything

clear; close all;