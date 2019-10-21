%% Setup
mini = 0;
maxi = 3;
framesnum = 160;

%% Control group
[dorsal_ratio, ventral_ratio, dorsal, ventral] = collect_DV_data(pwd, framesnum);
close all;
im_DV_data(dorsal_ratio, ventral_ratio, [mini maxi], [mini maxi], [1 0]);
%% Save control group data
save('Ctrl_Muscle_Dorsal_Ventral.mat', 'dorsal', 'dorsal_ratio', 'ventral', 'ventral_ratio');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ATR group
[dorsal_ratio, ventral_ratio, dorsal, ventral] = collect_DV_data(pwd, framesnum);
close all;
im_DV_data(dorsal_ratio, ventral_ratio, [mini maxi], [mini maxi], [0 1]);
%% Save control group data
save('ATR_Muscle_Dorsal_Ventral.mat', 'dorsal', 'dorsal_ratio', 'ventral', 'ventral_ratio');