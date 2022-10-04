%% Clear everything

clear; close all;

%% Load curated curvature and activity data, and normalize signal

% Load curated curvature and normalize signal

f_curvature = uigetfile('*_curvature_curated.mat', ...
    'Select curated curvature');
if f_curvature == 0
    fprintf('user cancelled selection. \n');
else
    load(f_curvature);
    curv_norm = normalize_signal(curvdatafiltered); % This is already filtered
    fprintf(['curvature data for ' f_curvature ' loaded and normalized. \n']);
end

% Load activity and normalize signal

c = regexp(f_curvature, '\d*');
f_activity = uigetfile('*.mat', ...
    ['Select activity for ' f_curvature], [f_curvature(1:c(1)) '*.mat']);
if f_activity == 0
    fprintf('user cancelled selection. \n');
else
    load(f_activity);
    neuronactivity = signal{1,1}./signal_mirror{1,1};
    neuronactivity_smd = smoothdata(neuronactivity, 'rloess');
    neuronactivity_norm = normalize_signal(neuronactivity_smd);
    fprintf(['activity data for ' f_activity ' loaded and normalized. \n']);
end

%% Draw figure

figure;

yyaxis right;
plot(curv_norm, 'linewidth', 4, ...
    'color', [0 0.5 1]);
ylim([-0.1 1.1]);

yyaxis left;
plot(neuronactivity_norm, 'linewidth', 4, ...
    'color', [1 0.5 0]);

% axis equal;
ylim([-0.1 1.1]);
set(gca, 'visible', 'off');

%% Save figure

savefig([f_activity(1:end-4) '.fig']);

%% Calculate cross-correlation

[r, lag] = xcorr(neuronactivity_norm, curv_norm, 'normalized');
plot(lag, r, 'k', 'linewidth', 4);
