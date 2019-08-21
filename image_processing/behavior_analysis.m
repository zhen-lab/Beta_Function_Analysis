%% Segment the images with behavior recordings

behavior_segmentation;
fprintf('segmentation of recording is done. \n');

%% Load tif file after segmentation is completed

setup_proof_reading;

%% Calculate speed by correcting stage movement

[vel_anterior_sign, vel_posterior_sign, ...
    dorsal_data, ventral_data, centerline_data_spline] = ...
    calculate_speed_behavior(filename, pathname, data);
fprintf('speed calculation is done. \n');

%% Curate speed profile by GUI

smallarea = 1; % If the images were cropped in segmetation into 1/4, set this to 1; else set this to 0
behavior_correct_speed(...
    imagelist, filename, 1, 2000, ...
    dorsal_data, ventral_data, ...
    centerline_data_spline, ...
    vel_anterior_sign, vel_posterior_sign,...
    smallarea);

%% Draw figures

totalfrm = size(imagelist,1);
windowlength = 9;

% brightness = cellfun(@(x) mean(x,'all'), imagelist(:,1), 'UniformOutput', false);
brightness = zeros(totalfrm,1);
for i = 1:totalfrm
    brightness(i) = mean(imagelist{i,1}, 'all');
end

velocity_corrected_mean = ...
    mean([velocity_corrected_anterior velocity_corrected_posterior], 2);
velocity_corrected_mean_smd = smoothdata(velocity_corrected_mean, 'movmedian', windowlength);

% figure; hold on;
% plot(vel_anterior_sign, 'g');
% plot(vel_posterior_sign, 'k');
% plot(velocity_corrected_mean, 'm');
% plot(brightness);
% plot([1 totalfrm], [0 0], ':k');

s = figure; hold on;
plot([1 totalfrm], [0 0], 'k', 'linewidth', 1);
plot(brightness-min(brightness), 'g', 'linewidth', 1);
plot(velocity_corrected_mean_smd, 'k', 'linewidth', 2);

%% Save curated speed profile

save([filename(1:end-4) '_velocity_curated.mat'], ...
    'velocity_corrected_anterior', 'velocity_corrected_posterior',...
    'velocity_corrected_mean', 'velocity_corrected_mean_smd',...
    'centerline_corrected');
savefig(s, [filename(1:end-4) '_velocity_curated.fig']);
saveas(s, [filename(1:end-4) '_velocity_curated.tif'], 'tiffn');

fprintf(['curated speed data and figures saved for ' filename '.\n']);

%% Clear everything

clear; close all

%% Obtain stimulation epochs and pool data from on/off stimulations
% after analysis is done

%% Extract optogenetics data and pool

% Select tiff file
setup_proof_reading;

% Generate brightness data
totalfrm = size(imagelist,1);
brightness = zeros(totalfrm,1);
for i = 1:totalfrm
    brightness(i) = mean(imagelist{i,1}, 'all');
end
brightnessSubtract = brightness - min(brightness);
plot(brightnessSubtract);

fprintf('brightness is calculated. \n');

%% Caculate stimuation epochs

% Set parameters
prompt = {'Median filter size', ...
    'Tolerance for boundaries', ...
    'Mid percent level (higher->taller)', ...
    'State levels (0 or 1)', ...
    'Number of stimutations'};
dlgtitle = ['Parameter for stimulations of ' filename];
dims = [1 35];
definput = {'10', '20', '50', '0', '4'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
medfiltsize = str2double(answer{1});
tolerancebrightness = str2double(answer{2});
midpercentreflevel = str2double(answer{3});
statelevel = str2double(answer{4});
numstim = str2double(answer{5});

% Median filter to remove jitters
brightnessMedFilt = medfilt1(brightnessSubtract, ...
    medfiltsize, 'truncate');
% plot(brightnessMedFilt);

close all;
% Calculate crossess
if statelevel==0
    
    [W, INITCROSS, FINALCROSS] = ...
        pulsewidth(brightnessMedFilt, ...
        'tolerance', tolerancebrightness, ...
        'MidPercentReferenceLevel', midpercentreflevel);
    if size(W,1)==numstim
        fprintf('pulse segmentation is correct! \n');
    else
        fprintf('the filter size and/or tolerance needs adjustment. \n');
    end
    
    % Plot result
    pulsewidth(brightnessMedFilt, ...
        'tolerance', tolerancebrightness, ...
        'MidPercentReferenceLevel', midpercentreflevel);
    
else % If the state levels are the reason for incorrect pulse segmentation
    
    prompt = {'High state level', 'Low state level'};
    dlgtitle = ['State levels for ' filename];
    dims = [1 35];
    definput = {'330', '10'};
    levels = inputdlg(prompt,dlgtitle,dims,definput);
    highlevel = str2double(levels{1});
    lowlevel = str2double(levels{2});
    [W, INITCROSS, FINALCROSS] = ...
        pulsewidth(brightnessMedFilt, ...
        'tolerance', tolerancebrightness, ...
        'MidPercentReferenceLevel', midpercentreflevel, ...
        'StateLevels', [lowlevel highlevel]);
    if size(W,1)==numstim
        fprintf('pulse segmentation is correct! \n');
    else
        fprintf('the filter size, tolerance, and/or boundaries needs adjustment. \n');
    end
    
    % Plot result
    pulsewidth(brightnessMedFilt, ...
        'tolerance', tolerancebrightness, ...
        'MidPercentReferenceLevel', midpercentreflevel, ...
        'StateLevels', [lowlevel highlevel]);
    
end

save([filename '_velocity_crosses.mat'], ...
    'W', 'INITCROSS', 'FINALCROSS');

%% Pool data

% Select the associated curated velocity data mat file
velfile = uigetfile('*curated.mat', ['Select curated velocity profile for ' filename]);
load(velfile, 'velocity_corrected_posterior');

% Smoothen velocity
hampelfiltsize = 27;
velHampel = ...
    hampel(velocity_corrected_posterior, hampelfiltsize);

% Pool velocity data
% Minus 1 to match the velocity data shift (velocity data matches frame plus 1)
init = [ceil(INITCROSS)-1; totalfrm]; 
final = [0; floor(FINALCROSS)-1];
stim = [init(1:end-1) final(2:end)];
nonstim = [final+1 init-1];
velStim = cell(size(stim,1),1);
velNonstim = cell(size(nonstim,1),1);
frmAll = 1:totalfrm-1; % Match the velocity data
frmStim = [];

for i = 1:numstim
    velStim{i} = velHampel(stim(i,1):stim(i,2));
    velNonstim{i} = velHampel(nonstim(i,1):nonstim(i,2));   
    if i==numstim
        velNonstim{i+1} = velHampel(nonstim(i+1,1):nonstim(i+1,2));
    end   
    frmStim = cat(2, frmStim, stim(i,1):stim(i,2));    
end

% Match the velocity data point number
frmNonstim = setdiff(frmAll, frmStim);
velStimCat = velHampel(frmStim); velStimCatMean = nanmean(velStimCat);
velNonstimCat = velHampel(frmNonstim); velNonstimCatMean = nanmean(velNonstimCat);

fprintf('velocity has been pooled. \n');

% Draw figures

close all;

% Velocity and stimulations
s(1) = figure(1);
hold on;
plot(brightnessSubtract, 'g');
plot(velHampel, 'k');

% Velocity and smoothened stimulations
set(gca, 'xticklabel', []);
s(2) = figure(2);
hold on;
plot(brightnessMedFilt, 'g');
plot(velHampel, 'k');
set(gca, 'xticklabel', []);

% Pooled and aligned velocities
s(3) = figure(3);
subplot(121); hold on;
plot(velStimCat, 'g'); 
plot(velNonstimCat, 'k');
subplot(122); hold on;
cellfun(@(x) plot(x, 'g'), velStim);
cellfun(@(x) plot(x, 'k'), velNonstim);

fprintf(['stimulations average velocity is: ' ...
    num2str(velStimCatMean) '\n']);
fprintf(['non-stimulations average velocity is: ' ...
    num2str(velNonstimCatMean) '\n']);

% Save data and figures

save([filename '_velocity_pooled.mat'], ...
    'velStim', 'velNonstim', ...
    'velStimCat', 'velNonstimCat', ...
    'velStimCatMean', 'velNonstimCatMean');

savefig(s, [filename '_velocity_stim_pooled_aligned.fig']);
saveas(s(1), [filename '_velocity_stimulations.tif'], 'tiffn');
saveas(s(2), [filename '_velocity_stimulations_filt.tif'], 'tiffn');
saveas(s(3), [filename '_velocity_pooled_aligned.tif'], 'tiffn');

fprintf(['data and figures have been saved for '...
    filename '. \n']);

%% Draw marked figure

lim = [-500 500];
s = figure; 
hold on; 
for i = 1:numstim
    fill([stim(i,1) stim(i,2) stim(i,2) stim(i,1)], [lim(1) lim(1) lim(2) lim(2)], [0 1 0], 'edgecolor', 'none');
end
plot(velHampel, 'k');
plot([1 length(velocity_corrected_posterior)],[0 0],'k');
set(gca, 'visible', 'off')

savefig(s, [filename '_velocity_stimulations_filt_marked.fig']);
saveas(s, [filename '_velocity_stimulations_filt_marked.tif'], 'tiffn');

%% Clear everything

clear; close all;

%% Collect data

f = dir(fullfile(pwd, '*pooled.mat'));
numsamples = numel(f);
velStimTotal = zeros(numsamples,1);
velNonstimTotal = zeros(numsamples,1);

for i = 1:numsamples
    pathname = fullfile(pwd, f(i).name);
    load(pathname, 'velStimCatMean', 'velNonstimCatMean');
    velStimTotal(i) = velStimCatMean;
    velNonstimTotal(i) = velNonstimCatMean;
end

velTotal = [velNonstimTotal velStimTotal];

parts = strsplit(pathname, '\');
filename = [parts{end-2} '_' parts{end-1} ...
    '_vel_stim_nonstim_mean.mat'];
save(filename, 'velTotal');
fprintf('velocity data saved. \n');