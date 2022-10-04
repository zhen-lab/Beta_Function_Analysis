%% Clear everything

clear; close all;

%% Set up inputs

filename = uigetfile('*.mat');
load(filename);
input_d = dorsal_ratio;
input_v = ventral_ratio;

expORctrl = questdlg('Experimental or Control group?', 'Exp/Ctrl', ...
    'Experimental', 'Control', 'Experimental');

prompt = {'Minimum of y axis','Maximum of y axis'};
dlgtitle = 'Y axis range';
dims = [1 35];
definput = {'0','3'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
ymin = str2double(answer{1});
ymax = str2double(answer{2});

s_exp = [];
s_ctrl = [];

switch expORctrl
    case 'Experimental'
        % Plot green lines for all-optical data
        close all; 
        s_exp = im_DV_data(input_d, input_v, [ymin ymax], [ymin ymax], [0 1]);
    case 'Control'
        % Plot black lines for all-optical data
        close all; 
        s_ctrl = im_DV_data(input_d, input_v, [ymin ymax], [ymin ymax], [1 0]);
end

%% Save figures for all-optical data

parts = strsplit(pwd, '\');
% regexp(parts, 'Alpha\w*', 'match');
answer = questdlg('Subfolder exists?', 'Subfolder', 'Yes', 'No', 'Yes');
switch answer
    case 'Yes'
        subfd = 1;
    case 'No'
        subfd = 0;
end
data_path_fig = fullfile(parts{1,1:end-2-subfd}, 'Alpha_Data_Plot', parts{1,end-subfd});
warning('off'); mkdir(data_path_fig);

if ~isempty(s_exp)
    s_save = s_exp;
    pref = 'Experimental_';
elseif ~isempty(s_ctrl)
    s_save = s_ctrl;
    pref = 'Control_';
else
    fprintf('no figures have been generated yet.\n');
end
% s_save = ~isempty(s_exp)*s_exp + ~isempty(s_ctrl)*s_ctrl;
% pref = ~isempty(s_exp)*'Experimental_' + ~isempty(s_ctrl)*'Control_';
savefig(s_save, fullfile(data_path_fig, [pref parts{end} '.fig']));
saveas(s_save(1), fullfile(data_path_fig, [pref parts{end} '_dorsal_overlapped.emf']), 'meta');
saveas(s_save(2), fullfile(data_path_fig, [pref parts{end} '_ventral_overlapped.emf']), 'meta');
saveas(s_save(3), fullfile(data_path_fig, [pref parts{end} '_dorsal_separated.emf']), 'meta');
saveas(s_save(4), fullfile(data_path_fig, [pref parts{end} '_ventral_separated.emf']), 'meta');
fprintf(['figures saved for ' parts{end-1} '_' parts{end} ' ' pref 'group. \n']);
% end
% else
%     if ~isempty(s_ctrl)
%         pref = 'Control_';
%         savefig(s_ctrl, fullfile(data_path_fig, [pref parts{end} '.fig']));
%         saveas(s_ctrl(1), fullfile(data_path_fig, [pref parts{end} '_dorsal_overlapped.emf']), 'meta');
%         saveas(s_ctrl(2), fullfile(data_path_fig, [pref parts{end} '_ventral_overlapped.emf']), 'meta');
%         saveas(s_ctrl(3), fullfile(data_path_fig, [pref parts{end} '_dorsal_separated.emf']), 'meta');
%         saveas(s_ctrl(4), fullfile(data_path_fig, [pref parts{end} '_ventral_separated.emf']), 'meta');
%         fprintf('figures saved for CONTROL group. \n');
%     end
% end

% if ~isempty(s_exp)
%     pref = 'Experimental_';
%     savefig(s_exp, fullfile(data_path_fig, [pref parts{end} '.fig']));
%     saveas(s_exp(1), fullfile(data_path_fig, [pref parts{end} '_dorsal_overlapped.tif']), 'tiffn');
%     saveas(s_exp(2), fullfile(data_path_fig, [pref parts{end} '_ventral_overlapped.tif']), 'tiffn');
%     saveas(s_exp(3), fullfile(data_path_fig, [pref parts{end} '_dorsal_separated.tif']), 'tiffn');
%     saveas(s_exp(4), fullfile(data_path_fig, [pref parts{end} '_ventral_separated.tif']), 'tiffn');
%     fprintf('figures saved for EXPERIMENTAL group. \n');
% else
%     if ~isempty(s_ctrl)
%         pref = 'Control_';
%         savefig(s_ctrl, fullfile(data_path_fig, [pref parts{end} '.fig']));
%         saveas(s_ctrl(1), fullfile(data_path_fig, [pref parts{end} '_dorsal_overlapped.tif']), 'tiffn');
%         saveas(s_ctrl(2), fullfile(data_path_fig, [pref parts{end} '_ventral_overlapped.tif']), 'tiffn');
%         saveas(s_ctrl(3), fullfile(data_path_fig, [pref parts{end} '_dorsal_separated.tif']), 'tiffn');
%         saveas(s_ctrl(4), fullfile(data_path_fig, [pref parts{end} '_ventral_separated.tif']), 'tiffn');
%         fprintf('figures saved for CONTROL group. \n');
%     end
% end

%% Plot together for non-all-optical data

ec_dv = questdlg('Plot exp/ctrl or d/v plots?', 'Plot type', 'Exp/Ctrl', 'D/V', 'Exp/Ctrl');

switch ec_dv
    
    case 'Exp/Ctrl'
       
        close all;
        file_c = uigetfile('.mat', 'Select control group data');
        load(file_c);
        dorsal_c = dorsal; ventral_c = ventral;
        file_h = uigetfile('.mat', 'Select experimental group data');
        load(file_h);
        dorsal_h = dorsal; ventral_h = ventral;

        dorsal_c(dorsal_c==0) = NaN;
        ventral_c(ventral_c==0) = NaN;
        dorsal_h(dorsal_h==0) = NaN;
        ventral_h(ventral_h==0) = NaN;
        
        dorsal_c(dorsal_c==inf) = NaN;
        ventral_c(ventral_c==inf) = NaN;
        dorsal_h(dorsal_h==inf) = NaN;
        ventral_h(ventral_h==inf) = NaN;

        % Plot experimental and control groups together
        ymin = floor(min([dorsal_c ; ventral_c ; dorsal_h ; ventral_h], [], 'all'));
        ymax = ceil(max([dorsal_c ; ventral_c ; dorsal_h ; ventral_h], [], 'all'));
        color_ctrl = [0.8 0.8 0.8; 0 0 0];
        color_exp = [1 0.6 0.6; 1 0 0];
%         color_ctrl = [0.8 0.8 0.8; 0 0 0];
%         color_exp = [1 0.6 0.3; 1 0.3 0];
%         color_exp = [0.5 1 0.5; 0 0.9 0];

        s(1) = figure(1); hold on
        plot(dorsal_c', 'color', color_ctrl(1,:), 'linewidth', 2);
        plot(nanmean(dorsal_c), 'color', color_ctrl(2,:), 'linewidth', 4);
        plot(dorsal_h', 'color', color_exp(1,:), 'linewidth', 2);
        plot(nanmean(dorsal_h), 'color', color_exp(2,:), 'linewidth', 4);
        set(gca, 'visible', 'off');
        set(gcf, 'color', 'w');
%         set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 0.5 1]);
        ylim([ymin ymax]); xlim([1 length(dorsal_c)]);

        s(2) = figure(2); hold on
        plot(ventral_c', 'color', color_ctrl(1,:), 'linewidth', 2);
        plot(nanmean(ventral_c), 'color', color_ctrl(2,:), 'linewidth', 4);
        plot(ventral_h', 'color', color_exp(1,:), 'linewidth', 2);
        plot(nanmean(ventral_h), 'color', color_exp(2,:), 'linewidth', 4);
        ylim([ymin ymax]); xlim([1 length(ventral_c)]);
        set(gca, 'visible', 'off');
        set(gcf, 'color', 'w');
%         set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 0.5 1]);
        fileused = file_c;

    case 'D/V'
        
        close all;
        filedv = uigetfile('.mat', 'Select data for dorsal and ventral muscle');
        load(filedv);
        
        dorsal(dorsal==0) = NaN;
        ventral(ventral==0) = NaN;
        
        ymin = floor(1.2*min([dorsal; ventral], [], 'all'));
        ymax = ceil(1.2*max([dorsal; ventral], [], 'all'));
        color_d = [0.8 0.8 0.8; 0 0 0];
        color_v = [1 0.6 0.6; 1 0 0];
%         color_d = [1 0.6 1; 0.9 0 0.9];
%         color_v = [0.6 1 0.6; 0 0.9 0];
        
        % Plot dorsal and ventral together
        s = figure; hold on
        plot(dorsal', 'color', color_d(1,:), 'linewidth', 2);
        plot(nanmean(dorsal), 'color', color_d(2,:), 'linewidth', 4);
        plot(ventral', 'color', color_v(1,:), 'linewidth', 2);
        plot(nanmean(ventral), 'color', color_v(2,:), 'linewidth', 4);
        
        set(gca, 'visible', 'off')
        set(gcf, 'color', 'w') 
        
        fileused = filedv;
end

%% Save figures for non-all-optical data

answer = questdlg('Subfolder exists?', 'Subfolder', 'Yes', 'No', 'No');
switch answer
    case 'Yes'
        subfd = 1;
    case 'No'
        subfd = 0;
end
parts = strsplit(pwd, '\');
data_path_fig = fullfile(parts{1,1:end-2-subfd}, 'Alpha_Data_Plot', parts{1,end-subfd});
warning('off'); mkdir(data_path_fig); 

switch length(s)
    case 1
        data_path_name_fig = fullfile(data_path_fig, ['dorsal-ventral_' fileused(end-22:end)]);
        savefig(s, [data_path_name_fig(1:end-4) '_dorsal-ventral.fig']);
        saveas(s, [data_path_name_fig(1:end-4) '_dorsal-ventral.emf'], 'meta');
    case 2
        data_path_name_fig = fullfile(data_path_fig, ['dorsal-ventral_' fileused(end-22:end)]);
        savefig(s, [data_path_name_fig(1:end-4) '_exp-ctrl.fig']);
        saveas(s(1), [data_path_name_fig(1:end-4) '_exp-ctrl_dorsal.emf'], 'meta');
        saveas(s(2), [data_path_name_fig(1:end-4) '_exp-ctrl_ventral.emf'], 'meta');
end

fprintf('figures saved for non-all-optical data. \n');

%% Output mean value for each non-all-optical data point

avg_dorsal_c = nanmean(dorsal_c, 2);
avg_ventral_c = nanmean(ventral_c, 2);
avg_dorsal_h = nanmean(dorsal_h, 2);
avg_ventral_h = nanmean(ventral_h, 2);

save('Mean_Muscle_Interneurons_Ablated_unc-25_Ablated_Ctrl', 'avg_dorsal_c', 'avg_dorsal_h', 'avg_ventral_c', 'avg_ventral_h');