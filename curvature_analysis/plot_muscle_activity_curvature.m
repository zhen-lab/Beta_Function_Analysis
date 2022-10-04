%% Clear everything

clear; close all;

%% Determine curvature direction if necessary

filename = uigetfile('.mat');

if filename==0
    
    fprintf('No file is chosen. \n');

else
    
    load(filename);
    close all;
    curvdataBody = curvdatafiltered;

    figure; 
    imnumstart = 1; 
%     imnumend = 12;
    imnumend = size(curvdatafiltered,2);
    subplot(2,4,[2,3]); hold on;
    plot(dorsal_data{2*imnumstart-1,1}, dorsal_data{2*imnumstart,1}, ':r');
    plot(ventral_data{2*imnumstart-1,1}, ventral_data{2*imnumstart,1}, ':b');
    plot(centerline_data_spline(:,2*imnumstart-1), centerline_data_spline(:,2*imnumstart), ':k');
    plot(centerline_data_spline(1,2*imnumstart-1), centerline_data_spline(1,2*imnumstart), ':og');
    plot(centerline_data_spline(end,2*imnumstart-1), centerline_data_spline(end,2*imnumstart), ':oy');
    plot(dorsal_data{2*imnumend-1,1}, dorsal_data{2*imnumend,1}, 'r');
    plot(ventral_data{2*imnumend-1,1}, ventral_data{2*imnumend,1}, 'b');
    plot(centerline_data_spline(:,2*imnumend-1), centerline_data_spline(:,2*imnumend), 'k');
    plot(centerline_data_spline(1,2*imnumend-1), centerline_data_spline(1,2*imnumend), 'og');
    plot(centerline_data_spline(end,2*imnumend-1), centerline_data_spline(end,2*imnumend), 'oy');
    set(gca, 'ydir', 'reverse');
    axis equal;
    subplot(2,4,5); imagesc(curvdataBody); caxis([-10 10]); title('Curvature unchanged');
    subplot(2,4,6); imagesc(-curvdataBody); caxis([-10 10]); title('Curvature flipped');
    subplot(2,4,7); imagesc(dorsal_smd./dorsal_smd_r); title('Dorsal activity');
    subplot(2,4,8); imagesc(ventral_smd./ventral_smd_r); title('Ventral activity');
    set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
    
    answer = questdlg(['Flip curvature for ' filename ' ?'], 'Sign of curvature', ...
        'No', 'Yes', 'No');
    
    if isequal(answer,'Yes')
        curvdataBody = -curvdatafiltered;
        fprintf(['Curvature data is flipped for ' filename '\n']);
%         subplot(142); imagesc(curvdataBody); title('Curvature after flipping');
        data_path_name = fullfile(pwd, ['cur-corrected_' filename(1:end-4) '.mat']);
        save(data_path_name, 'curvdataBody');
        fprintf(['data saved for ' filename '\n']);
    elseif isequal(answer,'No')
        fprintf(['Curvature data remains unchanged for ' filename '\n']);
        data_path_name = fullfile(pwd, ['cur-corrected_' filename(1:end-4) '.mat']);
        save(data_path_name, 'curvdataBody');
        fprintf(['data saved for ' filename '\n']);
    end
    
end

close all;

%% Select files if curvature already curated
% Only excute this section if necessary!!!!!

filenamea = uigetfile('*.mat', 'Select activity file');
load(filenamea);
filenamec = uigetfile('*.mat', 'Select curvature file');
load(filenamec);
filename = filenamea;

%% Generate figures

prompt = {'Minimum of activity','Maximum of activity', ...
    'Minimum of curvature','Maximum of curvature'};
inputdlgtitle = 'Range and channel order';
dims = [1 35];
definput = {'0','3','-10', '10'};
answer = inputdlg(prompt,inputdlgtitle,dims,definput);
body = 35;
hbline = [1 size(dorsal_smd,2); body body];
curv = curvdataBody; % define this everytime, otherwise it might be transposed from previous trials
cloc = 'EastOutside';

quest = 'GFP and RFP channel flipped?';
questdlgtitle = 'Channel order'; 
flp = questdlg(quest, questdlgtitle, 'Yes', 'No', 'No');
switch flp
    case 'Yes'
        ratio_d = dorsal_smd_r./dorsal_smd;
        ratio_v = ventral_smd_r./ventral_smd;
    case 'No'
        ratio_d = dorsal_smd./dorsal_smd_r;
        ratio_v = ventral_smd./ventral_smd_r;
end

quest = 'Transpose the plots?';
questdlgtitle = 'Transposition'; 
trp = questdlg(quest, questdlgtitle, 'Yes', 'No', 'No');
if isequal(trp, 'Yes')
    ratio_d = ratio_d';
    ratio_v = ratio_v';
    curv = curv';
    hbline = flipud(hbline);
    cloc = 'NorthOutside';
end

minia = str2double(answer{1});
maxia = str2double(answer{2});
minic = str2double(answer{3});
maxic = str2double(answer{4});

close all;
cmp = colormap(flipud(gray));
% deci = 1;

s(1) = figure(1);
imagesc(ratio_d); 
caxis([minia maxia]); 
c = colorbar;
set(c, 'ticklabels', [], 'location', cloc);
% tks = get(c, 'limits');
% tkslin = linspace(tks(1),tks(2),3);
% labels = arrayfun(@(x) sprintf(['%.' num2str(deci) 'f'],x), tkslin, 'un', 0);
% set(c, 'ticks', tkslin, 'ticklabels', labels);
colormap(magma);
% colormap(cmp);
set(gca, 'ticklength', [0 0], 'xticklabel', [], 'yticklabel', []);
hold on;
plot(hbline(1,:), hbline(2,:), ':w', 'linewidth', 2);

s(2) = figure(2);
imagesc(ratio_v); 
caxis([minia maxia]);
c = colorbar;
set(c, 'ticklabels', [], 'location', cloc);
% tks = get(c, 'limits');
% tkslin = linspace(tks(1),tks(2),3);
% labels = arrayfun(@(x) sprintf(['%.' num2str(deci) 'f'],x), tkslin, 'un', 0);
% set(c, 'ticks', tkslin, 'ticklabels', labels);
colormap(magma);
% colormap(cmp);
set(gca, 'ticklength', [0 0], 'xticklabel', [], 'yticklabel', []);
hold on;
plot(hbline(1,:), hbline(2,:), ':w', 'linewidth', 2);

s(3) = figure(3);
imagesc(curv); 
caxis([minic maxic]);
c = colorbar;
set(c, 'ticklabels', [], 'location', cloc);
% tks = get(c, 'limits');
% tkslin = linspace(tks(1),tks(2),3);
% labels = arrayfun(@(x) sprintf(['%.' num2str(deci)-1 'f'],x), tkslin, 'un', 0);
% set(c, 'ticks', tkslin, 'ticklabels', labels);
% colormap(viridis);
colormap(cmp);
set(gca, 'ticklength', [0 0], 'xticklabel', [], 'yticklabel', []);
hold on;
plot(hbline(1,:), hbline(2,:), ':w', 'linewidth', 2);

%% Add contour to the figures if necessary

hold on;
thr = 11; 
M = contour(curvdataBody, [thr thr]);
plot(M(1,2:end), M(2,2:end), ':k', 'linewidth', 3);

%% Save figures

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
data_path_name_fig = fullfile(data_path_fig, [filename(1:end-4) '.mat']);
savefig(s, [data_path_name_fig(1:end-4) '_act-cur']);
saveas(s(1), [data_path_name_fig(1:end-4) '_act_dorsal'], 'meta');
saveas(s(2), [data_path_name_fig(1:end-4) '_act_ventral'], 'meta');
saveas(s(3), [data_path_name_fig(1:end-4) '_cur'], 'meta');
fprintf('figures saved. \n');
