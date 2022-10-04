%% Clear everything

clear; close all;

%% Load data and draw figure

allfiles = {'Exp','Ctrl'};
filenum = size(allfiles,2);
ratio_norm_collect = cell(filenum,1);

for idx = 1:filenum

    filename = uigetfile('*.mat',['Select ' allfiles{idx} ' file']);
    load(filename);
    ratio_norm_collect{idx} = ratio_norm;
    
end

prompt = {'Minimum of y axis','Maximum of y axis'};
dlgtitle = 'Y axis range';
dims = [1 35];
definput = {'0','3'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
ymin = str2double(answer{1});
ymax = str2double(answer{2});

im_neurons_data_overlaid(ratio_norm_collect, [ymin ymax], 0.7);

%% Load data in case there is only experimental group

allfiles = 'Exp';
filename = uigetfile('*.mat',['Select ' allfiles ' file']);
load(filename);
ratio_norm_collect = cell(1,1);
ratio_norm_collect{:} = ratio_norm;

prompt = {'Minimum of y axis','Maximum of y axis'};
dlgtitle = 'Y axis range';
dims = [1 35];
definput = {'0','5'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
ymin = str2double(answer{1});
ymax = str2double(answer{2});

im_neurons_data_overlaid(ratio_norm_collect, [ymin ymax], 0.7);

%% Save figures

parts = strsplit(pwd, '\');

if ~isempty(s_exp)
    savefig(s_exp, ['Experimental_' parts{end} '.fig']);
    for i = 1:size(s_exp, 2)
        saveas(s_exp(i), ['Experimental_' num2str(i) '_' parts{end} '.emf'], 'meta');
    end
    fprintf('figures saved for EXPERIMENTAL group. \n');
else
    if ~isempty(s_ctrl)
        savefig(s_ctrl, ['Control_' parts{end} '.fig']);
        for i = 1:size(s_ctrl, 2)
            saveas(s_ctrl(i), ['Control_' num2str(i) '_' parts{end} '.emf'], 'meta');
        end
        fprintf('figures saved for CONTROL group. \n');
    end
end