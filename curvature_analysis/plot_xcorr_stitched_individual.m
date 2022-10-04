%% Clear everything 

clear; close all

%% Load data

filename = uigetfile('*.mat', 'Select xcorr data', '*_corr.mat');
if filename==0
    fprintf('user cancelled selection.\n');    
else
    pref = regexp(filename, '.*(?=[.]mat)', 'match');
    load(filename);
    if ~exist('lags_collect_for','var')
        fprintf('no lags and correlation data exist.\n');
    else
        prompt = {'Enter the first neuron to display:',...
            'Total number of neurons to display',...
            'Enter the range on x axis:'};
        dlgtitle = 'Specs';
        dims = [1 35];
        definput = {'9','2','100'};
        answer = inputdlg(prompt,dlgtitle,dims,definput);
        neuronnum = str2double(answer{1,1});
        neurontotal = str2double(answer{2,1});
        xrange = str2double(answer{3,1});
    end
end

%% Draw plot

close all;
f = figure;
for idx = 1:neurontotal
    thisneuron = neuronnum+idx-1;
    % Plot xcorr for forward phase
    subplot(neurontotal,2,idx*2-1);
    plot(lags_collect_for{1,thisneuron}, ...
        r_collect_for{1,thisneuron}, 'color', 'k', 'linewidth', 4);
    hold on;
    plot([0 0], [0 1], 'color', 'k');
    plot(phaseshift_collect_for(thisneuron), max(r_collect_for{1,thisneuron}), 'o', ...
                    'markersize', 12, 'markerfacecolor', [1 0 0], 'markeredgecolor', [1 0 0]);
    xlim([-xrange xrange]); ytickformat('%.1f');
    h = gca; h.XAxis.Visible = 'off'; box off
    % Plot xcorr for backward phase
    subplot(neurontotal,2,idx*2);
    plot(lags_collect_back{1,thisneuron}, r_collect_back{1,thisneuron}, 'color', 'k', 'linewidth', 4);
    hold on;
    plot([0 0], [0 1], 'color', 'k');
    plot(phaseshift_collect_back(thisneuron), max(r_collect_back{1,thisneuron}), 'o', ...
                    'markersize', 12, 'markerfacecolor', [1 0 0], 'markeredgecolor', [1 0 0]);
    xlim([-xrange xrange]); ytickformat('%.1f');
    h = gca; h.XAxis.Visible = 'off'; box off    
end

%% Save plot

if ~exist('pref','var')||~exist('f','var')
    fprintf('please load data and draw plot first.\n');
else
    filenamefig = [pref{1} '_xcorr.fig'];
    filenameemf = [pref{1} '_xcorr.emf'];
    savefig(f, filenamefig);
    saveas(f, filenameemf, 'meta');
end
