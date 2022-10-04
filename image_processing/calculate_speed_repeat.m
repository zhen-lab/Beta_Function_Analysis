%% Repeat this step for each recording until all analyses are done

clear; close all
calculate_speed;

%% Collect forward and backward data from all data

% Navigate to the folder that contains all velocity data
numneurons = 2;
[forward, backward, ~, ~] = collect_FB_data(pwd, numneurons);

%% Draw pooled signals in forward and backward movements respectively

answer = questdlg('Are neurons of the same type?', 'Neuron types', 'Yes', 'No', 'Yes');
switch answer
    case 'Yes'
        forward_plot = cat(1, forward(:));
        backward_plot = cat(1, backward(:));
    case 'No'
        forward_plot = forward;
        backward_plot = backward;
end
plot_forward_backward(forward_plot, backward_plot);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To draw plots separately from analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[f_v, ~] = uigetfile('*velocity.mat');
load(f_v);
min_consecutive = 2;
forward = double(vel_ap_sign_smd_mean > 0); 
forward_buffered = [0; double(forward); 0]; % To enclose pulses on both ends
backward = double(vel_ap_sign_smd_mean < 0); 
backward_buffered = [0; double(backward); 0]; % To enclose pulses on both ends
[W_forward,INITCROSS_forward,FINALCROSS_forward] = pulsewidth(double(forward_buffered));
[W_backward,INITCROSS_backward,FINALCROSS_backward] = pulsewidth(double(backward_buffered));

% Remove blocks with width no greater than the threshold for forward
forward_pruned = forward;
forward_pruned((ceil(INITCROSS_forward(W_forward<=min_consecutive))-1)...
:(floor(FINALCROSS_forward(W_forward<=min_consecutive))-1)) = 0;
INITCROSS_forward(W_forward<=min_consecutive) = [];
FINALCROSS_forward(W_forward<=min_consecutive) = [];
W_forward(W_forward<=min_consecutive) = [];
if isequal(size(W_forward), size(INITCROSS_forward), size(FINALCROSS_forward))
fprintf(['Forward Events #: ' num2str(size(W_forward, 1)) '\n']);
else
fprintf(['Forward Unequal: ' num2str(size(W_forward)) ' ' num2str(size(INITCROSS_forward)) ' ' num2str(size(FINALCROSS_forward)) '\n']);
end
% Remove blocks with width no greater than the threshold for backward
backward_pruned = backward;
backward_pruned((ceil(INITCROSS_backward(W_backward<=min_consecutive))-1)...
:(floor(FINALCROSS_backward(W_backward<=min_consecutive))-1)) = 0;
INITCROSS_backward(W_backward<=min_consecutive) = [];
FINALCROSS_backward(W_backward<=min_consecutive) = [];
W_backward(W_backward<=min_consecutive) = [];
if isequal(size(W_backward), size(INITCROSS_backward), size(FINALCROSS_backward))
fprintf(['Backward Events #: ' num2str(size(W_backward, 1)) '\n']);
else
fprintf(['Forward Unequal: ' num2str(size(W_backward)) ' ' num2str(size(INITCROSS_backward)) ' ' num2str(size(FINALCROSS_backward)) '\n']);
end
blocks_forward = ~isempty(W_forward)*size(W_forward,1);
blocks_backward = ~isempty(W_backward)*size(W_backward,1);

%% Load activities from anterior and posterior neurons

[f_a, ~] = uigetfile('.mat');
load(f_a); signal_anterior = signal{1,1}./signal_mirror{1,1};
[f_p, ~] = uigetfile('.mat');
load(f_p); signal_posterior = signal{1,1}./signal_mirror{1,1};

%% Draw heatmap and perievent plots

close all;
ratio_all = [signal_anterior signal_posterior];
reversals = [INITCROSS_backward FINALCROSS_backward];
fn = regexp(f_v, '\w*\d', 'match'); filename = [fn{1,1} '_heatmap_perievent'];
plot_heatmap_perievent(ratio_all, filename, reversals);

%% Draw plot for activities only

colorchoice = questdlg('Cholinergic or GABAergic neurons?', 'Color', ...
    'Cholinergic', 'GABAergic', 'Cholinergic');

switch colorchoice
    case 'Cholinergic'
        cmp = lbmap(10, 'redblue');
        clr_ant = 3; clr_pos = 4;
    case 'GABAergic'
        cmp = lbmap(10, 'blue'); 
        clr_ant = 4; clr_pos = 9;
end

close all;
ymax = 1.2*max([signal_anterior signal_posterior],...
    [], 'all', 'omitnan');

figure;
if ~(blocks_backward==0)
    
    for blocks = 1:length(W_backward)
        hold on;
        st = ceil(INITCROSS_backward(blocks))-1; en = floor(FINALCROSS_backward(blocks))-1;
        fill([st en en st], 1.5*[0 0 ymax ymax], 0.8*[1 1 1], 'edgecolor', 'none');
    end
    plot(signal_anterior, 'color', cmp(clr_ant,:), 'linewidth', 2);
    plot(signal_posterior, 'color', cmp(clr_pos,:), 'linewidth', 2);
    axis tight; set(gca, 'xticklabel', [], 'visible', 'off'); hold off;

end

%% Draw plot for activities and curvature

curvfile = uigetfile('*_curvature_curated.mat', ...
    'Select curvature data for anterior neuron');
load(curvfile, 'curvdatafiltered');
curv_anterior = curvdatafiltered;
curvfile = uigetfile('*_curvature_curated.mat', ...
    'Select curvature data for posterior neuron');
load(curvfile, 'curvdatafiltered');
curv_posterior = curvdatafiltered;

signal_anterior_smd = smooth(signal_anterior,...
    0.05, 'rloess');
signal_posterior_smd = smooth(signal_posterior,...
    0.05, 'rloess');

clr = gray(40); % refer to plot_multi_activity_curvature_nminus1_nplus1.m
% smoothsigymax = 1.2*max([signal_anterior_smd signal_posterior_smd],...
%     [], 'all', 'omitnan');
% curvymax = 1; curvymin = -2;
%%

signal_anterior_smd_norm = normalize_signal(signal_anterior_smd);
signal_posterior_smd_norm = normalize_signal(signal_posterior_smd);
curv_anterior_norm = normalize_signal(curv_anterior);
curv_posterior_norm = normalize_signal(curv_posterior);

close all;
if ~(blocks_backward==0)
        
    % Plot anterior data
    subplot(2,1,1); set(gca, 'visible', 'off');   
    for blocks = 1:length(W_backward)
        hold on;
        st = ceil(INITCROSS_backward(blocks))-1; en = floor(FINALCROSS_backward(blocks))-1;
        fill([st en en st], [0 0 1 1], 0.8*[1 1 1], 'edgecolor', 'none');
    end
%     yyaxis left; 
     axis tight; hold on;
%     yyaxis right; 
    plot(curv_anterior_norm, ...
        'color', clr(15,:), 'linewidth', 4); % ylim([curvymin curvymax]);
    % choice of curvature color: refer to plot_multi_activity_curvature_nminus1_nplus1.m
    plot(signal_anterior_smd_norm, ...
        'color', [1 0.1 0], 'linewidth', 4);
    
    % Plot posterior data
    subplot(2,1,2); set(gca, 'visible', 'off');
    for blocks = 1:length(W_backward)
        hold on;
        st = ceil(INITCROSS_backward(blocks))-1; en = floor(FINALCROSS_backward(blocks))-1;
        fill([st en en st], [0 0 1 1], 0.8*[1 1 1], 'edgecolor', 'none');
    end
%     yyaxis left; 
    axis tight; hold on;
%     yyaxis right; 
    plot(curv_posterior_norm, ...
        'color', clr(15,:), 'linewidth', 4); % ylim([curvymin curvymax]);
    % choice of curvature color: refer to plot_multi_activity_curvature_nminus1_nplus1.m
    plot(signal_posterior_smd_norm, ...
        'color', [1 0.1 0], 'linewidth', 4);
    
end

%% Save figure

figurename = [f_v(1:end-13) '_activity_curvature'];
savefig(figurename);
saveas(gca, figurename, 'meta');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To draw plots separately from analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%