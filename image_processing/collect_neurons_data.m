function [GCaMP, RFP, ratio, ratio_norm] = collect_neurons_data(folder, frames)

f = dir(fullfile(folder, '*.mat'));
numsamples = numel(f);
% GCaMP = zeros(frames, numsamples);
% RFP = zeros(frames, numsamples);
% ratio = zeros(frames, numsamples);

% Determine data types
name = fullfile(folder, f(1).name);
load(name);
% Case of one neuron per data
if exist('signal','var')
    GCaMP = zeros(frames, numsamples);
    RFP = zeros(frames, numsamples);
    ratio = zeros(frames, numsamples);
    for idx = 1:numsamples
        name = fullfile(folder, f(idx).name);
        load(name, 'signal', 'signal_mirror');
        frame_length = min(size(signal{1,1},1), frames); % in case there are recordings with shorter time window
        GCaMP(1:frame_length, idx) = signal{1,1}(1:frame_length);
        RFP(1:frame_length, idx) = signal_mirror{1,1}(1:frame_length);
        ratio(1:frame_length, idx) = signal{1,1}(1:frame_length)./signal_mirror{1,1}(1:frame_length);
    %     elseif exist('gfp','var')
    %         frame_length = min(size(gfp,1), frames); % in case there are recordings with shorter time window
    %         for numneuron = 1:size(gfp,2)
    %             idxupdated = idx+numneuron-1;
    %             GCaMP(1:frame_length, idxupdated) = gfp(1:frame_length,numneuron);
    %             RFP(1:frame_length, idxupdated) = rfp(1:frame_length,numneuron);
    %             ratio(1:frame_length, idxupdated) = gfp(1:frame_length,numneuron)./rfp(1:frame_length,numneuron);
    %         end
    end
% Case of multiple neurons per data
elseif exist('sig_gfp','var')
%     colnum = 1;
    GCaMP = [];
    RFP = [];
    ratio = [];
    for idx = 1:numsamples
        name = fullfile(folder, f(idx).name);
        load(name, 'sig_gfp', 'sig_rfp');
%         neuronnum = size(sig_gfp,2);
        frame_length = min(size(sig_gfp,1), frames); % in case there are recordings with shorter time window
        GCaMP = cat(2, GCaMP, sig_gfp(1:frame_length,:));
        RFP = cat(2, RFP, sig_rfp(1:frame_length,:));
        ratio = cat(2, ratio, sig_gfp(1:frame_length,:)./sig_rfp(1:frame_length,:));
%         colnum = colnum+neuronnum;
    end
end

ratio_start = repmat(ratio(1, :), frames, 1);
ratio_norm = (ratio./ratio_start)';

GCaMP = GCaMP'; RFP = RFP'; ratio = ratio';

close all;
plot(ratio_norm');

end