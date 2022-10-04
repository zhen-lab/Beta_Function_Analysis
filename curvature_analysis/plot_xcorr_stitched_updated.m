%% Clear workspace

clear;

%% Clear figures and creat new figure

close all;
figure;

%%

% f = dir(fullfile(pwd, '*.mat'));
fnum = dir(fullfile(pwd, '*curated.mat'));
numsamples = numel(fnum);
count = 0;
r_collect_for = {};
lags_collect_for = {};
r_collect_back = {};
lags_collect_back = {};
phaseshift_collect_for = [];
phaseshift_collect_back = [];

thresh = 0.2;

for idx = 1:numsamples
    
    fnamecurvature = fnum(idx).name;
    fname = regexp(fnamecurvature, '\S*(?=_curvature_curated)', 'match');
    
    if ~isempty(fname)
        
        count = count + 1;
        fnameshort = regexp(fnamecurvature, 'temp\d', 'match');
%         fnameshort = regexp(fnamecurvature, '.*(?=_\D\d?_curvature)', 'match');

        fnamesignal = [fname{1,1} '.mat'];
        fnamevelocity = [fnameshort{1,1} '_velocity.mat'];
        load(fnamesignal, 'signal', 'signal_mirror');
        load(fnamecurvature, 'curvdatafiltered');
        load(fnamevelocity, 'vel_ap_sign_smd_mean');
        
        curation = questdlg(['Was ' fname{1,1} ' pruned for frames?'], ...
            'Curation', 'Yes', 'No', 'No');
        
        switch curation
            case 'No'
                deleted = [];
            case 'Yes'
                file_delfrm = uigetfile('*_deleted_frames.mat', 'Select deleted frames');
                if file_delfrm == 0 
                    fprintf('user cancelled selection. \n');
                    break
                else
                    load(file_delfrm, 'deleted');
                end
        end
        curated_frames = setdiff(1:size(signal{1,1},1), deleted);

        % Use curated frames
        signal_gfp_curated = signal{1,1}(curated_frames);
        signal_rfp_curated = signal_mirror{1,1}(curated_frames);
        curvdata_curated = curvdatafiltered(curated_frames, :);
        velocity_curated = vel_ap_sign_smd_mean(curated_frames(1:end-1), :);
        
%         signalratio = signal{1,1}./signal_mirror{1,1};
%         curvdatafiltered = curvdatafiltered - min(curvdatafiltered);
        signalratio = normalize_signal(smoothdata(signal_gfp_curated./signal_rfp_curated, 'rloess'));
        curvdatafiltered = normalize_signal(curvdata_curated);

%         subplot(131); plot(signalratio); 
%         hold on; plot(curvdatafiltered);

        forward = double(velocity_curated > thresh); 
        forward_buffered = [0; double(forward); 0]; % To enclose pulses on both ends
        backward = double(velocity_curated < -thresh); 
        backward_buffered = [0; double(backward); 0]; % To enclose pulses on both ends
        [~,INITCROSS_forward,FINALCROSS_forward] = pulsewidth(double(forward_buffered));
        [~,INITCROSS_backward,FINALCROSS_backward] = pulsewidth(double(backward_buffered));
        forwardseg = [INITCROSS_forward FINALCROSS_forward];
        backwardseg = [INITCROSS_backward FINALCROSS_backward];
        forwardtotal = []; backwardtotal = [];
        
        if ~isempty(forwardseg)
            for i = 1:size(forwardseg,1)
                forwardtotal = cat(2, forwardtotal, ceil(forwardseg(i,1)):floor(forwardseg(i,2)));
            end
            signalratioforward = signalratio(forwardtotal);
            curvforward = curvdatafiltered(forwardtotal);
            [rfor, lagsfor] = xcorr(signalratioforward, curvforward, 'normalized');
%             rnormfor = rfor/sqrt(sum(abs(signalratioforward).^2)*sum(abs(curvforward).^2));
            rnormfor = rfor;
            phaseshiftfor = lagsfor(rfor==max(rfor));
            r_collect_for{count} = rnormfor;
            lags_collect_for{count} = lagsfor;
            phaseshift_collect_for(count) = phaseshiftfor;
            subplot(211); hold on; 
            plot(lagsfor, rfor, 'color', 0.8*[1 1 1], 'linewidth', 2);
            plot(phaseshiftfor, max(rfor), 'o', ...
                'markersize', 8, 'markerfacecolor', [1 0 0], 'markeredgecolor', [1 0 0]);
        else 
            r_collect_for{count} = [];
            lags_collect_for{count} = [];
%             phaseshift_collect_for(count) = [];
        end
        
        if ~isempty(backwardseg)
            for i = 1:size(backwardseg,1)
                backwardtotal = cat(2, backwardtotal, ceil(backwardseg(i,1)):floor(backwardseg(i,2)));
            end
            signalratiobackward = signalratio(backwardtotal);
            curvbackward = curvdatafiltered(backwardtotal);
            [rback, lagsback] = xcorr(signalratiobackward, curvbackward, 'normalized');
%             rnormback = rback/sqrt(sum(abs(signalratiobackward).^2)*sum(abs(curvbackward).^2));
            rnormback = rback;
            phaseshiftback = lagsback(rback==max(rback));
            r_collect_back{count} = rnormback;
            lags_collect_back{count} = lagsback;
            phaseshift_collect_back(count) = phaseshiftback;
            subplot(212); hold on; 
            plot(lagsback, rback, 'color', 0.8*[1 1 1], 'linewidth', 2);
            plot(phaseshiftback, max(rback), 'o', ...
                'markersize', 8, 'markerfacecolor', [1 0 0], 'markeredgecolor', [1 0 0]);
        else
            r_collect_back{count} = [];
            lags_collect_back{count} = [];
%             phaseshift_collect_back(count) = [];
        end
   
    end
    
end

rx = 200;
subplot(211); plot([0 0], [0 1], 'k'); 
xlim([-rx rx]); ylim([0 1]);
subplot(212); plot([0 0], [0 1], 'k'); 
xlim([-rx rx]); ylim([0 1]);    

phaseshift_collect_for = phaseshift_collect_for';
phaseshift_collect_back = phaseshift_collect_back';

%%
path = fnum.folder;
pathpart = strsplit(path, '\');
save([pathpart{end} '_corr.mat'], ...
    'phaseshift_collect_for', 'phaseshift_collect_back', ...
    'r_collect_for', 'r_collect_back', ...
    'lags_collect_for', 'lags_collect_back', ...
    'forwardtotal', 'backwardtotal');
fprintf('xcorr data saved. \n');

% %%
% 
% figure; maxsize = 1800*2+1;
% paddingmatrix = zeros(maxsize, size(r_collect_for, 2));
% for i = 1:size(r_collect_for,2)
%     
%     subplot(121); 
%     hold on;
%     plot(lags_collect_for{i}, r_collect_for{i}, 'linewidth', 1, 'color', 0.8*[1 1 1]);
%     paddingmatrix(i,:) + 
%     
% end
% 
% 
% 
% for i = 1:size(r_collect_back,2)
%     subplot(122); hold on; 
%     plot(lags_collect_back{i}, r_collect_back{i}, 'linewidth', 1, 'color', 0.8*[1 1 1]);
%     
% end
%%

subplot(211); h = gca; h.XAxis.Visible = 'off';
subplot(212); h = gca; h.XAxis.Visible = 'off';
