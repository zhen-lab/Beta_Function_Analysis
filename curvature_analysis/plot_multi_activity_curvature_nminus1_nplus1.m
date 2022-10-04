%% Clear everything

clear; close all

%% Load ratio and curvature data

rfilename = uigetfile('*.mat', 'Select activity file');
if rfilename==0
    fprintf('user cancelled selection.\n');
else
    load(rfilename);
    if ~exist('ratio','var')
        fprintf('no ratio data exists.\n');
    else
        pref = regexp(rfilename, '.*(?=_Ratio)', 'match');
        if isempty(pref)
            pref = regexp(rfilename, '.*(?=[.]mat)', 'match');
        end
        cfilename = uigetfile('*.mat', 'Select curvature file', [pref{1,1} '*.mat']);
        load(cfilename);
        if ~exist('curv_cor','var')&&~exist('curvdata','var')
            fprintf('no curvature data exists.\n');
        else
            endframe = size(ratio,1);
            endneuron = size(ratio,2);
            prompt = {'Enter start frame:','Enter end frame:'...
                'Enter start neuron:', 'Enter end neuron:'};
            dlgtitle = 'Range';
            dims = [1 35];
            definput = {'1',num2str(endframe),'1',num2str(endneuron)};
            answer = inputdlg(prompt,dlgtitle,dims,definput);
            if isempty(answer)
                fprintf('user cancelled input.\n');
            else
                rangeframe = str2double(answer{1,1}):str2double(answer{2,1});
                rangeneuron = str2double(answer{3,1}):str2double(answer{4,1});
            end
        end
    end
end

type = 1;

%%

prompt = {'Enter start of neuron:','Enter end of neuron:',...
    'Enter start of frame:','Enter end of frame:'};
dlgtitle = 'Specs for files';
dims = [1 35];
definput = {'1','2','1','955'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
rangeneuron = str2double(answer{1,1}):str2double(answer{2,1});
rangeframe = str2double(answer{3,1}):str2double(answer{4,1});
ratio = zeros(length(rangeframe),length(rangeneuron));
curv_cor = zeros(length(rangeframe),length(rangeneuron));

for idx = rangeneuron
    afile = uigetfile('*.mat',['Select ' num2str(idx) ' activity file']);
    if afile==0
        fprintf('user cancelled selection.\n');
        break
    else
        load(afile,'signal','signal_mirror');
        ratio(:,idx) = signal{1,1}(rangeframe)./signal_mirror{1,1}(rangeframe);
    end
end
fprintf('activity files loaded.\n');

for idx = rangeneuron(1):rangeneuron(end)+2
    cfile = uigetfile('*.mat',['Select ' num2str(idx) ' curvature file']);
    if cfile==0
        fprintf('user cancelled selection.\n');
        break
    else
        load(cfile);
        curv_cor(:,idx) = curvdatafiltered(rangeframe);
    end
end
fprintf('curvature files loaded.\n');

type = 2;

%% 
% DA #20x_10
rev = [492 777];

% DB #12_34-1800
% rev = [576 605; 907 955];

% neuronnum = 6;
% f_curv = uigetfile('*Curvature.mat', 'Select curvature data');
% f_activity = uigetfile('*Ratio.mat', 'Select activity data');
% load(f_curv);
% load(f_activity);

neuronnum = length(rangeneuron);
curv_norm = smoothdata(normalize_signal(curv_cor(rangeframe,:)), 'rloess');
ratio_norm = smoothdata(normalize_signal(ratio(rangeframe,:)), 'rloess');

%% Draw plot

numcomp = 3; stepcomp = 10;
% cmp = [0 0.5 1; 0.3 0.7 1; 0.9 1 1];
cmp = gray((numcomp+1)*stepcomp);
lm = [0 1];
loop = rangeneuron(1)+1:rangeneuron(end)+1;
subrow = neuronnum;

f = figure;
for idx = loop
    subplot(subrow, 1, idx-1);
        hold on;
    for blk = 1:size(rev,1)
        fill([rev(blk,1) rev(blk,2) rev(blk,2) rev(blk,1)],...
                [lm(1) lm(1) lm(2) lm(2)], 0.8*[1 1 1], 'edgecolor', 'none');
    end
    for curnnum = -1:1
        plot(curv_norm(:, idx+curnnum), ...
            'color', cmp((curnnum+2.5)*stepcomp, :), 'linewidth', 2); 
    end
    plot(ratio_norm(:, idx-1), 'color', [1 0.1 0], 'linewidth', 4);  
    set(gca, 'visible', 'off');
end

%% Save plot

if ~exist('pref','var')||~exist('f','var')
    fprintf('please load data and draw plot first.\n');
else
    filenamefig = [pref{1} '_gray.fig'];
    filenameemf = [pref{1} '_gray.emf'];
    savefig(f, filenamefig);
    saveas(f, filenameemf, 'meta');
end
