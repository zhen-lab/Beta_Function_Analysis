%% Clear everything

clear; close all;

%% Load data

filename = uigetfile('*.mat', 'Select activity file');
if filename==0
    fprintf('user cancelled selection.\n');
else
    load(filename);
    if ~exist('ratio','var')
        fprintf('no ratio data exists.\n');
    elseif ~exist('signal','var') % for the cases of DA and DB data
        signalall = ratio;
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
    else % for the case of DD data
        endframe = size(signal{1,1},1);
        prompt = {'Enter start frame:','Enter end frame:'...
            'Enter number of neurons for analysis:'};
        dlgtitle = 'Range';
        dims = [1 35];
        definput = {'1',num2str(endframe),'2'};
        answer = inputdlg(prompt,dlgtitle,dims,definput);
        if isempty(answer)
            fprintf('user cancelled input.\n');
        else
            rangeframe = str2double(answer{1,1}):str2double(answer{2,1});
            rangeneuron = 1:str2double(answer{3,1});
        end
        signalall = zeros(length(rangeframe),length(rangeneuron));
        for neuronnum = rangeneuron
            rfilename = uigetfile('*.mat', 'Select signal file');
            if rfilename==0
                fprintf('user cancelled selection.\n');
            else
                load(rfilename);
                signalall(:,neuronnum) = signal{1,1}./signal_mirror{1,1};
            end
        end
    end
end

%% Draw plot

f = figure;
neuronactivity_smd = smoothdata(signalall(rangeframe,rangeneuron), 'rloess');
% neuronactivity_norm = normalize_signal(neuronactivity_smd);
imagesc(neuronactivity_smd'); 
set(gca, 'xticklabel', [], 'yticklabel', [], 'ticklength', [0 0]);
colormap(flipud(gray));

clim = caxis;
prompt = {'Enter color axis min:','Enter color axis max:'};
dlgtitle = 'Caxis';
dims = [1 35];
definput = {num2str(clim(1)),num2str(clim(2))};
colorrange = inputdlg(prompt,dlgtitle,dims,definput);
if isempty(colorrange)
    fprintf('user cancelled input, caxis remains default.\n');
elseif str2double(colorrange{1,1})<str2double(colorrange{2,1})
    caxis([str2double(colorrange{1,1}) str2double(colorrange{2,1})])
    fprintf(['caxis updated:' colorrange{1,1} ' ' colorrange{2,1} '\n']);
else
    fprintf('please input a smaller number for min and larger number for max.\n');
end

c = colorbar;
deci = 1;
numtick = 3;
t = get(c, 'limits');
T = linspace(t(1), t(2), numtick);
TL = arrayfun(@(x) sprintf(['%.' num2str(deci) 'f'], x), T, 'un', 0);
set(c, 'Ticks', T, 'TickLabels', TL, 'FontSize', 12);

%% Save plot

pref = regexp(filename, '.*(?=[.]mat)', 'match');
filenamefig = [pref{1} '_gray.fig'];
filenameemf = [pref{1} '_gray.emf'];
savefig(f, filenamefig);
saveas(f, filenameemf, 'meta');
