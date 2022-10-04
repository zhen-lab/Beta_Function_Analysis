%% Clear everything

clear; close all;

%% Generate curvature if not done

generate_curvature_from_pts;

%% Select ratio and curvature data files for DAs

ratiofile = uigetfile('.mat', 'Select ratio data');
load(ratiofile);

pref = regexp(ratiofile, '.*(?=_Ratio)', 'match');
if ~isempty(pref)
    curvfile = uigetfile('.mat', 'Select curv data', [pref{1,1} '*.mat']);
else
    curvfile = uigetfile('.mat', 'Select curv data', '*.mat');
end
load(curvfile);

if exist('curvdata', 'var')
    curv_cor = curvdata;
    if ~isequal(size(ratio), size(curv_cor)) % this happens when curvature is calculated by angles between neurons
        ratio = ratio(:,2:end-1); % angle data do not exist for the first and last neurons in this case
    end
elseif exist('curv_rev', 'var')
    curv_cor = curv_rev;
    ratio = ratio_rev;
end

% Flip curvature data if necessary!!!!!!
curvflip = questdlg('Curvature data flipped?', 'Flip', 'Yes', 'No', 'No');
if isequal(curvflip, 'Yes')
    curv_cor = -curv_cor;
end

ratioall = ratio;
curvall = curv_cor;

%% Select ratio and curvature data files for DBs

ratiofile = uigetfile('.mat', 'Select ratio data');
load(ratiofile);
ratioall = signal{1,1}./signal_mirror{1,1};

pref = regexp(ratiofile, '.*(?=_\w*d*)', 'match');
n = regexp(ratiofile, '\d*(?=[.]mat)', 'match'); 
neuronnm = str2double(n{1,1});
curvall = zeros(size(ratioall,1),3);
for iii = 1:3 
    curvfile = uigetfile('.mat', ...
        ['Select curv data #: ' num2str(iii-2+neuronnm)], [pref{1,1} '*_curated.mat']);
    if curvfile==0
        break
    else
        load(curvfile);
        curvall(:,iii) = curvdatafiltered;
    end
end

if ~isequal(size(curvdatafiltered,1), size(neuron_position_data,1))
    fprintf('Calcium signal and curvature data mismatch at size. \n');
end

%% Pruning frames by manual input

p1 = [423:445 603:799 857:872 953:958];
p2 = [1011:1021];
p3 = [];
p4 = [];
deleted = union_several(p1,p2,p3,p4)';

frphasepre = [78:94 642:654 676:730 757:771 792:934 995:1109];
allphase = 1:length(ratioall);
frphasepre = setdiff(allphase, frphasepre);

frphasepruned = setdiff(frphasepre, intersect(frphasepre, deleted));

fprintf('frames manually pruned.\n');

%% Prune frames by loading existing data

prefsave = regexp(ratiofile, '.*(?=[.]mat)', 'match');
loadfilename = [prefsave{1} '_xcorr.mat'];
if exist(loadfilename, 'file')
    load(loadfilename);
    fprintf('existing pruned frames loaded.\n');
else
    fprintf('pruned frames have not been generated yet.\n');
end

%% Calculate cross-correlation in either forward or backward phase

close all;
[correall, lagsall, lagtimeall] = ...
    xcorr_pairs_special(ratioall, curvall, frphasepruned, 0.1);
fprintf('cross-correlation calculated.\n');
fprintf('do not forget to save data......\n');

%% Save data

prefsave = regexp(ratiofile, '.*(?=[.]mat)', 'match');
save([prefsave{1} '_xcorr.mat'], 'correall', 'lagsall', 'lagtimeall', ...
    'deleted', 'frphasepre', 'frphasepruned');
fprintf([prefsave{1} '_xcorr.mat saved.\n']);
savefig(gcf, [prefsave{1} '_xcorr.fig']); 
fprintf([prefsave{1} '_xcorr.fig saved.\n']);
