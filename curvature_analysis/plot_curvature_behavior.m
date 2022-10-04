%% Clear everything

clear; close all;

%% Draw plot

filename = uigetfile('*.mat', 'Select curvature file');
body = 35;

if filename==0
    fprintf('user cancelled selection.\n');
else
    load(filename);
    if ~exist('curvdatafiltered','var')
        fprintf('no curvdature data exist.\n');
    else
        pref = regexp(filename, '.*(?=[.]mat)', 'match');    
        f = figure;
        imagesc(curvdatafiltered);
        set(gca, 'xticklabel', [], 'yticklabel', [], 'ticklength', [0 0]);
        colorbar;
        colormap(flipud(gray));
        hold on;
        plot([1 size(curvdatafiltered,2)],[body body],':w','linewidth',2);
        pbaspect([3 1 1]);
    end
end

%% Save plot

if ~exist('f','var')||~exist('pref','var')
    fprintf('please draw plot first.\n');
else
    filenamesave = [pref{1} '_curvature_gray'];
    savefig(f, pref{1});
    saveas(f, pref{1}, 'meta');
end
fprintf(['figure saved for ' pref{1} '.\n']);
