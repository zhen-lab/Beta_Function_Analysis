%% Load data
load('traces.mat');
%% Cluster by kmeans
nm = 8;
idx = kmeans(gcamp', nm);
% Map clusters to neurons numbers
ord = 1:length(idx);
num = [idx ord'];
numsorted = sortrows(num);
gcamp_smoothed = smooth(gcamp);
gcamp_smoothed = reshape(gcamp_smoothed, size(gcamp,1),[]);
%% Draw figure
fig = imagesc(gcamp_smoothed(:,numsorted(:,2))');
colormap jet; 
c = colorbar;
caxis([0 200]);
set(c, 'location', 'eastoutside', 'ticklabels', []);
ticks = 0:100:500;
set(gca, 'xtick', ticks, 'xticklabel', ticks/5);
%% Save data and figures
save('kmeans', 'idx', 'numsorted');
[~, fname, ~] = fileparts(pwd);
savefig(fname);
saveas(fig, fname, 'tif')