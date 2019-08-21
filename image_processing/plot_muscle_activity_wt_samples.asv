fig = gcf; 
dataObjs = findobj(fig,'-property','CData');
d = dataObjs(1).CData;
d = -d;
imagesc(d);
colormap(plasma)
colorbar
set(gca, 'xticklabel', [], 'yticklabel', [], 'ticklength', [0 0]);
caxis([0 2])
set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 0.5]);