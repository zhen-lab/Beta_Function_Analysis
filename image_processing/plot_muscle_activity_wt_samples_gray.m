%%
colormap(flipud(gray));
%%
pbaspect([1.5 1 1]);
caxis([0.5 2]);
c = colorbar;
cticks = get(c, 'ticks');
labels = arrayfun(@(x) sprintf('%.1f',x), cticks, 'Un', 0);
set(c, 'ticklabels', labels);