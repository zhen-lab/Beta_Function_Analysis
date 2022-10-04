
numtraces = size(dorsal,1);
nrow = 5;
ncol = ceil(numtraces/nrow); 
half = zeros(numtraces, 2);

close all;
f = figure;

for i = 1:numtraces
   
    d = dorsal(i,:);
    v = ventral(i,:);
    d(isinf(d)) = NaN;
    v(isinf(v)) = NaN;
    
    heightd = (max(d)+min(d))/2;
    heightv = (max(v)+min(v))/2;
    
    halfd = find(d>=heightd, 1, 'first'); 
    halfv = find(v>=heightv, 1, 'first'); 
    half(i,:) = [halfd halfv];
    
    subplot(ncol,nrow,i);    
    hold on; 
    plot(d, 'm', 'linewidth', 3); 
    plot(v, 'g', 'linewidth', 3); 
    plot(halfd, heightd, 'om', 'markersize', 6, 'markerfacecolor', 'm');
    plot(halfv, heightv, 'og', 'markersize', 6, 'markerfacecolor', 'g');
    set(gca, 'visible', 'off'); 

end

% set(gcf, 'color', 'w', 'units', 'normalized',...
%     'outerposition', [0 0 1 1]);

parts = strsplit(pwd, '\');
folder = regexp(parts, 'Alpha\w*', 'match');
idx = find(~cellfun(@isempty, folder));
filename = parts(idx+1); fname = filename{1};
save(['Halftime_' fname '.mat'], 'half');
savefig(f, ['Halftime_' fname '.fig']);
saveas(f, ['Halftime_' fname '.tif'], 'tiffn');