listing = dir('*.mat');
figure(1); hold on;
for i = 1:length(listing)
    rawdata = matfile(listing(i).name);
    gfp = rawdata.signal;
    rfp = rawdata.signal_mirror;
    plot(smooth(gfp{1,1}./rfp{1,1}), 'linewidth', 1.5, 'color', 0.8*[1 0.2 0]);
end
% hold off;