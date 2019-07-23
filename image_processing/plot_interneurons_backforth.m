figure; 
hold on; 
maxi = 4; 
lwidth = 3;

for i = 1:size(backframes,1)
    st = backframes(i, 1); en = backframes(i, 2);
    h = fill([st en en st], [0 0 maxi maxi], 0.8*[1 1 1], 'edgecolor', 'none');
end

clr = lbmap(10, 'RedBlue');

load('temp36_AVA.mat')
plot(signal{1,1}./signal_mirror{1,1}, 'linewidth', lwidth, 'color', clr(2,:));
load('temp36_AVB.mat')
plot(signal{1,1}./signal_mirror{1,1}, 'linewidth', lwidth, 'color', clr(end-1,:));
ylim([0 maxi]);
set(gca, 'visible', 'off')