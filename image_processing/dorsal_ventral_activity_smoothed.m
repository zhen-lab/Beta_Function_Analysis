function [] = dorsal_ventral_activity_smoothed(dorsal_smd, ventral_smd, fps, caxis_limit)

numframes = length(dorsal_smd);
timefilter = 2; bodyfilter = 10;
h = fspecial('average', [timefilter bodyfilter]); 

dorsal_smd_fil = imfilter(dorsal_smd,  h , 'replicate');
dorsal_smd_fil = smooth(dorsal_smd_fil);
dorsal_smd_fil = reshape(dorsal_smd_fil, [], numframes);
figure(1); imagesc(dorsal_smd_fil); colorbar; caxis([0 caxis_limit]);
title('Acitivity Diagram');
set(gca, 'YDir', 'normal'); % To display the head-tail direction from bottom to top
set(gca,'YTICK',[1 20 40 60 80 100]);
set(gca,'YTICKLABEL',[0 0.2 0.4 0.6 0.8 1]);
set(gca,'XTICK',1:10*fps:numframes);
x_tick=get(gca,'XTICK');
set(gca,'XTICKLABEL',(x_tick-1)/fps);
ylabel('Fractional distance along the centerline/head=0,tail=1');
xlabel('Time/s');
colormap(hot);

ventral_smd_fil = imfilter(ventral_smd,  h , 'replicate');
ventral_smd_fil = smooth(ventral_smd_fil);
ventral_smd_fil = reshape(ventral_smd_fil, [], numframes);
figure(2); imagesc(ventral_smd_fil); colorbar; caxis([0 caxis_limit]);
title('Acitivity Diagram');
set(gca, 'YDir', 'normal'); % To display the head-tail direction from bottom to top
set(gca,'YTICK',[1 20 40 60 80 100]);
set(gca,'YTICKLABEL',[0 0.2 0.4 0.6 0.8 1]);
set(gca,'XTICK',1:10*fps:numframes);
x_tick=get(gca,'XTICK');
set(gca,'XTICKLABEL',(x_tick-1)/fps);
ylabel('Fractional distance along the centerline/head=0,tail=1');
xlabel('Time/s');
colormap(hot);

end