function [] = image_format(data, curv_or_activity, not_smd, maxi, dorsal_smd, ventral_smd)

if not_smd == 1
    numframes = size(data, 2);
    timefilter = 2; bodyfilter = 10;
    h = fspecial('average', [timefilter bodyfilter]);
    data = imfilter(data,  h , 'replicate');
    data = smooth(data);
    data = reshape(data, [], numframes);
end;

    imagesc(data); colorbar; 

if curv_or_activity == 1
    caxis([-10 10]);
%     title('Cuvature Diagram');
else
    caxis([0 maxi]);
%     title('Activity Diagram');
    colormap(hot);
end;
    
% fps =20; numframes = size(data, 2);
    
set(gca, 'YDir', 'normal'); % To display the head-tail direction from bottom to top
set(gca,'YTICK',[1 20 40 60 80 100]);
set(gca,'YTICKLABEL',[]);

% set(gca,'XTICK',1:10*fps:numframes);
% x_tick=get(gca,'XTICK');
set(gca,'XTICKLABEL',[]);

% ylabel('Fractional distance along the centerline/head=0,tail=1');
% xlabel('Time/s');

 [mean(mean(dorsal_smd(35:end, :))); mean(mean(ventral_smd(35:end, :)))]

end