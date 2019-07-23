%% Import tiff files
setup_proof_reading;
fprintf('tiff files loading finished. \n');

%% Save ratiometric images

thresh = 200;

for i = 1:length(imagelist)
    image = imagelist_g{i,1}./imagelist_r{i,1};
    shape = imagelist_r{i,1}>thresh;
    imageratio = image.*uint16(shape);
    if i==1
        imwrite(imageratio, [filename(1:end-4) '_ratio.tif']);
    else
        imwrite(imageratio, [filename(1:end-4) '_ratio.tif'], 'writemode', 'append');
    end
end

%%
close all;
mini = 0; maxi = 150; 
r = 10;
% backperiod = [109 193];
ftl = 68;
fps = 1800/ftl; ff = 0.3;
% neuronnames = {'AVB', 'AVA/AVE'}; neuronnum = size(neuronnames, 2);
lightstart = 1; lightend = 80;
cmp = colormap(hot); cmp(1,:) = [0 0 0];
cutoff = mean2(imagelist{1,1});
rotangle = -70;

for frm = lightstart:lightend
    
    figure(1);
%     set(gcf,'units','normalized','outerposition',[0 0 1 1], 'color', 'k');
    set(gcf, 'color', 'w', 'pos', [10 10 500 600]);
    hold off;
%     imagedraw = imagelist{frm,1};
    imagedraw = imrotate((imagelist{frm,1}), rotangle, 'bilinear', 'crop');
    imagelogic = imagedraw>cutoff;
    imagedraw = double(imagedraw).*double(imagelogic);
    imagesc(imagedraw); % Have not added borders!!!!
    hold on;
    text(10, 10, [num2str((frm-lightstart)/(1800/ftl), '%.2f') 's' ], 'Color', 'w', 'fontsize', 10);
    
    axis equal;
    colormap(cmp);
    c = colorbar('southoutside');
    set(c, 'xcolor', 'none', 'ycolor', 'none', 'ticklabels', '');
    cpos = get(c, 'Position');
%     cpos(1) = 0.41; cpos(3) = cpos(3)/3.5; 
    cpos(1) = 0.41; cpos(3) = cpos(3)/3;
    cpos(4) = cpos(4)/3;    
    set(c, 'Position', cpos);        
    caxis([mini maxi]);
    set(gca, 'visible', 'off');
    
    hold on
%     for nn = 1:neuronnum
%         plot(dual_position_data{frm,nn}(1)-n,dual_position_data{frm,nn}(2),'ow', 'markersize', 4*r);
%         text(dual_position_data{frm,nn}(1)-n-2*r, dual_position_data{frm,nn}(2)+2*r, neuronnames(nn), 'color', 'w', 'fontsize', 10);
%     end
      
%     text(10, 40, [num2str((frm-lightstart)/fps, '%.2f') 's' ], 'Color', 'w', 'fontsize', 10);
    
%     if frm<backperiod(1) || frm>backperiod(2)
%         delete(findall(gcf,'type','annotation'));
%         annotation('arrow', [0.41 0.39], [0.5 0.5], 'linewidth', 10, 'headwidth', 20, 'color', 'w');
%     else
%         delete(findall(gcf,'type','annotation'));
%         annotation('arrow', [0.61 0.63], [0.5 0.5], 'linewidth', 10, 'headwidth', 20, 'color', 'w');
%     end
    
    % Create gif file
    drawnow;
    frame = getframe(1);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);
    if frm == lightstart
      imwrite(imind, cm, [filename '.gif'], 'gif', 'LoopCount', inf, 'DelayTime', ff/fps);
    else
      imwrite(imind, cm, [filename '.gif'], 'gif', 'WriteMode', 'append', 'DelayTime', ff/fps);
    end

end