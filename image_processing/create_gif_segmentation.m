ff = 1;
fps = 1800/68;

for i = 1:size(img_stack,3)   % Create gif file
    
    figure(1);
    hold off;
    imshow(imagelist{i});    
    caxis([0 3000]);
    hold on;
    plot(dorsal_data{2*i-1,1}, dorsal_data{2*i,1}, 'r');
    plot(ventral_data{2*i-1,1}, ventral_data{2*i,1}, 'b');
    plot(centerline_data_spline(1,2*i-1), centerline_data_spline(1, 2*i), 'og', 'markersize', 4, 'markerfacecolor', 'g');
    plot(centerline_data_spline(end,2*i-1), centerline_data_spline(end, 2*i), 'oy', 'markersize', 4, 'markerfacecolor', 'y');
    plot(centerline_data_spline(:, 2*i-1), centerline_data_spline(:, 2*i), 'w', 'linewidth', 1.5);
    
    drawnow;
    frame = getframe(1);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);
    if i == 1
      imwrite(imind, cm, [filename '_segmented.gif'], 'gif', 'LoopCount', inf, 'DelayTime', ff/fps);
    else
      imwrite(imind, cm, [filename '_segmented.gif'], 'gif', 'WriteMode', 'append', 'DelayTime', ff/fps);
    end
    
 end