function [avgint] = create_gif_behavior(filename, imagelist, close_para, fill_para, thresh, adpt_max, adpt_min)
    
%     fps = 20;
%     spline_p = 0.001;
    istart = 1;
    iend = size(imagelist,1);
%     sigma_gauss = 1.2; % Gaussian filter for small area
    spf = 1/26;
    filename = [filename(1:end-4) '.gif'];
%     thresh = 1;
%     bgd = 300;
%     do_multitif = 1;


% numframes = iend - istart + 1; % Number of frames analyzed
% numcurvpts = 100; % Number of segments for a worm

% ventral_data = cell(2*numframes, 1);
% dorsal_data = cell(2*numframes, 1);
% centerline_data = zeros(numcurvpts, 2*numframes); 
% centerline_data_spline = zeros(numcurvpts, 2*numframes); % Stores position info for all segments in all frames
% curvdata = zeros(numframes,numcurvpts);

% Determine upper and lower bounds of image intensities
fprintf('initiating calculation of upper and lower bounds \n');
avgint = zeros(length(imagelist), 1);
for i = 1:length(imagelist)
    avgint(i) = mean(mean(imagelist{i,1}));
end
avgint_min = min(avgint);
avgint_max = max(avgint);
fprintf('finished calculation of upper and lower bounds \n');
%%%%%%%%%%%%%%%%%%%%%%%Main Loop%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for j = istart:iend

    img = imagelist{j,1};
    
    % Chooose adaptive threshold based on average image intensity
    avgint_img = mean(mean(imagelist{j,1}));
    [~, max_or_min] = min(abs([avgint_img; avgint_img]-[avgint_max; avgint_min]));
    adpt = (max_or_min==1)*adpt_min + (max_or_min==2)*adpt_max;
    
    % Adjust contrast and filter out noise for GFP+RFP image
    adj = imbinarize(img, 'adaptive', 'ForegroundPolarity', 'dark', 'Sensitivity', adpt);

    % Fill out data matrices with -1 for "no data"
%     corr_pts = -ones(numcurvpts, 6);
%     corr_data = -ones(numcurvpts, 3);
    
    %%%%%%%%%%Image processing%%%%%%%%%%%%%%%%%%%%%%%%

    sigma = sqrt(2); % Standard deviation of the Gaussian filter
    [~, threshold] = edge(adj, 'canny'); % Canny method for edge detection with automatic threshold, [low high]
    BW_edge = edge(adj, 'canny', thresh*threshold, sigma); 
    BW_edge = bwareaopen(BW_edge, 4);
    
    se1 = strel('disk', close_para);
    closeBW_edge = imclose(BW_edge, se1); 

    bw = imfill(closeBW_edge, 'holes');
    se2 = strel('disk', fill_para);
    bw = imopen(bw, se2);
    imagesc(bw); colormap(flipud(gray));
    set(gca, 'visible', 'off');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    drawnow;
    frame = getframe(1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if j == 1
      imwrite(imind,cm,filename,'gif', 'Loopcount',inf, 'delaytime', spf);
    else
      imwrite(imind,cm,filename,'gif','WriteMode','append', 'delaytime', spf);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
end % End of calculations for 1 frame

hold off;

end % End of everything