function [dorsal_data, ventral_data, centerline_data, centerline_data_spline, curvdata, curvdatafiltered] = extract_centerline_vd(imagelist, close_para, fill_para, thresh)
    
    fps = 20;
    spline_p = 0.001;
    istart = 1;
    iend = size(imagelist,1);
%     thresh = 10;
%     do_multitif = 1;

    % Enter the frame rate, start and end frame selected from the imaging sequences
    
%     if exist('do_multitif', 'var')
%         answer = inputdlg({'Fps', 'Start frame', 'End frame', 'Spline parameter', 'Threshold(1-3)', 'Multi-Tiff files'}, 'Cancel to clear previous', 1, ...
%             {num2str(fps),num2str(istart),num2str(iend),num2str(spline_p),num2str(thresh),num2str(do_multitif)});
%     else
%         answer = inputdlg({'Fps', 'Start frame', 'End frame', 'Spline parameter', 'Threshold(1-3)', 'Multi-Tiff files'}, '', 1, ...
%             {num2str(fps),num2str(istart),num2str(iend),num2str(spline_p),num2str(thresh),'1'});
%     end
%     
%     
%     if isempty(answer) 
%         answer = inputdlg({'Fps', 'Start frame', 'End frame', 'Spline parameter', 'Threshold(1-3)', 'Multi-Tiff files'}, '', 1, ...
%             {num2str(fps),num2str(istart),num2str(iend),num2str(spline_p),num2str(thresh),'1'});
%     end
    
%     fps = str2double(answer{1});
%     istart = str2double(answer{2});
%     iend = str2double(answer{3});
%     spline_p = str2double(answer{4});
%     thresh = str2double(answer{5});
%     do_multitif = str2double(answer{6});
    
% end

numframes = iend - istart + 1; % Number of frames analyzed
numcurvpts = 100; % Number of segments for a worm

ventral_data = cell(2*numframes, 1);
dorsal_data = cell(2*numframes, 1);
centerline_data = zeros(numcurvpts, 2*numframes); 
centerline_data_spline = zeros(numcurvpts, 2*numframes); % Stores position info for all segments in all frames
curvdata = zeros(numframes,numcurvpts);
angle_data = zeros(numframes, numcurvpts+1);
%vulva_indices = zeros(numframes,1);

boundary_segment1 = zeros(numframes,1);
boundary_segment2 = zeros(numframes,1);
bcols = [];
%max_pt = 100;
%min_pt = 1;

%%%%%%%%%%%%%%%%%%%%%%%Main Loop%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for j=istart:iend

    % GFP channel and RFP channel separation  
%     rb = imagelist_g{j, 1}; % GFP channel
%     c = imagelist_r{j, 1}; % RFP channel
    
    % Adjust contrast and filter out noise for GFP+RFP image
%     adj = 2*im2bw(imagelist{j, 1});
%     adj = imadjust(adj, [0.006; 0.9], [], 0.1);
%     adj = wiener2(adj, [5 5]);
%     adj = imclearborder(adj);
%     adj = imcomplement(adj);
    adj = imagelist{j, 1};
    
    % Fill out data matrices with -1 for "no data"
    corr_pts = -ones(numcurvpts, 6);
    corr_data = -ones(numcurvpts, 3);
%     corr_data_average = -ones(numcurvpts, 3);

    % Head and tail manual estimate in the starting frame
    if j == istart
        
        figure (1);
        imagesc(adj); colormap(gray); % Adjust low contrast gray image
        axis equal; axis off; hold on;
        
        text(10, 20, 'Click on head', 'Color', 'white');
        [headx, heady] = ginput(1);
        hold on; plot(headx, heady, 'og');
%         headx0 = headx; heady0 = heady;

        text(10, 40, 'Click on tail', 'Color', 'white');
        [tailx, taily] = ginput(1);
        hold on; plot(tailx, taily, 'oy');
%         tailx0 = tailx; taily0 = taily;
        
        
%         text(10, 60, 'Select 1st background point', 'Color', 'white');
%         [cropx1_bkg, cropy1_bkg] = ginput(1);
%         cropx1_bkg= floor(cropx1_bkg);
%         cropy1_bkg  = floor(cropy1_bkg);
%         hold on, plot([1 n], [cropy1_bkg cropy1_bkg], '-w');
%         hold on, plot([cropx1_bkg cropx1_bkg], [1 m], '-w');
%         
%         text(10, 80, 'Select 2nd background point', 'Color','white');
%         [cropx2_bkg, cropy2_bkg] = ginput(1);
%         cropx2_bkg = floor(cropx2_bkg);
%         cropy2_bkg = floor(cropy2_bkg);
%         hold on, plot([1 n], [cropy2_bkg cropy2_bkg], '-w');
%         hold on, plot([cropx2_bkg cropx2_bkg], [1 m/2], '-w');

%         text(10,100,'click two points separated by worm diameter','Color', 'white');
%         tmp1 = ginput(1); 
%         plot(tmp1(1),tmp1(2), 'ow');
%         
%         text(10,120,'click two points separated by worm diameter','Color', 'white');
%         tmp2 = ginput(1);
%         plot(tmp2(1),tmp2(2), 'ow');
%         
%         worm_diam = norm(tmp1-tmp2);
%         
%         filsize=0.25;
%         
%         if mod(round(filsize*worm_diam),2)==1
%            filradius = round(filsize*worm_diam/2);
%         else
%            filradius = round(filsize*worm_diam/2)+1;
%         end
%      
%         fil = fspecial('disk', filradius);
        
%         text(10,120,'define threshold value for computing the skeleton', 'Color', 'black');
%         text(10,10,'define threshold value for computing the muscle intensity', 'Color', 'white');     %user defines threshold value for finding worm boundary (this is global)
%         [threshx, threshy] = ginput(1);
%         lvl_original = mean(mean(I(round(threshy)-2:round(threshy)+2, round(threshx)-2:round(threshx)+2))); 
%         text(10,10,'define threshold value for computing the muscle intensity', 'Color', 'black');
    
        
%         bkg_rb=mean(mean(rb(cropy1_bkg:cropy2_bkg,cropx1_bkg:cropx2_bkg)));
%         bkg_c=mean(mean(c(cropy1_bkg:cropy2_bkg,cropx1_bkg:cropx2_bkg))); 
    
        
        text(10, 60, 'Click on ventral side', 'Color', 'white');     
        [vulvax, vulvay] = ginput(1);
        
    end
    
    figure (1);
    drawnow; hold off;
    imagesc(adj); colormap(gray); axis equal; axis off; hold on;
    text(10,  10,  num2str(j),  'Color',  'white');
%     text(10,  10,  num2str(j + istart - 1),  'Color',  'white');
    
    
    %%%%%%%%%%Image processing%%%%%%%%%%%%%%%%%%%%%%%%

%     rb_cross_adj=rb-bkg_rb;
%     I_cross_adj=c-bkg_c; % -0.153*rb_cross_adj;
%     ratio_cross=ratio(cropy1:cropy2,cropx1:cropx2);

    sigma = sqrt(2); % Standard deviation of the Gaussian filter
    
    [~, threshold] = edge(adj, 'canny'); % Canny method for edge detection with automatic threshold, [low high]
    BW_edge = edge(adj, 'canny', thresh*threshold, sigma); 
    BW_edge = bwareaopen(BW_edge, 4);
    
    se1 = strel('disk', close_para);
    closeBW_edge = imclose(BW_edge, se1); 

    bw = imfill(closeBW_edge, 'holes');
    se2 = strel('disk', fill_para);
    bw = imopen(bw, se2);
   
    STATS = regionprops(logical(bw), 'Area', 'Centroid');
    if size(STATS,1) == 0    
        disp('Error: no worm found');
        break;
    end    
    num=1;
    area=0;
    if (length(STATS)>1)
        for k=1:length(STATS)            
            if STATS(k).Area>area                
                area=STATS(k).Area;
                num=k;
            end            
        end        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [B, ~] = bwboundaries(bw, 'noholes');
    
    B1 = B{num};
    B1_size = size(B1, 1);
    ksep = ceil(B1_size/20);

    B1_plus = circshift(B1, [ksep 0]);
    B1_minus = circshift(B1, [-ksep 0]);
    
    % AA and BB are vectors between a point on boundary and neighbors +/- ksep away
    
    AA = B1 - B1_plus;  
    BB = B1 - B1_minus;

    cAA = AA(:, 1) + sqrt(-1) * AA(:, 2);
    cBB = BB(:, 1) + sqrt(-1) * BB(:, 2);

    B1_angle1 = unwrap(angle(cBB ./ cAA));
    
    %%%%%%%%%%%Find head and tail%%%%%%%%%%%%%%
    
    % First minimum
    min1 = find(B1_angle1 == min(B1_angle1),1); % Find point on boundary w/ minimum angle between AA, BB
    B1_angle2 = circshift(B1_angle1, -min1); % Shift to have tail as the starting point
    
    % Second minimum
    min2_tmp = round(.25*B1_size)-1+find(B1_angle2(round(.25*B1_size):round(0.75*B1_size))==min(B1_angle2(round(.25*B1_size):round(0.75*B1_size))),1);  % Find minimum in other half
    min2 = 1 + mod(min2_tmp + min1-1, B1_size);
    
    tmp = circshift(B1, [1-min1 0]);
    end1 = 1 + mod(min2 - min1-1, B1_size);
    path1 = tmp(1:end1,:);
    path2 = tmp(end:-1:end1,:);

    if norm(path1(1,:) - [heady headx]) > norm(path1(end,:) - [heady headx]) % if min1 is at tail, reverse both paths
        tmp = path1;
        path1 = path2(end:-1:1,:);
        path2 = tmp(end:-1:1,:);
    end
    
    heady = path1(1,1);
    headx = path1(1,2);
%     taily = path1(end,1);
%     tailx = path1(end,2);
    
    %%%%%%%%%%%Generate skeleton%%%%%%%%%%%%%%
    
    path_length = numcurvpts;

    path1_rescaled = zeros(path_length,2);
    path2_rescaled = zeros(path_length,2);
    path1_rescaled2 = zeros(path_length,2);
    path2_rescaled2 = zeros(path_length,2);
    
    % Interp the paths to the number of segments we want (default 100)
    path1_rescaled(:,1) = interp1(0:size(path1,1)-1, path1(:,1), (size(path1,1)-1)*(0:path_length-1)/(path_length-1), 'linear');
    path1_rescaled(:,2) = interp1(0:size(path1,1)-1, path1(:,2), (size(path1,1)-1)*(0:path_length-1)/(path_length-1), 'linear');
    path2_rescaled(:,1) = interp1(0:size(path2,1)-1, path2(:,1), (size(path2,1)-1)*(0:path_length-1)/(path_length-1), 'linear');
    path2_rescaled(:,2) = interp1(0:size(path2,1)-1, path2(:,2), (size(path2,1)-1)*(0:path_length-1)/(path_length-1), 'linear');

    for kk = 1:path_length
        path_temp=path2_rescaled(max(1,kk-15):min(path_length,kk+15),:);
        l=length(path_temp);
        tmp1 = repmat(path1_rescaled(kk,:), [l,1]) - path_temp;
        tmp2 = sqrt(tmp1(:,1).^2 + tmp1(:,2).^2);
        path2_rescaled2(kk,:) = path_temp(find(tmp2==min(tmp2),1),:);
    end
    
    for kk = 1:path_length
        path_temp=path1_rescaled(max(1,kk-15):min(path_length,kk+15),:);
        l=length(path_temp);
        tmp1 = repmat(path2_rescaled(kk,:), [l,1]) - path_temp;
        tmp2 = sqrt(tmp1(:,1).^2 + tmp1(:,2).^2);
        path1_rescaled2(kk,:) = path_temp(find(tmp2==min(tmp2),1),:);
    end
    
    weight_fn = ones(path_length,1);
    tmp = round(path_length*0.1);
    weight_fn(1:tmp)=(0:tmp-1)/tmp;
    weight_fn(end-tmp+1:end)=(tmp-1:-1:0)/tmp;
    weight_fn_new = [weight_fn weight_fn];
  
    midline2a = 0.5*(path1_rescaled+path2_rescaled2);
    midline2b = 0.5*(path1_rescaled2+path2_rescaled);
    midline_mixed = midline2a .* weight_fn_new + midline2b .* (1-weight_fn_new);
    [line, cv2, sp_curv_num] = spline_line(midline_mixed, spline_p, numcurvpts); % Interp and smooth the midline

%     figure(1);
%     drawnow;
    % Plot ventral (blue) and dorsal (red) edges, and skeleton(white)
    if j==istart
        d1_square=min(sum((path1_rescaled-repmat([vulvay vulvax],path_length,1)).^2,2));
        d2_square=min(sum((path2_rescaled-repmat([vulvay vulvax],path_length,1)).^2,2));
    end
        
    if d1_square<d2_square
        ventral=path1_rescaled;
        dorsal=path2_rescaled;
    else
        ventral=path2_rescaled;
        dorsal=path1_rescaled;
    end
 
    %d_p=csaps(dorsal(:,2),dorsal(:,1),spline_p);
    %v_p=csaps(ventral(:,2),ventral(:,1),spline_p);    
    %seg_length=sp_curv_num(end)/path_length;
    %d_path = spline_line(dorsal,spline_p,2*numcurvpts);
    %v_path = spline_line(ventral,spline_p,2*numcurvpts);
    hold on;
    plot(dorsal(:,2), dorsal(:,1), '-r', 'LineWidth', 1); 
    plot(ventral(:,2), ventral(:,1), '-b', 'LineWidth', 1); 
    plot(midline_mixed(:,2), midline_mixed(:,1), '-w', 'LineWidth', 1);
    
    % Plot head and tail
    plot(path1_rescaled(1,2), path1_rescaled(1,1), 'o', 'markerfacecolor', 'g');
    plot(path2_rescaled(end,2), path2_rescaled(end,1), 'o', 'markerfacecolor', 'y'); 
    
    % Store ventral and dorsal edges, and skeleton information in all frames
    ventral_data{2*(j-istart+1)-1} = ventral(:, 2); 
    ventral_data{2*(j-istart+1)} = ventral(:, 1);
    dorsal_data{2*(j-istart+1)-1} = dorsal(:, 2);
    dorsal_data{2*(j-istart+1)} = dorsal(:, 1);
    centerline_data(:, 2*(j-istart+1)-1) = midline_mixed(:, 1);
    centerline_data(:, 2*(j-istart+1)) = midline_mixed(:, 2);
    centerline_data_spline(:, 2*(j-istart+1)-1) = line(:, 1);
    centerline_data_spline(:, 2*(j-istart+1)) = line(:, 2);
    
    % Interpolate to equally spaced length units
    cv2i = interp1(sp_curv_num+.00001*(0:length(sp_curv_num)-1), cv2, (0:(sp_curv_num(end)-1)/(numcurvpts+1):(sp_curv_num(end)-1)));
    df2 = diff(cv2i, 1, 1);
    atdf2 = unwrap(atan2(-df2(:,2), df2(:,1)));
    angle_data((j-istart+1), :) = atdf2';
    
    % Collect curvature information for midline
    curv = unwrap(diff(atdf2, 1)); 
    corr_data(:, 1) = curv;
    
    % Save curvature data
    corr_data(end,2:3)=corr_data(end-1,2:3);
    corr_data(1,2:3)=corr_data(2,2:3);
    corr_pts_all = zeros(numframes, size(corr_pts, 1), size(corr_pts, 2)); % Frame number, segment number, 3 points
    corr_pts_all((j-istart+1), :, :) = corr_pts;
    curvdata((j-istart+1), :) = corr_data(:, 1);
    %ventral_brightness_data(j,:) = corr_data(:,3);
    %dorsal_brightness_data (j,:) = corr_data(:,2);
    %ratio_data(j,:,:)=ratio;
    
    %pause;
  
end % End of calculations for 1 frame

hold off;

% end % End of calculations for all frames

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%Delete invalid frames%%%%%%%

curvdata = curvdata'; % After transposition, columns are frame numbers, rows are segment numbers

figure(2); clf;

    imagesc(curvdata(:, :)*100); caxis([-10, 10]);
    %imagesc(dorsal_brightness_data_filtered(:,:)); colormap(jet); colorbar; caxis([0.8 5]);
    title('Click on invalid columns, and then press return');
    badcols = ginput; % The first column stores the frame numbers, the second stores the segment numbers
  
    if ~isempty(badcols)
        bcols = [bcols'; ceil(badcols(:,1))]; % The first column stores the frame number, the second stores the segment number
        bcols = unique(bcols); % Delete repeated elements
    end
    
     for j = 1:length(bcols)
         
         col = bcols(j);
         disp(strcat('Invalid col #', num2str(j), '=', num2str(col)));
         
         % Replace first frame if it is invalid
         if col == 1
             curvdata(:, 1) = curvdata(:, 2);
             corr_pts_all(1, :, :)=corr_pts_all(2, :, :);
%              ventral_brightness_data(1,:)=ventral_brightness_data(2,:);
%              dorsal_brightness_data(1,:)=ventral_brightness_data(2,:);
         end
         
         % Replace last frame if it is invalid
         if col == numframes
             curvdata(:, numframes) = curvdata(:, numframes-1);
             corr_pts_all(numframes, :, :) = corr_pts_all(numframes-1, :, :);
%              ventral_brightness_data(numframes,:)=ventral_brightness_data(numframes-1,:);
%              dorsal_brightness_data(numframes,:)=ventral_brightness_data(numframes-1,:);
         end
         
         % Replace invalid frame in between with its flanking frames
         if (col > 1) && (col < numframes)
             %m=size(badcols,1);
             %pre=col-j;
             %post=col+m-j+1;
             %if pre<1 pre=1;
             %end
             %if post>numframes
             %    post=numframes; 
             %end
             pre = col - 1;
             while ~isempty(find(bcols == pre, 1)) && (pre > 1)
                 pre = pre-1;
             end
             
             post = col + 1;
             while ~isempty(find(bcols == post, 1)) && (post < numframes)
                 post = post+1;
             end
             
             curvdata(:, col) = 0.5 * (curvdata(:, pre) + curvdata(:, post)); % Mean of previous and post frames
             boundary_segment1(col) = 0.5 * (boundary_segment1(pre) + boundary_segment1(post));
             boundary_segment2(col) = 0.5 * (boundary_segment2(pre) + boundary_segment2(post));
             corr_pts_all(col, :, :) = 0.5 * (corr_pts_all(pre, :, :) + corr_pts_all(post, :, :));
%              ventral_brightness_data(col,:)=0.5*(ventral_brightness_data(pre,:)+ventral_brightness_data(post,:));
%              dorsal_brightness_data(col,:)=0.5*(dorsal_brightness_data(pre,:)+dorsal_brightness_data(post,:));
             
         end
         
     end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%Filter curvature by averaging%%%%%%%

% answer = inputdlg({'Time filter', 'Body coord filter', 'Mean=0, Median=1'}, '', 1, {num2str(2), num2str(10), '0'});
% timefilter = str2double(answer{1});
% bodyfilter = str2double(answer{2});

timefilter = 2; bodyfilter = 10;
h = fspecial('average', [timefilter bodyfilter]);
curvdatafiltered = imfilter(curvdata*100,  h , 'replicate');
curvdatafiltered = smooth(curvdatafiltered);
curvdatafiltered = reshape(curvdatafiltered, [], numframes);
% ventral_brightness_data_filtered = imfilter(ventral_brightness_data, h, 'replicate');
% dorsal_brightness_data_filtered = imfilter(dorsal_brightness_data, h, 'replicate');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%2D illustration of curvature%%%%%%
figure(3); imagesc(curvdatafiltered); colorbar; caxis([-10 10]);

title('Cuvature Diagram');

set(gca, 'YDir', 'normal'); % To display the head-tail direction from bottom to top
set(gca,'YTICK',[1 20 40 60 80 100]);
set(gca,'YTICKLABEL',[0 0.2 0.4 0.6 0.8 1]);

set(gca,'XTICK',1:10*fps:numframes);
x_tick=get(gca,'XTICK');
set(gca,'XTICKLABEL',(x_tick-1)/fps);

ylabel('Fractional distance along the centerline/head=0,tail=1');
xlabel('Time/s');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%3D illustration of curvature%%%%%%
% figure(4); mesh(curvdatafiltered);
% 
% title('Cuvature Diagram');
% 
% set(gca,'YTICK',[1 20 40 60 80 100]);
% set(gca,'YTICKLABEL',[0 0.2 0.4 0.6 0.8 1]);
% 
% set(gca,'XTICK',1:10*fps:numframes);
% x_tick=get(gca,'XTICK');
% set(gca,'XTICKLABEL',(x_tick-1)/fps);
% 
% ylabel('Fractional distance along the centerline/head=0,tail=1');
% xlabel('Time/s');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figure; imagesc(ventral_brightness_data_filtered(:,:)); colormap(jet); colorbar; caxis([3 15]);
% 
% title('ventral ratio GC3/RFP');
% set(gca,'XTICK',[1 20 40 60 80 100]);
% set(gca,'XTICKLABEL',[0 0.2 0.4 0.6 0.8 1]);
% 
% set(gca,'YTICK',1:10*fps:numframes);
% y_tick=get(gca,'YTICK');
% set(gca,'YTICKLABEL',(y_tick-1)/fps);
% 
% xlabel('fractional distance along the centerline (head=0; tail=1)');
% ylabel('time (s)');
% 
% figure; imagesc(dorsal_brightness_data_filtered(:,:)); colormap(jet); colorbar; caxis([3 15]);
% 
% title('dorsal ratio GC3/RFP');
% set(gca,'XTICK',[1 20 40 60 80 100]);
% set(gca,'XTICKLABEL',[0 0.2 0.4 0.6 0.8 1]);
% 
% set(gca,'YTICK',1:10*fps:numframes);
% y_tick=get(gca,'YTICK');
% set(gca,'YTICKLABEL',(y_tick-1)/fps);
% 
% xlabel('fractional distance along the centerline (head=0; tail=1)');
% ylabel('time (s)');

% if length(questdlg('Save this data? '))==3
%     [fn savepathname]= uiputfile('*.mat', 'choose file to save', strcat(fname, '_',num2str(istart),'-',num2str(iend),'.mat'));
%     if length(fn) > 1
%         fnamemat = strcat(savepathname,fn);
%         save(fnamemat);
%     end
% end
end % End of everything