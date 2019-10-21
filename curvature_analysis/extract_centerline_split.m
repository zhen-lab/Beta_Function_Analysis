function [centerline_data, curvdata, curvdatafiltered, istart, iend] = extract_centerline_split(imagelist)

button = length(questdlg('Load new data?', '', 'Yes (TIFF)', 'Yes (MAT)', 'No', 'Yes (TIFF)') ) ;

% User chooses the type of file(s)

if button == 9
    
    [filename, pathname] = uigetfile({'*.mat'});
    load([pathname filename]);
    
elseif button == 10
    
    do_dialog = 1;

    if do_dialog

        [~, pathname] = uigetfile({'*.tif'});

        if exist('pathname', 'var') 
            try
                if isdir(pathname)
                cd(pathname);
                end
            catch err
                error('No such file exists');
            end
        end

        fps = 10;
        spline_p = 0.001;
        istart = 1;
        iend = size(imagelist,1);
        thresh = 0.2;
        do_multitif = 1;

        % Enter the frame rate, start and end frame selected from the imaging sequences
        if exist('do_multitif', 'var')
            answer = inputdlg({'Fps', 'Start frame', 'End frame', 'Spline parameter', 'Threshold(1-3)', 'Multi-Tiff files'}, 'Cancel to clear previous', 1, ...
                {num2str(fps),num2str(istart),num2str(iend),num2str(spline_p),num2str(thresh),num2str(do_multitif)});
        else
            answer = inputdlg({'Fps', 'Start frame', 'End frame', 'Spline parameter', 'Threshold(1-3)', 'Multi-Tiff files'}, '', 1, ...
                {num2str(fps),num2str(istart),num2str(iend),num2str(spline_p),num2str(thresh),'1'});
        end
        if isempty(answer) 
            answer = inputdlg({'Fps', 'Start frame', 'End frame', 'Spline parameter', 'Threshold(1-3)', 'Multi-Tiff files'}, '', 1, ...
                {num2str(fps),num2str(istart),num2str(iend),num2str(spline_p),num2str(thresh),'1'});
        end

        fps = str2double(answer{1});
        istart = str2double(answer{2});
        iend = str2double(answer{3});
        spline_p = str2double(answer{4});
        thresh = str2double(answer{5});
%         do_multitif = str2double(answer{6});
    
    end

numframes = iend - istart + 1; % Number of frames analyzed
numcurvpts = 100; % Number of segments for a worm
centerline_data = zeros(numcurvpts, 2*numframes); % Stores position info for all segments in all frames
curvdata = zeros(numframes,numcurvpts);
angle_data = zeros(numframes, numcurvpts+1);
bcols = [];

%%%%%%%%%%%%%%%%%%%%%%%Main Loop%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for j=istart:iend
    
    adj = imagelist{j, 1};
%     adj_rb_c = wiener2(adj, [5 5]);

    % Head and tail manual estimate in the starting frame
    if j == istart
        
        figure (1);
        imagesc(adj); colormap(gray); axis equal; axis off; hold on;
        
        text(10, 10, 'Click on head', 'Color', 'white');
        [headx, heady] = ginput(1);
        hold on; plot(headx, heady, 'or'); text(headx, heady, 'Head', 'Color', 'r');

        text(10, 20, 'Click on tail', 'Color', 'white');
        [tailx, taily] = ginput(1);
        hold on; plot(tailx, taily, 'oy'); text(tailx, taily, 'Tail', 'Color', 'y');

        text(10, 30,'Click on vulva', 'Color', 'white');     
        [vulvax, vulvay] = ginput(1);
        
%         text(10, 40, 'Select 1st background point', 'Color', 'white');
%         [cropx1_bkg, cropy1_bkg] = ginput(1);
%         cropx1_bkg = floor(cropx1_bkg); % Col number
%         cropy1_bkg  = floor(cropy1_bkg); % Row number
%         text(10, 50, 'Select 2nd background point', 'Color','white');
%         [cropx2_bkg, cropy2_bkg] = ginput(1);
%         cropx2_bkg = floor(cropx2_bkg); % Col number
%         cropy2_bkg = floor(cropy2_bkg); % Row number
%         % Use the mean of square between two background points
%         bkg_adj = mean(mean(adj(cropy1_bkg:cropy2_bkg, cropx1_bkg:cropx2_bkg)));
            
    end
    
%%%%%%%%%%Image processing%%%%%%%%%%%%%%%%%%%%%%%%
    
    figure (1), hold off;
    
    adj_rb_c = adj;

    imagesc(adj_rb_c); colormap(gray); axis equal; axis off; hold on;
    
    text(10,  10,  num2str(j),  'Color',  'white'); % Frame number

    sigma = sqrt(2); % Standard deviation of the Gaussian filter
    [~, threshold] = edge(adj_rb_c, 'canny'); % Canny method for edge detection with automatic threshold, [low high]
    BW_edge = edge(adj_rb_c, 'canny', thresh*threshold, sigma); 
    BW_edge = bwareaopen(BW_edge, 4);
    
    se1 = strel('disk', 16);
    closeBW_edge = imclose(BW_edge,se1); 
    bw = imfill(closeBW_edge,'holes');
    se2 = strel('disk', 8);
    bw = imopen(bw,se2);
   
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
    
    min1 = find(B1_angle1 == min(B1_angle1),1); % find point on boundary w/ minimum angle between AA, BB
    B1_angle2 = circshift(B1_angle1, -min1);
    min2_tmp = round(.25*B1_size)-1+find(B1_angle2(round(.25*B1_size):round(0.75*B1_size))==min(B1_angle2(round(.25*B1_size):round(0.75*B1_size))),1);  % find minimum in other half
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
    
    %%%%%%%%%%%Generate skeleton%%%%%%%%%%%%%%
    
    path_length = numcurvpts;

    path1_rescaled = zeros(path_length,2);
    path2_rescaled = zeros(path_length,2);
    path1_rescaled2 = zeros(path_length,2);
    path2_rescaled2 = zeros(path_length,2);
    
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

    figure(1);
    
    % Plot head and tail
    plot(path1_rescaled(1,2), path1_rescaled(1,1), 'om', 'MarkerFaceColor', 'm'); hold on;
    plot(path2_rescaled(end,2), path2_rescaled(end,1), 'oy', 'MarkerFaceColor', 'y'); hold on;

    if j==istart
        d1_square=min(sum((path1_rescaled-repmat([vulvay vulvax],path_length,1)).^2,2));
        d2_square=min(sum((path2_rescaled-repmat([vulvay vulvax],path_length,1)).^2,2));
    end
        
    if d1_square<d2_square
        ventral=path1;
        dorsal=path2;
    else
        ventral=path2;
        dorsal=path1;
    end
  
    % Plot ventral and dorsal edges
    plot(ventral(:,2), ventral(:,1), '-b', 'LineWidth', 1); hold on;
    plot(dorsal(:,2), dorsal(:,1), '-r', 'LineWidth', 1); hold on;
    
    % Plot skeleton
    [line, cv2, sp_curv_num] = spline_line(midline_mixed, spline_p, numcurvpts);
    plot(line(:,1), line(:,2), '-w', 'LineWidth', 2);
    
    % Store centerline information in all frames into centerline_data
    centerline_data(:, 2*(j-istart+1)-1) = line(:, 1); 
    centerline_data(:, 2*(j-istart+1)) = line(:, 2);

    % Interpolate to equally spaced length units
    cv2i = interp1(sp_curv_num+.00001*(0:length(sp_curv_num)-1), cv2, (0:(sp_curv_num(end)-1)/(numcurvpts+1):(sp_curv_num(end)-1)));
    df2 = diff(cv2i, 1, 1);
    atdf2 = unwrap(atan2(-df2(:,2), df2(:,1)));
    angle_data(j, :) = atdf2';
    
    % Collect curvature information for midline
    curv = unwrap(diff(atdf2, 1)); 
    curvdata((j-istart+1), :) = curv;
    curvdata(1,2:3)=curvdata(2,2:3);    
    curvdata(end,2:3)=curvdata(end-1,2:3);
  
end % End of calculations for 1 frame

hold off;

end % End of calculations for all frames

curvdata = curvdata'; % After transposition, columns are frame numbers, rows are segment numbers

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2); clf;

    imagesc(curvdata(:, :)*100); caxis([-10, 10]);
    title('Click on invalid columns, and then press return');
    badcols = ginput; % The first column stores the frame numbers, the second stores the segment numbers
  
    if ~isempty(badcols)
        bcols = [bcols'; ceil(badcols(:,1))]; % The first column stores the frame number, the second stores the segment number
        bcols = unique(bcols); % Delete repeated elements
    end
    
     for j = 1:length(bcols)
         
         col = bcols(j);
         disp(strcat('Invalid col #', num2str(j), '=', num2str(col)));
         
         % Replace first and last frame if it is invalid
         if col == 1
             curvdata(:, 1) = curvdata(:, 2);
         end

         if col == numframes
             curvdata(:, numframes) = curvdata(:, numframes-1);
         end
         
         % Replace invalid frame in between with its flanking frames
         if (col > 1) && (col < numframes)

             pre = col - 1;
             while ~isempty(find(bcols == pre, 1)) && (pre > 1)
                 pre = pre-1;
             end
             
             post = col + 1;
             while ~isempty(find(bcols == post, 1)) && (post < numframes)
                 post = post+1;
             end
             
             curvdata(:, col) = 0.5 * (curvdata(:, pre) + curvdata(:, post)); % Mean of previous and post frames
             
         end
         
     end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%Filter curvature by averaging%%%%%%%

h = fspecial('average', [2 10]);
curvdatafiltered = imfilter(curvdata*100,  h , 'replicate');
curvdatafiltered = smooth(curvdatafiltered);
curvdatafiltered = reshape(curvdatafiltered, numcurvpts, []);

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

end