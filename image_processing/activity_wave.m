function [intensity_total_fil, M] = activity_wave(imagelist_ratio, edge_data, centerline_data_spline, curvdatafiltered, thresh_size, one_vertices_zero_circle, one_avi_zero_gif, cutoff_or_not)

    % Prepare data    
    col_num = size(imagelist_ratio{1, 1}, 1);
    row_num = size(imagelist_ratio{1, 1}, 2);
    
    xq = 1:row_num;
    xq_rep = repmat(xq, col_num, 1);
    xq_rep_shape = reshape(xq_rep, [], 1);
    
    yq = 1:col_num;
    yq_rep = repmat(yq, row_num, 1);
    yq_rep_shape = reshape(yq_rep', [], 1);
    
    length_all = size(centerline_data_spline, 1); % Number of points for a worm: 100 points and 99 segments as default
    length_seg = length_all-1;
    frame_num = size(imagelist_ratio, 1);
    intensity_total = zeros(length_seg, frame_num);
    intensity_frame = zeros(length_seg, 1);
    center_seg = zeros(length_all-1, 2);
 	masksize = zeros(length_all-1, 1);
    color_code = zeros(length_all-1, 3);

    M(1:frame_num) = struct('cdata', [], 'colormap', []);
    fps = 20;
    
    segment_intensity_cutoff = 200; % Based on distribution of intensity for one segment
    j = 0; % Counting the oversized intensity_frame
    upbnd = 1.1;
    lowbnd = 0.9;
   	gray = 0.1; % Edge grayness for bubble plot
    mag_thresh = 5; % Magnification for the size of mask

%%%%%%%%%%%%%%%%%%%%%%%%%% Calculation of intensity in all segments in all frames %%%%%%%%%%%%%%%%%%%%%%%%%

    % Loop for every frame
    for k = 1:frame_num
    
%         hold off;
%         figure(1);
%         set(gcf, 'Color', [1 1 1]);
%         disp(k);
        frame = imagelist_ratio{k, 1};
        frame(isnan(frame)) = 0; % Replace all NaNs with zero
        
        path = [edge_data{2*k-1, 1} edge_data{2*k, 1}];
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
   
    % User choose the bubble representation
%             if ~one_vertices_zero_circle
% 
%                 hold off;
%                 imagesc(imagelist_ratio{k, 1});                
%                 title('Ratiometric Calcium Recording');
%                 set(gca, 'XTick', [], 'YTick', []);
%                 axis image; hold on; box off;
%                 
%                 text(10, 10, ['Frame # ' num2str(k)], 'Color', 'w');
%                 
%                 % Obtain points in all segments for one worm
%                 for i = 1:length_all-1 % Number of segments is 1 less than number of points
% 
%                     % Shape the segment
%                     xv = [path(i, 1) path(i+1, 1) centerline_data_spline(i+1, 2*k-1) centerline_data_spline(i, 2*k-1) path(i, 1)]; % Row number for polygon joints
%                     yv = [path(i, 2) path(i+1, 2) centerline_data_spline(i+1, 2*k) centerline_data_spline(i, 2*k) path(i, 2)]; % Col number for polygon joints
% 
%                     % Points inside one segment
%                     in = inpolygon(xq_rep_shape, yq_rep_shape, xv, yv); % The logical value for points inside the segment
%                     y_ind = yq_rep_shape(in); % Col number for points inside segment
%                     x_ind = xq_rep_shape(in); % Row number for points inside segment
% 
% %                     empty_image = frame*0;
% %                     for num = 1: length(y_ind)
% %                         empty_image(y_ind(num), x_ind(num)) = frame(y_ind(num), x_ind(num));
% %                     end
% %                     
% %                     intensity_frame(i) = nansum(nansum(empty_image));
% 
%                     % Calculate total intensity value for one segment
%                     intensity_frame(i) = sum(diag(frame(y_ind, x_ind))); % Diagonal of the submatrix contains all the points within polygon                      
%                     
%                     if isnan(intensity_frame(i)) % When no points are inside the segment
%                         intensity_frame(i) = 0;
%                     elseif intensity_frame(i) > segment_intensity_cutoff
%                         intensity_frame(i) = segment_intensity_cutoff; % Restrict the intensity value to the a posteriori cutoff
%                         j = j + 1;
%                     end
% 
%                     % Plot the circle centered in the middle of the segment
%                     center_seg(i, :) = [(path(i, 1)+path(i+1, 1))/2, (path(i, 2)+path(i+1, 2))/2];
%                    	masksize(i) = intensity_frame(i) / thresh_size + 1; % Prevent it from reaching zero
%                                         
%                     if (i == 1) || (intensity_frame(i) > upbnd*intensity_frame(i-1))
%                         color_code(i, :) = [1 0 0]; % Red
%                     elseif (intensity_frame(i) <= upbnd*intensity_frame(i-1)) && (intensity_frame(i) >= lowbnd*intensity_frame(i-1))
%                         color_code(i, :) = [1 1 0]; % Yellow
%                     else
%                         color_code(i, :) = [0 1 0]; % Green
%                     end
% 
% %                     plot(center_seg(i, 1), center_seg(i, 2), 'o', 'MarkerSize', masksize(i), 'MarkerFaceColor', color_code(i, :), 'MarkerEdgeColor', [1 1 1]*gray);
%                     
%                 end % End of all segments for one frame in bubble representation
%                 
%                 scatter(center_seg(:, 1), center_seg(:, 2), masksize.^2, color_code, 'filled', 'MarkerEdgeColor', [1 1 1]*gray);
%                 filename = 'Muscle_Activity_Bubble';
            
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%             else % User chose the wave-spectrum representation
                
                subplot(1, 2, 1);
                hold off;
                imagesc(imagelist_ratio{k, 1});
                axis image; hold on; box off;
                set(gca, 'XTick', [], 'YTick', []); xlabel('Ratiometric Calcium Recording');
                
                text(10, 10, ['Frame # ' num2str(k)], 'Color', 'w');
                
                % Obtain all points in segments for one worm
                for i = 1:length_seg % Number of segments is 1 less than number of points

                    % Shape the segment
                    xv = [path(i, 1) path(i+1, 1) centerline_data_spline(i+1, 2*k-1) centerline_data_spline(i, 2*k-1) path(i, 1)]; % Row number for polygon joints
                    yv = [path(i, 2) path(i+1, 2) centerline_data_spline(i+1, 2*k) centerline_data_spline(i, 2*k) path(i, 2)]; % Col number for polygon joints

                    % Points inside one segment
                    in = inpolygon(xq_rep_shape, yq_rep_shape, xv, yv); % The logical value for points inside the segment
                    y_ind = yq_rep_shape(in); % Col number for points inside segment
                    x_ind = xq_rep_shape(in); % Row number for points inside segment

                    % Obtain total intensity value for one segment
                    intensity_frame(i) = sum(diag(frame(y_ind, x_ind))); % Diagonal of the submatrix contains all the points within polygon
                                       
%                     if (cutoff_or_not * intensity_frame(i)) > segment_intensity_cutoff
%                         intensity_frame(i) = segment_intensity_cutoff;
%                         j = j + 1;
%                     end

                    
                end % End of all segments for one frame
                
    %                     masksize = ceil(intensity_frame(i)/thresh_size);
    %                     dirvec = path(i+1, 1:2)-path(i, 1:2);
    %                       dirvec = dirvec/norm(dirvec);
    %                     perpvec = [dirvec(2), -dirvec(1)];
    %                       perpvec = perpvec/norm(perpvec);
    %                     ver1 = path(i, 1:2) + masksize*(perpvec);
    %                     ver2 = path(i+1, 1:2) + masksize*(perpvec);
    %                    p = patch([path(i, 1) ver1(1) ver2(1) path(i+1, 1) path(i, 1)], [path(i, 2) ver1(2) ver2(2) path(i+1, 2) path(i, 2)], [1 0.5 0]);
    %                    set(p, 'FaceAlpha', 0.6);
                
%                 subplot(1,2,2);                 
%                 hold off;
%                 bar(intensity_frame/thresh_size, 'FaceColor', [1 0.6 0], 'EdgeColor', 'none');
%                 ylim([0 mag_thresh * thresh_size]);
%                 hold on;
%                 imagesc(curvdatafiltered(1:end-1, k)');
%                 set(gca, 'XTick', []); box off;
%                 xlabel('Activity and Curvature for Each Segment');
                
%                 filename = 'Muscle_Activity_Wave';

%             end % End of if-else for bubble or wave-spectrum format
                      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%         disp(num2str(sum(intensity_frame)));
        % Store total intensity for one frame
        intensity_total(:, k) = intensity_frame;
        
%%%%%%%%%%%%%%%%%%%%%%%%%% Animation for all frames %%%%%%%%%%%%%%%%%%%%%%%%%

%         % Choose format to save, gif or avi    
%         if ~one_avi_zero_gif % Choose gif
% 
%             % Store the frame for gif
%             drawnow;
%             frame = getframe(1);
%             im = frame2im(frame);
%             [imind, cm] = rgb2ind(im, 256);
% 
%             if k == 1
%               imwrite(imind, cm, [filename '.gif'], 'gif', 'LoopCount', Inf, 'DelayTime', 1/fps);
%             else
%               imwrite(imind, cm, [filename '.gif'], 'gif', 'WriteMode', 'append', 'DelayTime', 1/fps);
%             end
% 
%         else % Choose avi
%             
%             M(k) = getframe(gcf); % Leaving gcf out crops the frame in the movie.
% 
%         end % End of storing data into animation
    
    end % End of all frames
    
%     disp(['Total number of oversized masks is ' num2str(j) '.']);
    timefilter = 2; bodyfilter = 10; 
    h = fspecial('average', [timefilter bodyfilter]);
    numframes = size(intensity_total, 2);
    intensity_total_fil = imfilter(intensity_total,  h , 'replicate');
    intensity_total_fil = smooth(intensity_total_fil);
    intensity_total_fil = reshape(intensity_total_fil, [], numframes);
    
    beep;

end % End of function
