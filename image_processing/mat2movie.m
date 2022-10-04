% function [] = mat2movie(filename, imagelist_g, imagelist_r, gfp, rfp, np, thresh, mini, maxi, r)
%% 
close all; clear;
%%
% setup_proof_reading;
% fprintf('tiff files loading finished. \n');
% imagelist_g = imagelist(1:2:end,1);
% imagelist_r = imagelist(2:2:end,1);
setup_proof_reading;
updated = 0;
imagelist_g = cellfun(@double, imagelist_g, 'uniformoutput', 0);
imagelist_r = cellfun(@double, imagelist_r, 'uniformoutput', 0);

fprintf('tiff loading completed. \n');

%% For processed movies, register channels with existing parameters for

% fpara = uigetfile('*.mat', ['Select a file for parameters of ' filename]);
% load(fpara);
% 
% figure;
% subplot(121);
% imtomove = imwarp(imagelist_r{1},tform,'OutputView',imref2d(size(imagelist_r{1})));
% imshowpair(imtomove, imagelist_g{1}); title('Left channel moved');
% subplot(122);
% imtomove = imwarp(imagelist_g{1},tform,'OutputView',imref2d(size(imagelist_g{1})));
% imshowpair(imagelist_r{1}, imtomove); title('Right channel moved');

channeltomove = questdlg('Which channel to move?', 'Channel', 'Left', 'Right', 'Right');
switch channeltomove
    case 'Left'
        [img_updated, movingRegistered, tform] = image_registration_tform_muscle...
            (imagelist_g, imagelist_r, channeltomove);
    case 'Right'            
        [img_updated, movingRegistered, tform] = image_registration_tform_muscle...
            (imagelist_r, imagelist_g, channeltomove);
end
% %% Register RFP and GFP channels
% figure; 
% subplot(1,2,1); imshowpair(imagelist_g{1,1}, imagelist_r{1,1});
% image_registration_tform; 
% subplot(1,2,2); imshowpair(imagelist_g{1,1}, imagelist_r{1,1});

%% Update the channel that has been moved

% Remove edges from registration
imagelist_moved = movingRegistered;
for i = 1:length(imagelist_g)
    img = movingRegistered{i};
    img = img + mean(img,[1 2])*double(img==0);
    imagelist_moved{i,1} = img;
end
fprintf('edges removed. \n');

% Update the channel
switch channeltomove
    case 'Left'
        imagelist_r = imagelist_moved;
    case 'Right'
        imagelist_g = imagelist_moved;
end
updated = 1;
fprintf('channel updated. \n');

% % Overlay two channels if necessary
% imagelist_use = imagelist_moved;
% for i = 1:length(imagelist)
%     imagelist_use{i,1} = imagelist_g{i,1}+imagelist_r{i,1};
% end
% fprintf('channels overlaid. \n');

%% Find proper adaptive threshold for the numerator image (normally for the GCaMP side)

if updated~=1
    fprintf('please update the registered channel first.\n');
else
    % Find proper adaptive threshold
    th = .3;
    % Apply threshold and filtering to all numerator images
    gsfilt = .6; 
    % gpsf = fspecial('gaussian',9,9);
    imagelist_gseg = imagelist_g;
    imagelist_segratio = imagelist_g;
    imagelist_segratiofilt = imagelist_g;
    % imagelist_segratiodecon = imagelist_g;

    for ni = 1:length(imagelist_g)
        imusedseg = imgaussfilt(imagelist_g{ni,1}+imagelist_r{ni,1},gsfilt);
        Itoseg = uint16(imusedseg);
        T = adaptthresh(Itoseg, th);
        segmask = double(imbinarize(Itoseg, T));
        if ni==1, close all; imshow(segmask), end
        imagelist_gseg{ni,1} = segmask.*imagelist_g{ni,1};
        imagelist_segratio{ni,1} = imagelist_gseg{ni,1}./imagelist_r{ni,1};
        imagelist_segratiofilt{ni,1} = imgaussfilt(imagelist_segratio{ni,1},gsfilt);
    %     imagelist_gsegdecon = deconvlucy(imagelist_gseg{ni,1},gpsf);
    %     imagelist_segratiofilt{ni,1} = imagelist_gseggfilt./imagelist_r{ni,1};
    %     imagelist_segratiodecon{ni,1} = imagelist_gsegdecon./imagelist_r{ni,1};
    end
end

%% Animation

close all;
bump_finder = bump_detector([1,1], [3,3], [51 51]);
normalization = 0.7 * sum_all(bump_finder(bump_finder>0)*255);
dc = @(x) deconvlucy(x, fspecial('gaussian', 5, 2), 1);
flt = @(x, h) imfilter(double(dc(x)), h, 'circular');
flt_and_normalize = @(x) uint8(flt(x, bump_finder)/normalization*256-5);

if ~(size(imagelist_g) == size(imagelist_r))  
    disp('Size mismatching');
else
    cropstart = 1; backstart = [49 550; 840 957]; relativebs = backstart-cropstart;
    framenum = length(imagelist_g);
    
    np_cropped = np_adj(cropstart:(cropstart+framenum-1),:);
    gfp_cropped = gfp(cropstart:(cropstart+framenum-1),:);
    rfp_cropped = rfp(cropstart:(cropstart+framenum-1),:);
    
    thresh = 3; 
    r = 8; mini = 0; maxi = 100;
    canvas_x = 3; canvas_y = size(np_cropped, 2)/2;
    cmp = lbmap(canvas_y, 'redblue');
    cmproi = colormap(inferno); cmproi(1,:) = [0 0 0];
    ratio = gfp_cropped./rfp_cropped; 
    ratio_norm = ratio./repmat(max(ratio), size(ratio,1) ,1);
    fps = 10; ff = 0.4;
    adj = 0.1;
    
    zeroS = zeros(r+1, r+1); % Add border to avoid bleeding of ROIs
    neuronnames = {'AVB', 'AVA/AVE'}; neuronnum = canvas_y;
    
    for i = 1:framenum
        
        %Denoise GFP image
%         imagelist_g_frame = imadjust(imagelist_g{i,1}, [thresh; 0.9], [], 0.1);
%         imagelist_g_logical = imagelist_g_frame > 0;
%         imagelist_g_logical = imagelist_g{i, 1} > thresh * mean2(imagelist_g{i, 1});
%         imagelogic = uint16(imagelist_g_logical) .* imagelist_g{i, 1};
%         flt_g = flt_and_normalize(imagelist_g{i,1});
%         flt_r = flt_and_normalize(imagelist_r{i,1});

%         flt_g = imagelist_g{i,1};
%         flt_r = imagelist_r{i,1};
%         imagelogic = flt_r > thresh*mean2(flt_r);
        
        %Create ratiometric image
%         imageratio = double(imagelist_g_denoised)./double(imagelist_r{i, 1});
%         imagedraw = double(imagelist_g{i, 1})./double(imagelist_r{i, 1});
%         imagedraw = blkdiag(zeroS, imageratio, zeroS);
%         imagedraw = blkdiag(zeroS, imagelist_g{i, 1}, zeroS);

%         imagedraw = double(imagelist_g{i,1}).*double(imagelogic)./double(imagelist_r{1,1});%         imagedraw = double(flt_g).*double(imagelogic)./double(flt_r);
        imagedraw = imagelist{i,1};        
        imagedraw = blkdiag(zeroS, imagedraw, zeroS);
        
        % Draw ratiometric image %%%%%%%%%%%%%%%%%%%%%%%%%
        figure(1);
%         set(gcf, 'color', 'w');
%         set(gcf,'position',[230 200 1000 400]);
        set(gcf,'units','normalized','outerposition',[0 0 1 1], 'color', 'k');
        
        hold off;
        % Example of plot
        subplot(canvas_y, canvas_x, 1:canvas_x:canvas_x*canvas_y);
        imagesc(imagedraw); colormap(cmproi);
        axis tight;
        set(gca, 'visible', 'off', 'ydir', 'reverse');
        
        c = colorbar('South');
        set(c, 'xcolor', 'w', 'ycolor', 'w');
        cpos = get(c, 'Position');
        cpos(4) = cpos(4) / 2;
        set(c, 'Position', cpos, 'FontSize', 8);        
        caxis([mini maxi]);
%         title(['Time ' num2str(i/fps, '%.2f') 's'], 'color', 'w');
%         hold on;
        text(10, 10, [num2str(i/fps, '%.2f') 's' ], 'Color', 'w');

        
        if i<relativebs(1) || i>relativebs(2)
            delete(findall(gcf,'type','annotation'));
            annotation('arrow', [0.16 0.14], [0.5 0.5], 'linewidth', 10, 'headwidth', 20, 'color', 'w');
        else
            delete(findall(gcf,'type','annotation'));
            annotation('arrow', [0.32 0.34], [0.5 0.5], 'linewidth', 10, 'headwidth', 20, 'color', 'w');
        end
%         zoom(1.4);
        % Add ROIs on image
        for j = 1:canvas_y
            
%             if np(i,2*j-1) > r && np(i,2*j) > r
                
                rectangle('Curvature', [0.5,0.5], 'Position', [np_cropped(i,2*j-1), np_cropped(i,2*j), 2*r, 2*r], 'EdgeColor', cmp(j,:), 'Linestyle', '-', 'Linewidth', 1);
                text(np_cropped(i,2*j-1)-2*r, np_cropped(i,2*j)+2*r, neuronnames(j), 'color', 'w', 'fontsize', 8);
            
%             else
                
%                 rmini = floor(min(np(i,2*j-1), np(i,2*j)));
%                 rectangle('Curvature', [0,0], 'Position', [np(i,2*j-1)-rmini, np(i,2*j)-rmini, 2*rmini, 2*rmini], 'EdgeColor', 'w', 'Linestyle', ':', 'Linewidth', 1);
%                 text(np(i,2*j-1)-rmini, np(i,2*j)-rmini, num2str(j), 'color', 'w');
%                 
%             end
        
        end
             

        hold off;
        
        % Draw zoomed in neuron images %%%%%%%%%%%%%%%%%%%%%%%%%
        for j = 1:canvas_y
            
            % Draw zoomed-in ROIs
            subplot(canvas_y, canvas_x, (canvas_y-j)*canvas_x+2);
            imagesc(imagedraw(round(np_cropped(i,2*j)):round(np_cropped(i,2*j))+2*r, round(np_cropped(i,2*j-1)):round(np_cropped(i,2*j-1))+2*r));
            caxis([mini maxi]); 
            axis image;
            set(gca, 'xticklabel', '', 'yticklabel', '', 'xtick', [], 'ytick', [], 'ydir', 'reverse');
            title(neuronnames(j), 'color', 'w');   
            
            ax = gca; ax.XColor = cmp(j,:); ax.YColor = cmp(j,:);
%             subplot(canvas_y, canvas_x, (j-1)*canvas_x+2);
%             b = bar3(imagedraw(round(np_cropped(i,2*j)):round(np_cropped(i,2*j))+2*r, round(np_cropped(i,2*j-1)):round(np_cropped(i,2*j-1))+2*r));
% %             for k = 1:length(b)
% %                 zdata = b(k).ZData;
% %                 b(k).CData = zdata;
% %                 b(k).FaceColor = 'interp';
% %             end
%             set(b, 'edgecolor', 'none'); 
% %             colormap(cmproi);
% %             axis normal; 
%             xlim([0 2*r]); ylim([0 2*r]); zlim([mini maxi]); 
%             set(gca, 'visible', 'off');
%             zoom(1.8);
            
            % Draw signals
            subplot(canvas_y, canvas_x, (canvas_y-j)*canvas_x+3);
            hold off;         
%             line([i i],[0 1],'color','w','linewidth',1);
            if i>relativebs(1) && i<relativebs(2)
                f = fill([relativebs(1) i i relativebs(1)], [0 0 1 1], 0.8*[1 1 1]);
                set(f, 'edgecolor', 'none', 'facealpha', 0.5);
            end; hold on; 
            plot(ratio_norm(1:i, j), 'color', cmp(j,:), 'linewidth', 2);
%             plot(gfp(1:i, j), 'color', cmp(j,:), 'linewidth', 2);
            ylim([0 1]); xlim([1 framenum]);
            clr = get(gcf, 'color');
            set(gca, 'xticklabel', '', 'yticklabel', '', 'xtick', [], 'color', clr);
            title(neuronnames(j), 'color', 'w'); 
            box off;
            
        end
        
        % Create gif file
        drawnow;
        frame = getframe(1);
        im = frame2im(frame);
        [imind, cm] = rgb2ind(im, 256);

        if i == 1
          imwrite(imind, cm, [filename '.gif'], 'gif', 'LoopCount', 1, 'DelayTime', ff/fps);
        else
          imwrite(imind, cm, [filename '.gif'], 'gif', 'WriteMode', 'append', 'DelayTime', ff/fps);
        end
        
%         % Store the frame
%         M(i) = getframe(gcf); % Leaving gcf out crops the frame in the movie.
        
    end

 end