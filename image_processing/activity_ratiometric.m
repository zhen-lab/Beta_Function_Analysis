%% Clear everything

clear; close all;

%% Load tiff files

setup_proof_reading;
updated = 0;
imagelist_g = cellfun(@double, imagelist_g, 'uniformoutput', 0);
imagelist_r = cellfun(@double, imagelist_r, 'uniformoutput', 0);

fprintf('tiff loading completed. \n');

%% Registration of GFP and RFP channels

if ~(filename==0)
    
    fprintf('tiff files loading finished. \n');
    channeltomove = questdlg('Move left or right channel?',...
        'Channel to move', 'Left', 'Right', ...
        'Right');
    
    % Register RFP and GFP channels
    switch channeltomove
        case 'Left'
            figure;
            subplot(1,2,1); imshowpair(imagelist_r{1}, imagelist_g{1});
            [img_updated, movingRegistered, tform] = image_registration_tform_muscle...
                (imagelist_g, imagelist_r, channeltomove);
            subplot(1,2,2); imshowpair(movingRegistered{1}, imagelist_g{1});
        case 'Right'            
            figure;
            subplot(1,2,1); imshowpair(imagelist_r{1}, imagelist_g{1});
            [img_updated, movingRegistered, tform] = image_registration_tform_muscle...
                (imagelist_r, imagelist_g, channeltomove);
            subplot(1,2,2); imshowpair(imagelist_r{1}, movingRegistered{1});
    end
    
end

fprintf(['registration for the first frame completed for \n' filename '\n']);

% %% (Master Data Struct) Check using analyzed data 
% % For those processed before calculated registration
% % This part is for master data file structures
% masdata = uigetfile('.mat', 'Select master data file');
% load(masdata);
% recnumtotal = size(mst_curvdatafiltered,1);
% prompt = 'DV data from recording #';
% inputdlgtitle = 'Recording #';
% dims = [1 35];
% definput = {'1'};
% recnum = inputdlg(prompt,inputdlgtitle,dims,definput);
% recnumnum = str2double(recnum);
% centerline_data_spline = mst_centerline_data_spline{recnumnum};
% curvdatafiltered = mst_curvdatafiltered{recnumnum};
% dorsal_data = mst_dorsal_data{recnumnum};
% ventral_data = mst_ventral_data{recnumnum};
% hold on;
% plot(dorsal_data{1,1}, dorsal_data{2,1}, 'r');
% plot(ventral_data{1,1}, ventral_data{2,1}, 'b');
% fprintf('contour overlaid. \n');

%% Update the channel that has been moved

% Remove edges from registration
imagelist_moved = movingRegistered;
for i = 1:(range(2)-range(1)+1)
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
    th = .42;
    % Apply threshold and filtering to all numerator images
    gsfilt = .6; 
    % gpsf = fspecial('gaussian',9,9);
    imagelist_gseg = imagelist_g;
    imagelist_segratio = imagelist_g;
    imagelist_segratiofilt = imagelist_g;
    % imagelist_segratiodecon = imagelist_g;

    for ni = 1:(range(2)-range(1)+1)
%         % Top/Bottom hat
%         img = double(imagelist_g{ni,1});
%         se = strel('disk',100);
%         imgsub = ...
%             imsubtract(img,imbothat(img,se));
%         imagelist_gseg{ni,1} = imgsub;
%         if ni==1, close all; subplot(121); imagesc(img); subplot(122); imagesc(imgsub), end
%         imagelist_segratio{ni,1} = imagelist_gseg{ni,1}./imagelist_r{ni,1};
%         imagelist_segratiofilt{ni,1} = imgaussfilt(imagelist_segratio{ni,1},gsfilt);
        % Adaptive thresholding
        imusedseg = imgaussfilt(imagelist_g{ni,1}+imagelist_r{ni,1},gsfilt);
        Itoseg = uint16(imusedseg);
        T = adaptthresh(Itoseg, th);
        segmask = double(imbinarize(Itoseg, T));
        if ni==1, close all; imshow(segmask), end
        imagelist_gseg{ni,1} = segmask.*double(imagelist_g{ni,1});
        imagelist_segratio{ni,1} = imagelist_gseg{ni,1}./imagelist_r{ni,1};
        imagelist_segratiofilt{ni,1} = imgaussfilt(imagelist_segratio{ni,1},gsfilt);
%         % Deconvolution
%         imagelist_gsegdecon = deconvlucy(imagelist_gseg{ni,1},gpsf);
%         imagelist_segratiofilt{ni,1} = imagelist_gseggfilt./imagelist_r{ni,1};
%         imagelist_segratiodecon{ni,1} = imagelist_gsegdecon./imagelist_r{ni,1};
    end
end

%% Test segmented images and filtering
close all;
clrtheme = 'inferno';
mini = 0; 
maxi = 3; 
subplot(131); imagesc(imagelist_g{end,1}); caxis([0 3000]); title('GCaMP'); axis equal; colormap(clrtheme);
subplot(132); imagesc(imagelist_segratio{end,1}); caxis([mini maxi]); title('Segmented ratio'); axis equal; colormap(clrtheme);
subplot(133); imagesc(imagelist_segratiofilt{end,1}); caxis([mini maxi]); title('Segmented filtered ratio'); axis equal; colormap(clrtheme);

%%
close all;

if ~exist('maxi','var')
    fprintf('please test segmented images first.\n');
% elseif ~exist('curvdatafiltered','var')
%         fprintf('please load analyzed data first. \n');
else
    subplot(121); imagesc(imagelist_segratiofilt{1,1});
    subplot(122); imagesc(imagelist_segratiofilt{end,1});
    prompt = {'Total time length (second):','Total frame number:', ...
        'Play speed adjustment:',...
        'Rotation (clockwise +, counterclockwise -):',...
        'Flip? (up/down 1, left/right -1, no 0):'};
    inputdlgtitle = 'Parameters for GIF';
    dims = [1 35];
    definput = {'68',num2str(size(imagelist_g,1)),'0.3','0','0'}; % for muscle activity on Fig. 1: '438','10','2','east'
    answer = inputdlg(prompt,inputdlgtitle,dims,definput);
    ftl = str2double(answer{1,1});
    frmnum = str2double(answer{2,1});
    fps = frmnum/ftl;
    ff = str2double(answer{3,1});
    rotangle = str2double(answer{4,1}); % clockwise: positive value, counterclosewise: negative value
    flipside = answer{5,1};
    % r = 10;
    % backperiod = [109 193];
    % neuronnames = {'AVB', 'AVA/AVE'}; neuronnum = size(neuronnames, 2);
    % mini = 0; 
    % maxi = 1; 
    % ftl = 90;
    % fps = 1800/ftl; 
    % ff = 0.3;
    lightstart = 1; 
    lightend = range(2)-range(1)+1;
    cmp = colormap(inferno); % cmp(1,:) = [0 0 0];
    % cutoff = mean2(imagelist{1,1});
    % rotangle = 0; % clockwise: positive value, counterclosewise: negative value
    close all;
    
    for frm = lightstart:lightend

        % Plotting muscles
        
        figure(1);
%         subplot(121); % Actual recording
    %     set(gcf,'units','normalized','outerposition',[0 0 1 1], 'color', 'k');
        set(gcf, 'color', 'w', 'pos', [10 10 500 600]);
        hold off;
    %     imagedraw = imagelist{frm,1};
    %     imagedraw = imrotate(rot90(imagelist_segratiofilt{frm,1},2), -rotangle, 'bilinear', 'crop');
        imagedraw = imrotate(imagelist_segratiofilt{frm,1}, -rotangle, 'bilinear', 'crop');
        switch flipside
            case '1'
                imagedraw = flipud(imagedraw);
            case '-1'
                imagedraw = fliplr(imagedraw);
        end
    %     imagelogic = imagedraw>cutoff;
    %     imagedraw = double(imagedraw).*double(imagelogic);
        imagesc(imagedraw); % Have not added borders!!!!
        hold on;
        text(10, 10, [num2str((frm-lightstart)/(frmnum/ftl), '%.2f') 's' ], 'Color', 'w', 'fontsize', 10);

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
        
%         imagex = size(imagedraw,2); imagey = size(imagedraw,1);
%         ll = 0.1*min(imagex, imagey);
%         line([0.9*imagex 0.9*imagex], [0.8*imagey 0.8*imagey+ll], 'color','w', 'linewidth', 1);
%         line([0.9*imagex-ll/2 0.9*imagex+ll/2], [0.8*imagey+ll/2 0.8*imagey+ll/2], 'color','w', 'linewidth', 1);        
%         text(0.9*imagex, 0.8*imagey+ll, 'V', 'Color', 'w', 'fontsize', 10, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
%         text(0.9*imagex, 0.8*imagey, 'D', 'Color', 'w', 'fontsize', 10, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
%         text(0.9*imagex-ll/2, 0.8*imagey+ll/2, 'A', 'Color', 'w', 'fontsize', 10, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
%         text(0.9*imagex+ll/2, 0.8*imagey+ll/2, 'P', 'Color', 'w', 'fontsize', 10, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');

        % Ploting neurons
        hold on;
        for nn = 1:neuronnum
            plot(dual_position_data{frm,nn}(1)-n,dual_position_data{frm,nn}(2),'ow', 'markersize', 4*r);
            text(dual_position_data{frm,nn}(1)-n-2*r, dual_position_data{frm,nn}(2)+2*r, neuronnames(nn), 'color', 'w', 'fontsize', 10);
        end

        text(10, 40, [num2str((frm-lightstart)/fps, '%.2f') 's' ], 'Color', 'w', 'fontsize', 10);

        if frm<backperiod(1) || frm>backperiod(2)
            delete(findall(gcf,'type','annotation'));
            annotation('arrow', [0.41 0.39], [0.5 0.5], 'linewidth', 10, 'headwidth', 20, 'color', 'w');
        else
            delete(findall(gcf,'type','annotation'));
            annotation('arrow', [0.61 0.63], [0.5 0.5], 'linewidth', 10, 'headwidth', 20, 'color', 'w');
        end
            
%     subplot(122); % Segmentation to visualize orientations
%     hold on;
%     plot(dorsal_data{2*frm-1,1}, dorsal_data{2*frm,1}, ':r');
%     plot(ventral_data{2*frm-1,1}, ventral_data{2*frm,1}, ':b');
%     plot(centerline_data_spline(:,2*frm-1), centerline_data_spline(:,2*frm), ':k');
%     plot(centerline_data_spline(1,2*frm-1), centerline_data_spline(1,2*frm), ':og');
%     plot(centerline_data_spline(end,2*frm-1), centerline_data_spline(end,2*frm), ':oy');
%     set(gca, 'ydir', 'reverse');
%     axis equal;
    
        % Create gif file
        drawnow;
        frame = getframe(1);
        im = frame2im(frame);
        [imind, cm] = rgb2ind(im, 256);
        if frm == lightstart
          imwrite(imind, cm, [filename '_regsegfiltrat.gif'], 'gif', 'LoopCount', inf, 'DelayTime', ff/fps);
        else
          imwrite(imind, cm, [filename '_regsegfiltrat.gif'], 'gif', 'WriteMode', 'append', 'DelayTime', ff/fps);
        end

    end

    save([filename '_movpara.mat'], 'tform', 'th', 'gsfilt', 'rotangle', 'ftl', 'ff', 'mini', 'maxi');

end

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