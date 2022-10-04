%% Clear everything

clear; close all;

%% Load tiff files

updated = 0;
nummov = 4;
imagelist_all = cell(nummov,1);
imagelist_g_all = cell(nummov,1);
imagelist_r_all = cell(nummov,1);
imagelist_reg_all = cell(nummov,1);
range_all = cell(nummov,1);
filename_all = cell(nummov,1);

for nm = 1:nummov

    setup_proof_reading;
    imagelist_all{nm} = imagelist;
    imagelist_g_all{nm} = cellfun(@double, imagelist_g, 'uniformoutput', 0);
    imagelist_r_all{nm} = cellfun(@double, imagelist_r, 'uniformoutput', 0);
    filename_all{nm} = filename;
    range_all{nm} = range;
    
    % Clear other variables to enable loading for the next recording
    clearvars -except imagelist_all imagelist_g_all imagelist_r_all filename_all range_all nummov

end

fprintf('tiff loading completed. \n');

%% Registration of GFP and RFP channels

for nm = 1:nummov

    if ~(filename_all{nm,1}==0)
        
        imagelist_g = imagelist_g_all{nm,1};
        imagelist_r = imagelist_r_all{nm,1};
        filename = filename_all{nm,1};
%         filemovpara = uigetfile('.mat','Select movpara file');
%         load(filemovpara);
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
                figure;        
            case 'Right'            
                subplot(1,2,1); imshowpair(imagelist_r{1}, imagelist_g{1});
                [img_updated, movingRegistered, tform] = image_registration_tform_muscle...
                    (imagelist_r, imagelist_g, channeltomove);
                subplot(1,2,2); imshowpair(imagelist_r{1}, movingRegistered{1});
        end
        
    end
    
    imagelist_reg_all{nm,1} = movingRegistered;

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

    % Update the channel that has been moved
    
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

    imagelist_g_all{nm,1} = imagelist_g;
    imagelist_r_all{nm,1} = imagelist_r;
    fprintf('channel updated. \n');
    
    % % Overlay two channels if necessary
    % imagelist_use = imagelist_moved;
    % for i = 1:length(imagelist)
    %     imagelist_use{i,1} = imagelist_g{i,1}+imagelist_r{i,1};
    % end
    % fprintf('channels overlaid. \n');

end

fprintf('all channels updated. \n');

%% Find proper adaptive threshold for the numerator image (normally for the GCaMP side)

th_all = cell(nummov,1);
gsfilt_all = cell(nummov,1);
imagelist_segratiofilt_all = cell(nummov,1);

for nm = 1:nummov

    imagelist_g = imagelist_g_all{nm,1};
    imagelist_r = imagelist_r_all{nm,1};

    if updated~=1
        fprintf('please update the registered channel first.\n');
    else
        prompt = {'Enter threshold:','Enter Gaussian filter:'};
        dlgtitle = 'Input';
        dims = [1 35];
        definput = {'.34','.6'};
        answer = inputdlg(prompt,dlgtitle,dims,definput);
        th = str2double(answer{1,1});
        gsfilt = str2double(answer{2,1});
        
%         % Find proper adaptive threshold
%         th = .34;
%         % Apply threshold and filtering to all numerator images
%         gsfilt = .6; 
        % gpsf = fspecial('gaussian',9,9);
        imagelist_gseg = imagelist_g;
        imagelist_segratio = imagelist_g;
        imagelist_segratiofilt = imagelist_g;
        % imagelist_segratiodecon = imagelist_g;
    
        for ni = 1:length(imagelist_g)
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
            imagelist_gseg{ni,1} = segmask.*double(imgaussfilt(imagelist_g{ni,1},gsfilt));
            imagelist_segratio{ni,1} = imagelist_gseg{ni,1}./imagelist_r{ni,1};
            imagelist_segratiofilt{ni,1} = imgaussfilt(imagelist_segratio{ni,1},gsfilt);
    %         % Deconvolution
    %         imagelist_gsegdecon = deconvlucy(imagelist_gseg{ni,1},gpsf);
    %         imagelist_segratiofilt{ni,1} = imagelist_gseggfilt./imagelist_r{ni,1};
    %         imagelist_segratiodecon{ni,1} = imagelist_gsegdecon./imagelist_r{ni,1};
        end
    end

    th_all{nm,1} = th;
    gsfilt_all{nm,1} = gsfilt;
    imagelist_segratiofilt_all{nm,1} = imagelist_segratiofilt;

end
    
fprintf('all images thresholded and filtered. \n');

%% Test segmented images and filtering

close all;
maxi_all = cell(nummov);
mini_all = cell(nummov);

for nm = 1:nummov
    
    imagelist_g = imagelist_g_all{nm,1};
    imagelist_segratiofilt = imagelist_segratiofilt_all{nm,1};

    prompt = {'Enter maximum:','Enter minimum:'};
    dlgtitle = 'Input';
    dims = [1 35];
    definput = {'1','0'};
    answer = inputdlg(prompt,dlgtitle,dims,definput);
    maxi = str2double(answer{1,1});
    mini = str2double(answer{2,1});
    
    figure(nm)
    clrtheme = 'inferno';
    subplot(121); imagesc(imagelist_g{end,1}); caxis([0 3000]); title('GCaMP'); axis equal; colormap(clrtheme);
%     subplot(132); imagesc(imagelist_segratio{end,1}); caxis([mini maxi]); title('Segmented ratio'); axis equal; colormap(clrtheme);
    subplot(122); imagesc(imagelist_segratiofilt{end,1}); caxis([mini maxi]); title('Segmented filtered ratio'); axis equal; colormap(clrtheme);
    
    mini_all{nm,1} = mini; 
    maxi_all{nm,1} = maxi; 

end

fprintf('all lookup table set. \n');

%% Make movies

close all;
% imagelist_g = imagelist_g_all{1,1};
% lightstart = 1; 
% lightend = length(imagelist_g);

for nm = 1:nummov
    
    imagelist_g = imagelist_g_all{nm,1};
    imagelist_segratiofilt = imagelist_segratiofilt_all{nm,1};
    filename = filename_all{nm,1};
%     if nm>1
%         lightstart = lightstart + length(imagelist_g_all{nm-1,1});
%         lightend = lightend + length(imagelist_g);
%     end
    
    if ~exist('maxi','var')
        fprintf('please test segmented images first.\n');
    else
        [fileanalyzeddata, fileanalyzedpath] = uigetfile('.mat', ['Select analyzed data for ' filename]);
        load(fullfile(fileanalyzedpath, fileanalyzeddata));
        subplot(121); imagesc(imagelist_segratiofilt{1,1});
        subplot(122); imagesc(imagelist_segratiofilt{end,1});
        prompt = {'Total time length (second):',...
            'Total frame number:', ...
            'Play speed adjustment:'};
    %         'Flip? (up/down 1, left/right -1, no 0):'};
    %         'Rotation (clockwise +, counterclockwise -):',...
        inputdlgtitle = 'Parameters for GIF';
        dims = [1 35];
        definput = {'180',num2str(size(imagelist_g,1)),'0.3','0'}; % for muscle activity on Fig. 1: '438','10','2','east'
        answer = inputdlg(prompt,inputdlgtitle,dims,definput);
        ftl = str2double(answer{1,1});
        frmnum = str2double(answer{2,1});
        fps = frmnum/ftl;
        ff = str2double(answer{3,1});
    %     rotangle = str2double(answer{4,1}); % clockwise: positive value, counterclosewise: negative value
    %     flipside = answer{4,1};
        % r = 10;
        % backperiod = [109 193];
        % neuronnames = {'AVB', 'AVA/AVE'}; neuronnum = size(neuronnames, 2);
        % mini = 0; 
        % maxi = 1; 
        % ftl = 90;
        % fps = 1800/ftl; 
        % ff = 0.3;
        lightstart = 1; 
        lightend = length(imagelist_g);
        cmp = colormap(inferno); % cmp(1,:) = [0 0 0];
        % cutoff = mean2(imagelist{1,1});
        % rotangle = 0; % clockwise: positive value, counterclosewise: negative value
%         close all;
        
        % Determine rotation angle
        nump = 40;
        segsize = size(centerline_data_spline,1);
        antP = segsize*0.4; 
        posP = segsize*0.6;
        posPseries = uint16(linspace(posP,segsize,nump));
        antPseries = uint16(linspace(1,antP,nump));
        rotangleall = zeros(frmnum,1);
        for numangle = 1:frmnum
            centerline_frm = [centerline_data_spline(:,2*numangle-1) centerline_data_spline(:,2*numangle)];
            xside = centerline_frm(posPseries,1) - centerline_frm(antPseries,1);
            yside = centerline_frm(posPseries,2) - centerline_frm(antPseries,2);
            rotangleall(numangle) = mean(-atan2(yside, xside)); % radian
        end
        rotunwraprad = rad2deg(smoothdata(unwrap(rotangleall),1,'movmedian',uint16(frmnum/10)));
    
        % Stablize the centroids
        numpolygon = 6;
        polygonseries = uint16(linspace(1,segsize,numpolygon));
        centAll = zeros(uint16(segsize/2),2);
        for pi = 1:size(centerline_data_spline,2)/2
            polyx = centerline_data_spline(polygonseries,2*pi-1);
            polyy = centerline_data_spline(polygonseries,2*pi);
            polyin = polyshape(polyx,polyy); warning off;
            [centX,centY] = centroid(polyin);
            centAll(pi,:) = [centX centY];
        end
        centSmoothed = smoothdata(centAll,...
            1,'movmedian',uint16(frmnum/5));
    %     figure; plot(centAll);
    %     hold on; plot(centSmoothed)
    
        % Caculate body length to determine boundaries
        centdiffx = diff(centerline_data_spline(:,1:2:end),1,1);
        centdiffy = diff(centerline_data_spline(:,2:2:end),1,1);
        bodyhalflength = double(max(sum(sqrt(centdiffx.^2 + centdiffy.^2)))/2);
    %     bbsize = 2;
        
        padded = 0; padsize = 0; imagelist_segratiofiltuse = imagelist_segratiofilt;
        % Pad array if bounding box exceeds boundary
        if any(centSmoothed<bodyhalflength,'all')
            bboxdiff = ceil(bodyhalflength-centSmoothed);
            padsize = max(bboxdiff,[],'all');
            for pdi = 1:length(imagelist_segratiofilt)
                imagelist_segratiofiltuse{pdi,1} = padarray(imagelist_segratiofilt{pdi,1},...
                    [padsize padsize],0,'both');
            end
            padded = 1;
            fprintf('image padded.\n');
        else
            imagelist_segratiofiltuse = imagelist_segratiofilt;
        end
        
        movname = inputdlg(['Move name for ' filename]);
        for frm = lightstart:lightend
    
            % Plotting muscles
            
            figure(1);
            set(gcf,'units','normalized','outerposition',[0 0 0.5 1]);
            subplot(3,2,[1 3]); % Actual recording
    %         set(gcf, 'color', 'w', 'pos', [10 10 500 600]);
            hold off;
            imagedraw = imagelist_segratiofiltuse{frm,1}; 
            imagex = size(imagedraw,2); imagey = size(imagedraw,1);
            centerpoint_frm = centSmoothed(frm,:) + padded*padsize;
            centerline_frm = [centerline_data_spline(:,2*frm-1) centerline_data_spline(:,2*frm)] + padded*padsize;
            dorsal_data_frm = [dorsal_data{2*frm-1,1} dorsal_data{2*frm,1}] + padded*padsize;
            ventral_data_frm = [ventral_data{2*frm-1,1} ventral_data{2*frm,1}] + padded*padsize;
        %     imagedraw = imrotate(rot90(imagelist_segratiofilt{frm,1},2), -rotangle, 'bilinear', 'crop');
    %         imagedraw = imrotate(imagelist_segratiofilt{frm,1}, -rotangle, 'bilinear', 'crop');
    %         switch flipside
    %             case '1'
    %                 imagedraw = flipud(imagedraw);
    %                 centerline_frm(:,2) = -centerline_data_spline(:,2);
    %                 dorsal_data_frm(:,2) = -dorsal_data_frm(:,2);
    %                 ventral_data_frm(:,2) = -ventral_data_frm(:,2);
    %             case '-1'
    %                 imagedraw = fliplr(imagedraw);
    %                 centerline_frm(:,1) = -centerline_data_spline(:,1);
    %                 dorsal_data_frm(:,1) = -dorsal_data_frm(:,1);
    %                 ventral_data_frm(:,1) = -ventral_data_frm(:,1);
    %         end
    %         % Determine rotation angle
    %         xside = centerline_frm(midp,1) - ventral_data_frm(midp,1);
    %         yside = centerline_frm(midp,2) - ventral_data_frm(midp,2);
    % %         tanrotangle = xside/yside;
    %         rotangle = atan2d(xside, yside);
            % Rotate counterclockwise
            rotangle = rotunwraprad(frm);
            imagedrawrot = imrotate(imagedraw, -rotangle, 'bilinear');
            ctrOrig = [size(imagedraw,2);size(imagedraw,1)]/2;
            ctrRot = [size(imagedrawrot,2);size(imagedrawrot,1)]/2;        
    %         Rp = [cosd(-rotangle) -sind(-rotangle); sind(-rotangle) cosd(-rotangle)];
            R = [cosd(rotangle) -sind(rotangle); sind(rotangle) cosd(rotangle)];
    %         ctrWorm = [centerline_frm(uint16(segsize/2),1); centerline_frm(uint16(segsize/2),2)];
    %         ctrWormrot = (Rp*(ctrWorm-ctrOrig)+ctrRot)';
            vcPS = centerpoint_frm'; vcPSrot = (R*(vcPS-ctrOrig)+ctrRot)';
            vc = centerline_frm'; vcrot = (R*(vc-ctrOrig)+ctrRot)';
            vd = dorsal_data_frm'; vdrot = (R*(vd-ctrOrig)+ctrRot)';
            vv = ventral_data_frm'; vvrot = (R*(vv-ctrOrig)+ctrRot)';
            
            % Determine if flipping upside down is necessary
            segcent = uint16(segsize/2);
            cCenter = vcrot(segcent,:);
            vCenter = vvrot(segcent,:);
            dCenter = vdrot(segcent,:);
            vOrient = dot(vCenter-cCenter,[0;1]);
            dOrient = dot(dCenter-cCenter,[0;1]);
    %         if vvrot(posP,2)<vdrot(posP,2) || vvrot(antP,2)<vdrot(antP,2) % y direction upside down, therefore ventral is on top but then shown at bottom
            if vOrient<0 || dOrient>0
                imagedrawrot = flipud(imagedrawrot);
                vcPSrot(:,2) = 2*ctrRot(2)-vcPSrot(:,2);            
                vcrot(:,2) = 2*ctrRot(2)-vcrot(:,2);
                vdrot(:,2) = 2*ctrRot(2)-vdrot(:,2);
                vvrot(:,2) = 2*ctrRot(2)-vvrot(:,2); 
    %             fprintf('flipped ud \n');
            end     
            % Determine if flipping left to right is necessary
            if vcrot(posP,1)<vcrot(antP,1)
                imagedrawrot = fliplr(imagedrawrot);
                vcPSrot(:,1) = 2*ctrRot(1)-vcPSrot(:,1);
                vcrot(:,1) = 2*ctrRot(1)-vcrot(:,1);
                vdrot(:,1) = 2*ctrRot(1)-vdrot(:,1);
                vvrot(:,1) = 2*ctrRot(1)-vvrot(:,1);
    %             fprintf('flipped lr \n');
            end
            
        %     imagelogic = imagedraw>cutoff;
        %     imagedraw = double(imagedraw).*double(imagelogic);
    %         ctrWorm = vcPSrot;
    %         ctrImg = [imagex/2 imagey/2];
    %         imgshift = ctrImg - ctrWorm;
    %         imagedrawshift = imtranslate(imagedraw, [imgshift(1),imgshift(2)]);
            imagedrawcropped = imcrop(imagedrawrot,...
                [vcPSrot(1)-bodyhalflength vcPSrot(2)-bodyhalflength ...
                2*bodyhalflength 2*bodyhalflength]);
            imagesc(imagedrawcropped); % Have not added borders!!!!
            axis equal; 
            hold on;
            text(10, 10, [num2str((frm-lightstart)/(frmnum/ftl), '%.2f') 's' ], 'Color', 'w', 'fontsize', 10);
            colormap(cmp);
            c = colorbar('southoutside');
            set(c, 'xcolor', 'none', 'ycolor', 'none', 'ticklabels', '');
            cpos = get(c, 'Position');
            cpos(1) = 0.25; % x
            cpos(2) = 0.48; % y
            cpos(3) = cpos(3)/3; % w
            cpos(4) = cpos(4)/3; % h
            set(c, 'Position', cpos);        
            caxis([mini maxi]);
            set(gca, 'ticklength', [0 0], 'Color', 'k', 'XTickLabel', [], 'YTickLabel', []);
%             set(gca, 'visible', 'on');
            title(movname);
            
    %         imagex = size(imagedraw,2); imagey = size(imagedraw,1);
    %         ll = 0.1*min(imagex, imagey);
    %         line([0.9*imagex 0.9*imagex], [0.8*imagey 0.8*imagey+ll], 'color','w', 'linewidth', 1);
    %         line([0.9*imagex-ll/2 0.9*imagex+ll/2], [0.8*imagey+ll/2 0.8*imagey+ll/2], 'color','w', 'linewidth', 1);        
    %         text(0.9*imagex, 0.8*imagey+ll, 'V', 'Color', 'w', 'fontsize', 10, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
    %         text(0.9*imagex, 0.8*imagey, 'D', 'Color', 'w', 'fontsize', 10, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    %         text(0.9*imagex-ll/2, 0.8*imagey+ll/2, 'A', 'Color', 'w', 'fontsize', 10, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
    %         text(0.9*imagex+ll/2, 0.8*imagey+ll/2, 'P', 'Color', 'w', 'fontsize', 10, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
    
    %         % Ploting neurons
    %         hold on
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
                
    %     subplot(122); % Segmentation to visualize orientations
    %     hold on;
    %     plot(dorsal_data{2*frm-1,1}, dorsal_data{2*frm,1}, ':r');
    %     plot(ventral_data{2*frm-1,1}, ventral_data{2*frm,1}, ':b');
    %     plot(centerline_data_spline(:,2*frm-1), centerline_data_spline(:,2*frm), ':k');
    %     plot(centerline_data_spline(1,2*frm-1), centerline_data_spline(1,2*frm), ':og');
    %     plot(centerline_data_spline(end,2*frm-1), centerline_data_spline(end,2*frm), ':oy');
    %     set(gca, 'ydir', 'reverse');
    %     axis equal;
    
            subplot(3,2,[2 4]); % Muscle segmentations
            hold off;
            lnwidth = 1; msize = 5; headbody = 35;
            pc = plot(vcrot(:,1), vcrot(:,2), 'k', 'linewidth', lnwidth); hold on;
            plot([vdrot(headbody,1) vcrot(headbody,1)], [vdrot(headbody,2) vcrot(headbody,2)], ':k', 'linewidth', lnwidth);
            plot([vvrot(headbody,1) vcrot(headbody,1)], [vvrot(headbody,2) vcrot(headbody,2)], ':k', 'linewidth', lnwidth);
            pd = plot(vdrot(:,1), vdrot(:,2), 'color', [0 0.2 1], 'linewidth', lnwidth); 
            pv = plot(vvrot(:,1), vvrot(:,2), 'color', [1 0.2 0], 'linewidth', lnwidth);
            ph = plot(vcrot(1,1), vcrot(1,2), 'o', 'markerfacecolor', 'k', 'markeredgecolor', 'none', 'markersize', msize);
            pt = plot(vcrot(end,1), vcrot(end,2), 'o', 'markerfacecolor', 'w', 'markeredgecolor', 'k', 'markersize', msize);
            rectangle('position', ...
                [vcPSrot(1)-bodyhalflength vcPSrot(2)-bodyhalflength ...
                2*bodyhalflength 2*bodyhalflength], 'edgecolor', 'none');
    %         xpt = [centerline_data_spline(:,2*ni-1)'; dorsal_data{2*ni-1,1}'];
    %         ypt = [centerline_data_spline(:,2*ni)'; dorsal_data{2*ni,1}'];
    %         plot(xpt, ypt, ':k', 'linewidth', lnwidth/2); hold on;
    %         xpt = [centerline_data_spline(:,2*ni-1)'; ventral_data{2*ni-1,1}'];
    %         ypt = [centerline_data_spline(:,2*ni)'; ventral_data{2*ni,1}'];
    %         plot(xpt, ypt, ':k', 'linewidth', lnwidth/2);
    %         plot(dorsal_data{2*frm-1,1}, dorsal_data{2*frm,1}, 'color', [0 0.2 1], 'linewidth', lnwidth); hold on;
    %         plot(ventral_data{2*frm-1,1}, ventral_data{2*frm,1}, 'color', [1 0.2 0], 'linewidth', lnwidth);
    %         plot(centerline_data_spline(:,2*frm-1), centerline_data_spline(:,2*frm), ':k');
    %         plot(centerline_data_spline(1,2*frm-1), centerline_data_spline(1,2*frm), 'o', 'markerfacecolor', 'k', 'markeredgecolor', 'none', 'markersize', msize);
    %         plot(centerline_data_spline(end,2*frm-1), centerline_data_spline(end,2*frm), 'o', 'markerfacecolor', 'w', 'markeredgecolor', 'k', 'markersize', msize);
            legend([pd pc pv ph pt],{'Dorsal', 'Centerline', 'Ventral', 'Head', 'Tail'}); legend('boxoff');
            set(gca, 'ydir', 'reverse', 'visible', 'off'); 
    %         imagex = size(imagedraw,2); imagey = size(imagedraw,1);
    %         xlim([1 imagex]); ylim([1 imagey]);
            axis equal; % title('Segmentation');
            hold off;
            
            subplot(3,2,[5 6]); % Activity data
            hold off;
    %         maxd = max(dorsal_smd./dorsal_smd_r,[],'all');
    %         maxv = max(ventral_smd./ventral_smd_r,[],'all');
    %         maxy = 1.1*max(maxd, maxv);
    %         maxc = max(curvdataBody,[],'all');
    %         minc = min(curvdataBody,[],'all');
            bd = bar(dorsal_smd(:,frm)./dorsal_smd_r(:,frm)); hold on;
            set(bd, 'edgecolor', 'none', 'facealpha', 0.5, 'facecolor', 'b');
            bv = bar(ventral_smd(:,frm)./ventral_smd_r(:,frm)); hold on;
            set(bv, 'edgecolor', 'none', 'facealpha', 0.5, 'facecolor', 'r');
    %         plot((curvdataBody(:,frm)-minc)/maxc, ':k');
            plot([headbody headbody], [mini maxi], ':k');
    %         [~,hobj,~,~] = legend([bd bv], {'Dorsal', 'Ventral'}); legend('boxoff');
    %         hl = findobj(hobj,'type','line'); set(hl,'LineWidth',0.5);
            set(gca, 'ticklength', [0 0]); ylim([mini maxi]); 
            xlabel('Body segment'); ylabel('GCaMP/RFP'); title('Muscle activity along anterior-posterior axis');
            
            
            % Create gif file
            drawnow;
            frame = getframe(1);
            im = frame2im(frame);
            [imind, cm] = rgb2ind(im, 256);
            fullfilename = [filename_all{1,1} '_' filename_all{2,1} '_regsegfiltrat.gif'];
            if frm == 1 && nm ==1
              imwrite(imind, cm, fullfilename, 'gif', 'LoopCount', inf, 'DelayTime', ff/fps);
            else
              imwrite(imind, cm, fullfilename, 'gif', 'WriteMode', 'append', 'DelayTime', ff/fps);
            end
    
        end
    
%         save([filename '_movpara.mat'], 'tform', 'th', 'gsfilt', 'rotangle', 'ftl', 'ff', 'mini', 'maxi');
    
    end
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