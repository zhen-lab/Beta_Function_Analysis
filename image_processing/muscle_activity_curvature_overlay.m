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
    th = .25;
    % Apply threshold and filtering to all numerator images
    gsfilt = .6; 
    % gpsf = fspecial('gaussian',9,9);
    imagelist_gseg = imagelist_g;
    imagelist_rseg = imagelist_r;

    for ni = 1:length(imagelist_g)
        % Adaptive thresholding
        imusedseg = imgaussfilt(imagelist_g{ni,1}+imagelist_r{ni,1},gsfilt);
        Itoseg = uint16(imusedseg);
        T = adaptthresh(Itoseg, th);
        segmask = double(imbinarize(Itoseg, T));
        if ni==1, close all; imshow(segmask), end
        imagelist_gseg{ni,1} = segmask.*double(imgaussfilt(imagelist_g{ni,1},gsfilt));
        imagelist_rseg{ni,1} = segmask.*double(imgaussfilt(imagelist_r{ni,1},gsfilt));

    end
end

%% Test segmented images and filtering
close all;
clrtheme = 'inferno';
mini = 0; 
maxi = 2000; 
subplot(121); imagesc(imagelist_gseg{end,1}); caxis([mini maxi]); title('GCaMP'); axis equal; colormap(clrtheme);
subplot(122); imagesc(imagelist_rseg{end,1}); caxis([mini maxi]); title('RFP'); axis equal; colormap(clrtheme);

%%
close all;

if ~exist('maxi','var')
    fprintf('please test segmented images first.\n');
elseif ~exist('curvdatafiltered','var')
        fprintf('please load analyzed data first. \n');
else
    subplot(121); imagesc(imagelist_gseg{1,1});
    subplot(122); imagesc(imagelist_gseg{end,1});
    prompt = {'Frame selected',...
        'Expansion factor'};
%         'Flip? (up/down 1, left/right -1, no 0):'};
%         'Rotation (clockwise +, counterclockwise -):',...
    inputdlgtitle = 'Select frame';
    dims = [1 35];
    definput = {'1','1.5'}; % for muscle activity on Fig. 1: '438','10','2','east'
    answer = inputdlg(prompt,inputdlgtitle,dims,definput);
    frmselect = str2double(answer{1,1});
    expfactor = str2double(answer{2,1});
    frmnum = size(imagelist_g,1);        
    headbody = 35;
    cmp = colormap(inferno);
    close all;
    
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
    
    padded = 0; padsize = 0; 
    imagelist_gsegfiltuse = imagelist_gseg;
    imagelist_rsegfiltuse = imagelist_rseg;
    % Pad array if bounding box exceeds boundary
    if any(centSmoothed<bodyhalflength,'all')
        bboxdiff = ceil(bodyhalflength-centSmoothed);
        padsize = max(bboxdiff,[],'all');
        for pdi = 1:length(imagelist_gseg)
            imagelist_gsegfiltuse{pdi,1} = padarray(imagelist_gseg{pdi,1},...
                [padsize padsize],0,'both');
            imagelist_rsegfiltuse{pdi,1} = padarray(imagelist_rseg{pdi,1},...
                [padsize padsize],0,'both');
        end
        padded = 1;
        fprintf('image padded.\n');
    else
        imagelist_gsegfiltuse = imagelist_gseg;
        imagelist_rsegfiltuse = imagelist_rseg;
    end
        
    imagedrawg = imagelist_gsegfiltuse{frmselect,1}; 
    imagedrawr = imagelist_rsegfiltuse{frmselect,1}; 
%     imagex = size(imagedrawg,2); imagey = size(imagedrawg,1);
    centerpoint_frm = centSmoothed(frmselect,:) + padded*padsize;
    centerline_frm = [centerline_data_spline(:,2*frmselect-1) centerline_data_spline(:,2*frmselect)] + padded*padsize;
    dorsal_data_frm = [dorsal_data{2*frmselect-1,1} dorsal_data{2*frmselect,1}] + padded*padsize;
    ventral_data_frm = [ventral_data{2*frmselect-1,1} ventral_data{2*frmselect,1}] + padded*padsize;

    % Find angle for rotation of the selected frame
    rotangle = rotunwraprad(frmselect);
    % Rotate images
    Rimage = [cosd(-rotangle) -sind(-rotangle) 0; ...
            sind(-rotangle) cosd(-rotangle) 0; ...
            0 0 1];
    tformimage = affine2d(Rimage);
    imagedrawrotg = imwarp(imagedrawg,tformimage);
    imagedrawrotr = imwarp(imagedrawr,tformimage);
    % Rotate points (has to be negative of the angle that image rotation used)
    Rpoint = [cosd(rotangle) -sind(rotangle); ...
            sind(rotangle) cosd(rotangle)];
    ctrOrig = [size(imagedrawg,2);size(imagedrawg,1)]/2;
    ctrRot = [size(imagedrawrotg,2);size(imagedrawrotg,1)]/2;        
    vcPS = centerpoint_frm'; vcPSrot = (Rpoint*(vcPS-ctrOrig)+ctrRot)';
    vc = centerline_frm'; vcrot = (Rpoint*(vc-ctrOrig)+ctrRot)';
    vd = dorsal_data_frm'; vdrot = (Rpoint*(vd-ctrOrig)+ctrRot)';
    vv = ventral_data_frm'; vvrot = (Rpoint*(vv-ctrOrig)+ctrRot)';

    % Determine if flipping upside down is necessary
    segcent = uint16(segsize/2);
    cCenter = vcrot(segcent,:);
    vCenter = vvrot(segcent,:);
    dCenter = vdrot(segcent,:);
    vOrient = dot(vCenter-cCenter,[0;1]);
    dOrient = dot(dCenter-cCenter,[0;1]);
%         if vvrot(posP,2)<vdrot(posP,2) || vvrot(antP,2)<vdrot(antP,2) % y direction upside down, therefore ventral is on top but then shown at bottom
    if vOrient<0 || dOrient>0
        imagedrawrotg = flipud(imagedrawrotg);
        imagedrawrotr = flipud(imagedrawrotr);
        vcPSrot(:,2) = 2*ctrRot(2)-vcPSrot(:,2);            
        vcrot(:,2) = 2*ctrRot(2)-vcrot(:,2);
        vdrot(:,2) = 2*ctrRot(2)-vdrot(:,2);
        vvrot(:,2) = 2*ctrRot(2)-vvrot(:,2); 
%             fprintf('flipped ud \n');
    end     
    % Determine if flipping left to right is necessary
    if vcrot(posP,1)<vcrot(antP,1)
        imagedrawrotg = fliplr(imagedrawrotg);
        imagedrawrotr = fliplr(imagedrawrotr);
        vcPSrot(:,1) = 2*ctrRot(1)-vcPSrot(:,1);
        vcrot(:,1) = 2*ctrRot(1)-vcrot(:,1);
        vdrot(:,1) = 2*ctrRot(1)-vdrot(:,1);
        vvrot(:,1) = 2*ctrRot(1)-vvrot(:,1);
%             fprintf('flipped lr \n');
    end

    % Crop images according to new center
%     expfactor = 1.8;
    ctrCropped = [vcPSrot(1)-expfactor/2*bodyhalflength ...
        vcPSrot(2)-expfactor/2*bodyhalflength];
    ctrCroppedRep = repmat(ctrCropped,size(vcrot,1),1);
    imagedrawcroppedg = imcrop(imagedrawrotg,...
        [ctrCropped(1) ctrCropped(2) ...
        expfactor*bodyhalflength expfactor*bodyhalflength]);
    imagedrawcroppedr = imcrop(imagedrawrotr,...
        [ctrCropped(1) ctrCropped(2) ...
        expfactor*bodyhalflength expfactor*bodyhalflength]);
    % Shift the points by the new center
    vcrot = vcrot-ctrCroppedRep;
    vdrot = vdrot-ctrCroppedRep;
    vvrot = vvrot-ctrCroppedRep;
%     % Define shading area for the bodywall
%     pgon = polyshape([vvrot(headbody:end,1); vdrot(end:-1:headbody,1)],...
%         [vvrot(headbody:end,2); vdrot(end:-1:headbody,2)]);
    % Define shading area for the head
    pgon = polyshape([vcrot(headbody,1); vvrot(headbody:-1:1,1); vdrot(1:headbody,1)],...
        [vcrot(headbody,2); vvrot(headbody:-1:1,2); vdrot(1:headbody,2)]);

    % Draw figure GCaMP
    fg = figure(1);
    imagesc(imagedrawcroppedg); 
    caxis([mini maxi]); hold on;
    pg = plot(pgon);
    colormap(cmp);
    set(pg, 'facealpha', 0.7, 'edgecolor', 'none', 'facecolor', [1 1 1]);
    rectangle('position', ...
        [ctrCropped(1) ctrCropped(2) ...
        expfactor*bodyhalflength expfactor*bodyhalflength], 'edgecolor', 'none');
    imagex = expfactor*bodyhalflength; imagey = expfactor*bodyhalflength;
    ll = 0.1*min(imagex, imagey);
    line([0.9*imagex 0.9*imagex], [0.8*imagey 0.8*imagey+ll], 'color','w', 'linewidth', 1);
    line([0.9*imagex-ll/2 0.9*imagex+ll/2], [0.8*imagey+ll/2 0.8*imagey+ll/2], 'color','w', 'linewidth', 1);        
    text(0.9*imagex, 0.8*imagey+ll, 'V', 'Color', 'w', 'fontsize', 10, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
    text(0.9*imagex, 0.8*imagey, 'D', 'Color', 'w', 'fontsize', 10, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    text(0.9*imagex-ll/2, 0.8*imagey+ll/2, 'A', 'Color', 'w', 'fontsize', 10, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
    text(0.9*imagex+ll/2, 0.8*imagey+ll/2, 'P', 'Color', 'w', 'fontsize', 10, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
    set(gca, 'ydir', 'reverse', 'visible', 'off'); 
    axis square; 
    % Draw figure RFP
    fr = figure(2);
    imagesc(imagedrawcroppedr); 
    caxis([mini maxi]); hold on;
    pg = plot(pgon);
    colormap(cmp);
    set(pg, 'facealpha', 0.7, 'edgecolor', 'none', 'facecolor', [1 1 1]);
    rectangle('position', ...
        [ctrCropped(1) ctrCropped(2) ...
        expfactor*bodyhalflength expfactor*bodyhalflength], 'edgecolor', 'none');
    line([0.9*imagex 0.9*imagex], [0.8*imagey 0.8*imagey+ll], 'color','w', 'linewidth', 1);
    line([0.9*imagex-ll/2 0.9*imagex+ll/2], [0.8*imagey+ll/2 0.8*imagey+ll/2], 'color','w', 'linewidth', 1);        
    text(0.9*imagex, 0.8*imagey+ll, 'V', 'Color', 'w', 'fontsize', 10, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
    text(0.9*imagex, 0.8*imagey, 'D', 'Color', 'w', 'fontsize', 10, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    text(0.9*imagex-ll/2, 0.8*imagey+ll/2, 'A', 'Color', 'w', 'fontsize', 10, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
    text(0.9*imagex+ll/2, 0.8*imagey+ll/2, 'P', 'Color', 'w', 'fontsize', 10, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
    set(gca, 'ydir', 'reverse', 'visible', 'off'); 
    axis square; 
    
%     savefig(fg, [filename '_frm#' num2str(frmselect) '_headshaded_g']);
%     savefig(fr, [filename '_frm#' num2str(frmselect) '_headshaded_r']);   
    saveas(fg, [filename '_frm#' num2str(frmselect) '_headshaded_g.tif'], 'tiffn');
    saveas(fr, [filename '_frm#' num2str(frmselect) '_headshaded_r.tif'], 'tiffn');

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