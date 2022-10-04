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
    movpara.channeltomove = questdlg('Move left or right channel?',...
        'Channel to move', 'Left', 'Right', ...
        'Right');
    
    % Register RFP and GFP channels
    switch movpara.channeltomove
        case 'Left'
            figure;
            subplot(1,2,1); imshowpair(imagelist_r{1}, imagelist_g{1});
            [img_updated, movingRegistered, tform] = image_registration_tform_muscle...
                (imagelist_g, imagelist_r, movpara.channeltomove);
            subplot(1,2,2); imshowpair(movingRegistered{1}, imagelist_g{1});
        case 'Right'  
            figure;                  
            subplot(1,2,1); imshowpair(imagelist_r{1}, imagelist_g{1});
            [img_updated, movingRegistered, tform] = image_registration_tform_muscle...
                (imagelist_r, imagelist_g, movpara.channeltomove);
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
switch movpara.channeltomove
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
    movpara.th = .25;
    % Apply threshold and filtering to all numerator images
    movpara.gsfilt = .6; 
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
%         imagelist_segratiofilt{ni,1} = imgaussfilt(imagelist_segratio{ni,1},movpara.movpara.gsfilt);
        % Adaptive thresholding
        imusedseg = imgaussfilt(imagelist_g{ni,1}+imagelist_r{ni,1},movpara.gsfilt);
        Itoseg = uint16(imusedseg);
        T = adaptthresh(Itoseg, movpara.th);
        segmask = double(imbinarize(Itoseg, T));
        if ni==1, close all; imshow(segmask), end
        imagelist_gseg{ni,1} = segmask.*double(imgaussfilt(imagelist_g{ni,1},movpara.gsfilt));
        imagelist_segratio{ni,1} = imagelist_gseg{ni,1}./imagelist_r{ni,1};
        imagelist_segratiofilt{ni,1} = imgaussfilt(imagelist_segratio{ni,1},movpara.gsfilt);
%         % Deconvolution
%         imagelist_gsegdecon = deconvlucy(imagelist_gseg{ni,1},gpsf);
%         imagelist_segratiofilt{ni,1} = imagelist_gseggfilt./imagelist_r{ni,1};
%         imagelist_segratiodecon{ni,1} = imagelist_gsegdecon./imagelist_r{ni,1};
    end
end

%% Test segmented images and filtering

close all;
clrtheme = 'inferno';
movpara.mini = 0; 
movpara.maxi = 2; 
subplot(131); imagesc(imagelist_g{end,1}); caxis([0 3000]); title('GCaMP'); axis equal; colormap(clrtheme);
subplot(132); imagesc(imagelist_segratio{end,1}); caxis([movpara.mini movpara.maxi]); title('Segmented ratio'); axis equal; colormap(clrtheme);
subplot(133); imagesc(imagelist_segratiofilt{end,1}); caxis([movpara.mini movpara.maxi]); title('Segmented filtered ratio'); axis equal; colormap(clrtheme);

%% Loop for animation

close all;
if ~exist('movpara','var')
    fprintf('please test segmented images first.\n');
elseif ~exist('np','var')
        fprintf('please load analyzed data first. \n');
elseif size(np,2)/2==1
    fprintf('please choose a sample with more than one neuron analyzed.\n');
else
    canvas_x = 3;
    canvas_y = size(np,2)/2;
    neurontotal = canvas_y;
    ratio = gfp./rfp;
    typeRecord = questdlg(...
        'All-optical or calcium only recording?',...
        'All-optical/Calcium-only',...
        'All-optical', 'Calcium-only', 'Calcium-only');
    switch typeRecord
        case 'All-optical'
            ratio_norm = ratio./repmat(ratio(1,:),size(ratio,1),1);
        case 'Calcium-only'
            ratio_norm = ratio./repmat(max(ratio),size(ratio,1),1);
            forOrback = questdlg('Highligh backward or forward?',...
                'Backward or forward',...
                'Backward','Forward','Backward');
    end
    
    subplot(121); imagesc(imagelist_segratiofilt{1,1}); 
    subplot(122); imagesc(imagelist_segratiofilt{end,1}); 
    prompt = {...
        'Total time length (second):',...
        'Total frame number:', ...
        'Play speed adjustment:',...
        'Flip upside down? (Yes=1):',...
        'Flip left to right? (Yes=1):',...
        'Name of neurons:',...
        'Starting neuron number:',...
        'Radius of tracker:',...
        'Highlight boundaries, enter [start end]:'};
%         'Lower bound of activity:',...
%         'Upper bound of activity:'};
    inputdlgtitle = 'Parameters for GIF';
    dims = [1 35];
    definput = {'90',num2str(size(imagelist_g,1)),...
        '0.3','0','0',...
        'DB','5',num2str(7*ones(1,neurontotal)),...
        ''}; % type in backward frames boundaries, with no colons
    answer = inputdlg(prompt,inputdlgtitle,dims,definput);
    movpara.flt = str2double(answer{1,1});
    movpara.frmnum = str2double(answer{2,1});
    fps = movpara.frmnum/movpara.flt;
    movpara.ff = str2double(answer{3,1});
    movpara.flpud = str2double(answer{4,1});
    movpara.flplr = str2double(answer{5,1});
    movpara.neuronname = answer{6,1}; 
    movpara.neuronstart = str2double(answer{7,1}); 
    neuronseries = movpara.neuronstart:movpara.neuronstart+canvas_y-1;
%     rotangle = str2double(answer{5,1}); % clockwise: positive value, counterclosewise: negative value
    movpara.r = str2num(answer{8,1}); % radius of tracker
    movpara.backbound = str2num(answer{9,1});
    backboundReshape = reshape(movpara.backbound,2,[])'; 
    backeventnum = size(backboundReshape,1);

%     actLow = str2double(answer{9,1}); 
%     actTop = str2double(answer{10,1});
    % backperiod = [109 193];
    % neuronnames = {'AVB', 'AVA/AVE'}; neuronnum = size(neuronnames, 2);
    % movpara.mini = 0; 
    % movpara.maxi = 1; 
    % movpara.flt = 90;
    % fps = 1800/movpara.flt; 
    % movpara.ff = 0.3;
    lightstart = 1; 
    lightend = length(imagelist_g);
    cmp = colormap(inferno); % cmp(1,:) = [0 0 0];
    cmpnn = lbmap(canvas_y, 'redblue');
    % cutoff = mean2(imagelist{1,1});
    % rotangle = 0; % clockwise: positive value, counterclosewise: negative value
    close all;
    
    % Restore neuron positions tracked on the right channel
    switch framestruct
        case 'Split'
            if any(np(:,1:2:end)>size(img_stack,2)/2)
                np(:,1:2:end) = np(:,1:2:end)-size(img_stack,2)/2;
            end
        case 'Alternating'
            if any(np(:,1:2:end)>size(img_stack,2))
                np(:,1:2:end) = np(:,1:2:end)-size(img_stack,2);
            end
    end

    % Determine rotation angle
    frmtotal = size(np,1);
    vectotal = uint16(size(np,2)/2)-1;
    npX = diff(np(:,1:2:end),1,2);
    npY = diff(np(:,2:2:end),1,2);
%     Mx = np(:,1:2:end); nk = nchoosek(1:size(Mx,2), 2);
%     npX = Mx(:,nk(:,1)) - Mx(:,nk(:,2));
%     My = np(:,2:2:end); % nk = nchoosek(1:size(My,2), 2);
%     npY = My(:,nk(:,1)) - My(:,nk(:,2));
    thetaNeuronRadAll = npX;
    for numvec = 1:vectotal
        for numfrm = 1:frmtotal
%             thetaNeuronRadAll(numfrm,numvec) = ...
%                 acos(dot([npX(numfrm,numvec) npY(numfrm,numvec)],[1 0])...
%                 /sqrt(npX(numfrm,numvec)^2+npY(numfrm,numvec)^2));
            thetaNeuronRadAll(numfrm,numvec) = ...
                atan2(npY(numfrm,numvec),npX(numfrm,numvec));
        end
    end
    % First unwrap angles, then average for all
    stafctr = 5;
    thetaNeuronRadAllMean = mean(unwrap(thetaNeuronRadAll),2);
    thetaNeuronRadAllMeanSmDeg = rad2deg(smoothdata(...
        thetaNeuronRadAllMean,1,'movmedian',uint16(frmtotal/stafctr)));
    
    
    % Stablize the center which is the mean of X and Y coordinates
    centX = mean(np(:,1:2:end),2);
    centY = mean(np(:,2:2:end),2);
        centSmoothed = smoothdata([centX centY],...
        1,'movmedian',uint16(frmtotal/stafctr));
    % Calculate the maximum distance between center and neurons
    diffcnX = (repmat(centSmoothed(:,1),1,neurontotal)-np(:,1:2:end)).^2;
    diffcnY = (repmat(centSmoothed(:,2),1,neurontotal)-np(:,2:2:end)).^2;
    maxlengthCN = sqrt(max(diffcnX+diffcnY,[],'all'));
    % Calculate the maximum length of connected neurons
    diffnnX = diff(np(:,1:2:end),1,2).^2;
    diffnnY = diff(np(:,2:2:end),1,2).^2;
    maxlengthNN = sqrt(max(sum(diffnnX+diffnnY,2)));
    
    % Pad array to avoid bounding box bleeding
    % Pad size is defined by longer of the two plus tracker size:
    % 1 the longest distance between center and neurons
    % 2 the longest connected segments between neurons
    padsize = ceil(max(maxlengthCN,maxlengthNN)+max(movpara.r));
    imagelist_segratiofiltpad = imagelist_segratiofilt;
    for pdi = 1:length(imagelist_segratiofilt)
        imagelist_segratiofiltpad{pdi,1} = ...
            padarray(imagelist_segratiofilt{pdi,1},...
            [padsize padsize],0,'both');
    end
    % Update coordinates for center and neurons
    centSmoothed = centSmoothed + padsize;
    nppad = np + padsize;
        
%     % Stablize the center which is the mean of X and Y coordinates
%     nppad = np + padsize;
%     centX = mean(nppad(:,1:2:end),2);
%     centY = mean(nppad(:,2:2:end),2);
%     centSmoothed = smoothdata([centX centY],...
%         1,'movmedian',uint16(frmtotal/stafctr));
    
    for frm = lightstart:lightend

        figure(1); 
        
%%%%%%%%% Actual recording, first column of figure %%%%%%%%%%%%%
        subplot(canvas_y,canvas_x,1:canvas_x:canvas_x*canvas_y); 
        set(gcf,'units','normalized','outerposition',[0 0 1 0.75]);
        hold off;
        imagedraw = imagelist_segratiofiltpad{frm,1}; 
        imagex = size(imagedraw,2); imagey = size(imagedraw,1);
        maxlength = padsize; % longest length of connected segments between neurons
        centerpoint_frm = centSmoothed(frm,:);
        np_frm = [nppad(frm,1:2:end); nppad(frm,2:2:end)]'; 

        % Rotate image counterclockwise
        rotangle_frm = thetaNeuronRadAllMeanSmDeg(frm);
        imagedrawrot = imrotate(imagedraw, rotangle_frm, 'bilinear');
        % Points are rotated opposite as images are upside down
        R = [cosd(-rotangle_frm) -sind(-rotangle_frm);...
            sind(-rotangle_frm) cosd(-rotangle_frm)];
        ctrOrig = [size(imagedraw,2);size(imagedraw,1)]/2;
        ctrRot = [size(imagedrawrot,2);size(imagedrawrot,1)]/2;        
        vcPS = centerpoint_frm'; vcPSrot = (R*(vcPS-ctrOrig)+ctrRot)';
        vn = np_frm'; vnrot = (R*(vn-ctrOrig)+ctrRot)';
        
        % Determine if flipping upside down is necessary
        if logical(movpara.flpud)
            imagedrawrot = flipud(imagedrawrot);
            vcPSrot(:,2) = 2*ctrRot(2)-vcPSrot(:,2);
            vnrot(:,2) = 2*ctrRot(2)-vnrot(:,2);
        end
        % Determine if flipping left to right is necessary
        if logical(movpara.flplr)
            imagedrawrot = fliplr(imagedrawrot);
            vcPSrot(:,1) = 2*ctrRot(1)-vcPSrot(:,1);
            vnrot(:,1) = 2*ctrRot(1)-vnrot(:,1);
        end

        % Crop image by centering around the center of trackers
        imagedrawcropped ...
            = imcrop(imagedrawrot,...
            [vcPSrot(1)-maxlength vcPSrot(2)-maxlength ...
            2*maxlength 2*maxlength]);
        % Update the coordinates for trackers
        vnrotcropped = vnrot-repmat(vcPSrot-maxlength,size(vnrot,1),1);
        
        % Draw annotated recording after rotation+cropping
        imagesc(imagedrawcropped);
        axis equal; 
        hold on;
        for nn = 1:canvas_y
            radnn = movpara.r(nn);
            rectangle('Curvature', [0,0], ...
                'Position',...
                [vnrotcropped(nn,1)-radnn, vnrotcropped(nn,2)-radnn,...
                2*radnn, 2*radnn],...
                'EdgeColor', cmpnn(nn,:), 'Linestyle', '-', 'Linewidth', 1);
%             text(vnrotcropped(nn,1), vnrotcropped(nn,2)+movpara.r,...
%                 num2str(neuronseries(nn)),...
%                 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top',...
%                  'Color', 'w', 'fontsize', 10);
        end
        if isequal(typeRecord,'Calcium-only')
%             if mod(length(movpara.backbound),2)~=0
%                 break
%             else
                % Creat a set for backing frames
                bounddiff = diff(backboundReshape,1,2)+1;
                backfrm = zeros(sum(bounddiff),1);
                startp = 1;
                for bnum = 1:backeventnum
                    backchunk = backboundReshape(bnum,1):backboundReshape(bnum,2);
                    lengthchunk = length(backchunk);
                    backfrm(startp:(startp+lengthchunk-1)) = backchunk;
                    startp = startp+lengthchunk;
                end
                % Check if the current frame fall in the backing frames
                arwht = 0.5;
                if ismember(frm,backfrm)
                    delete(findall(gcf,'type','annotation'));
                    switch forOrback
                        case 'Backward'
                            annotation('arrow', [0.35 0.37], [arwht arwht], 'linewidth', 10, 'headwidth', 20, 'color', 'k');                          
                        case 'Forward'
                            annotation('arrow', [0.12 0.10], [arwht arwht], 'linewidth', 10, 'headwidth', 20, 'color', 'k');                            
                    end  
                else
                    delete(findall(gcf,'type','annotation'));
                    switch forOrback
                        case 'Backward'
                            annotation('arrow', [0.12 0.10], [arwht arwht], 'linewidth', 10, 'headwidth', 20, 'color', 'k');                            
                        case 'Forward'
                            annotation('arrow', [0.35 0.37], [arwht arwht], 'linewidth', 10, 'headwidth', 20, 'color', 'k');                          
                    end  
                end
        end
%         plot(vnrot(:,1),vnrot(:,2),'ow','markersize',movpara.r);
        text(20, 0, [num2str((frm-lightstart)/(frmtotal/movpara.flt), '%.2f') 's' ], 'Color', 'w', 'fontsize', 10);
        colormap(cmp);
        c = colorbar('southoutside');
        set(c, 'xcolor', 'none', 'ycolor', 'none', 'ticklabels', '');
        cpos = get(c, 'Position');
        cpos(1) = 0.20; % x
%         cpos(2) = 0.08; % y
        cpos(3) = cpos(3)/3; % w
        cpos(4) = cpos(4)/3; % h
        set(c, 'Position', cpos);        
        caxis([movpara.mini movpara.maxi]);
        set(gca, 'xtick', [], 'ytick', [], 'color', 'k');
%%%%%%%%% End of first column of figure %%%%%%%%%%%%%
        
%%%%%%%%% Zoomed-in neurons, second column of figure %%%%%%%%%%%%%
        for pn = 1:canvas_y
            
            radpn = movpara.r(pn);
            lnwd = 2;
            subplot(canvas_y,canvas_x,2+(pn-1)*canvas_x);
            imageN = imcrop(imagedrawcropped,...
                [vnrotcropped(pn,1)-radpn vnrotcropped(pn,2)-radpn...
                2*radpn 2*radpn]);
            imagesc(imageN);
            caxis([movpara.mini movpara.maxi]); 
            axis square;
            set(gca, 'xticklabel', '', 'yticklabel', '',...
                'xtick', [], 'ytick', [], 'ydir', 'reverse',...
                'xcolor', cmpnn(pn,:), 'ycolor', cmpnn(pn,:), 'linewidth', lnwd);  
            if ~isequal(movpara.neuronname,'DD')
                ylabel([movpara.neuronname num2str(neuronseries(pn))],'color','k','fontsize',10);   
            else
                ylabel(['NMJ' num2str(neuronseries(pn)) ...
                    ' (' movpara.neuronname num2str(neuronseries(pn)) ')'],...
                    'color','k','fontsize',10);   
            end
%%%%%%%%% End of second column of figure %%%%%%%%%%%%%

%%%%%%%%% Activity of neurons, third column of figure %%%%%%%%%%%%%
            subplot(canvas_y,canvas_x,pn*canvas_x);
            hold off;            
            yupbound = 1.1*max(ratio_norm,[],'all');
            switch typeRecord
                case 'All-optical'
                    plot([0 numfrm],[1 1],':k'); hold on;
                case 'Calcium-only'
                    for cn = 1:backeventnum
                        f = fill([backboundReshape(cn,1) backboundReshape(cn,2)...
                            backboundReshape(cn,2) backboundReshape(cn,1)],...
                            [0 0 yupbound yupbound], 0.8*[1 1 1]);
                        set(f, 'edgecolor', 'none', 'facealpha', 0.5);
                        hold on;
                    end
                    fw = fill([frm numfrm numfrm frm],...
                        [0 0 yupbound yupbound], [1 1 1]);
                    set(fw, 'edgecolor', 'none');
            end
            plot(ratio_norm(1:frm,pn),'color',cmpnn(pn,:),'linewidth',lnwd);
            xlim([1 numfrm]); ylim([0 yupbound]);
            numxtick = 4;
            xtickseries = floor(linspace(1,numfrm,numxtick));
            xticklabelseries = string(round((xtickseries-1)/fps,1));
            ytickseries = 0:1:yupbound;
            yticklabelseries = string(ytickseries);
            if pn==canvas_y
                set(gca,...
                    'xtick', xtickseries, 'ytick', ytickseries,...
                    'xticklabel', xticklabelseries, 'yticklabel', yticklabelseries,...
                    'ticklength', [0 0]);                
            else
                set(gca,...
                    'xtick', [], 'ytick', ytickseries,...
                    'xticklabel', [], 'yticklabel', yticklabelseries,...
                    'ticklength', [0 0]);    
            end   
            
        end
%%%%%%%%% End of third column of figure %%%%%%%%%%%%%

%         hold on
    %     for nn = 1:neuronnum
    %         plot(dual_position_data{frm,nn}(1)-n,dual_position_data{frm,nn}(2),'ow', 'markersize', 4*movpara.r);
    %         text(dual_position_data{frm,nn}(1)-n-2*movpara.r, dual_position_data{frm,nn}(2)+2*movpara.r, neuronnames(nn), 'color', 'w', 'fontsize', 10);
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
%     plot(dorsal_data{2*frm-1,1}, dorsal_data{2*frm,1}, ':movpara.r');
%     plot(ventral_data{2*frm-1,1}, ventral_data{2*frm,1}, ':b');
%     plot(centerline_data_spline(:,2*frm-1), centerline_data_spline(:,2*frm), ':k');
%     plot(centerline_data_spline(1,2*frm-1), centerline_data_spline(1,2*frm), ':og');
%     plot(centerline_data_spline(end,2*frm-1), centerline_data_spline(end,2*frm), ':oy');
%     set(gca, 'ydir', 'reverse');
%     axis equal;

%         subplot(3,2,[2 4]); % Muscle segmentations
%         hold off;
%         lnwidth = 1; msize = 5; headbody = 35;
%         pc = plot(vcrot(:,1), vcrot(:,2), 'k', 'linewidth', lnwidth); hold on;
%         plot([vdrot(headbody,1) vcrot(headbody,1)], [vdrot(headbody,2) vcrot(headbody,2)], ':k', 'linewidth', lnwidth);
%         plot([vvrot(headbody,1) vcrot(headbody,1)], [vvrot(headbody,2) vcrot(headbody,2)], ':k', 'linewidth', lnwidth);
%         pd = plot(vdrot(:,1), vdrot(:,2), 'color', [0 0.2 1], 'linewidth', lnwidth); 
%         pv = plot(vvrot(:,1), vvrot(:,2), 'color', [1 0.2 0], 'linewidth', lnwidth);
%         ph = plot(vcrot(1,1), vcrot(1,2), 'o', 'markerfacecolor', 'k', 'markeredgecolor', 'none', 'markersize', msize);
%         pt = plot(vcrot(end,1), vcrot(end,2), 'o', 'markerfacecolor', 'w', 'markeredgecolor', 'k', 'markersize', msize);
%         rectangle('position', ...
%             [vcPSrot(1)-bodyhalflength vcPSrot(2)-bodyhalflength ...
%             2*bodyhalflength 2*bodyhalflength], 'edgecolor', 'none');
% %         xpt = [centerline_data_spline(:,2*ni-1)'; dorsal_data{2*ni-1,1}'];
% %         ypt = [centerline_data_spline(:,2*ni)'; dorsal_data{2*ni,1}'];
% %         plot(xpt, ypt, ':k', 'linewidth', lnwidth/2); hold on;
% %         xpt = [centerline_data_spline(:,2*ni-1)'; ventral_data{2*ni-1,1}'];
% %         ypt = [centerline_data_spline(:,2*ni)'; ventral_data{2*ni,1}'];
% %         plot(xpt, ypt, ':k', 'linewidth', lnwidth/2);
% %         plot(dorsal_data{2*frm-1,1}, dorsal_data{2*frm,1}, 'color', [0 0.2 1], 'linewidth', lnwidth); hold on;
% %         plot(ventral_data{2*frm-1,1}, ventral_data{2*frm,1}, 'color', [1 0.2 0], 'linewidth', lnwidth);
% %         plot(centerline_data_spline(:,2*frm-1), centerline_data_spline(:,2*frm), ':k');
% %         plot(centerline_data_spline(1,2*frm-1), centerline_data_spline(1,2*frm), 'o', 'markerfacecolor', 'k', 'markeredgecolor', 'none', 'markersize', msize);
% %         plot(centerline_data_spline(end,2*frm-1), centerline_data_spline(end,2*frm), 'o', 'markerfacecolor', 'w', 'markeredgecolor', 'k', 'markersize', msize);
%         legend([pd pc pv ph pt],{'Dorsal', 'Centerline', 'Ventral', 'Head', 'Tail'}); legend('boxoff');
%         set(gca, 'ydir', 'reverse', 'visible', 'off'); 
% %         imagex = size(imagedraw,2); imagey = size(imagedraw,1);
% %         xlim([1 imagex]); ylim([1 imagey]);
%         axis equal; % title('Segmentation');
%         hold off;
        
%         subplot(3,2,[5 6]); % Activity data
%         hold off;
% %         maxd = max(dorsal_smd./dorsal_smd_r,[],'all');
% %         maxv = max(ventral_smd./ventral_smd_r,[],'all');
% %         maxy = 1.1*max(maxd, maxv);
% %         maxc = max(curvdataBody,[],'all');
% %         minc = min(curvdataBody,[],'all');
%         bd = bar(dorsal_smd(:,frm)./dorsal_smd_r(:,frm)); hold on;
%         set(bd, 'edgecolor', 'none', 'facealpha', 0.5, 'facecolor', 'b');
%         bv = bar(ventral_smd(:,frm)./ventral_smd_r(:,frm)); hold on;
%         set(bv, 'edgecolor', 'none', 'facealpha', 0.5, 'facecolor', 'movpara.r');
% %         plot((curvdataBody(:,frm)-minc)/maxc, ':k');
%         plot([headbody headbody], [movpara.mini movpara.maxi], ':k');
% %         [~,hobj,~,~] = legend([bd bv], {'Dorsal', 'Ventral'}); legend('boxoff');
% %         hl = findobj(hobj,'type','line'); set(hl,'LineWidth',0.5);
%         set(gca, 'ticklength', [0 0]); ylim([movpara.mini movpara.maxi]); 
%         xlabel('Body segment'); ylabel('GCaMP/RFP'); title('Muscle activity along anterior-posterior axis');
        
        
        % Create gif file
        drawnow;
        frame = getframe(1);
        im = frame2im(frame);
        [imind, cm] = rgb2ind(im, 256);
        if frm == lightstart
          imwrite(imind, cm, [filename '_regsegfiltrat.gif'], 'gif', 'LoopCount', inf, 'DelayTime', movpara.ff/fps);
        else
          imwrite(imind, cm, [filename '_regsegfiltrat.gif'], 'gif', 'WriteMode', 'append', 'DelayTime', movpara.ff/fps);
        end

    end

    save([filename '_movpara.mat'], 'tform', 'thetaNeuronRadAll', ....
        'movpara');

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