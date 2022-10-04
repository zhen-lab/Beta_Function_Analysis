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
for i = 1:length(imagelist)
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

% Find proper adaptive threshold
th = 0.10;
Itoseg = imagelist_r{1,1};
T = adaptthresh(Itoseg, th);
Iseg = imbinarize(Itoseg, T); 
imshow(Iseg);

% Apply threshold and filtering to all numerator images
gsfilt = .6; 
% gpsf = fspecial('gaussian',9,9);
imagelist_gseg = imagelist_g;
imagelist_segratio = imagelist_g;
imagelist_segratiofilt = imagelist_g;
% imagelist_segratiodecon = imagelist_g;

for ni = 1:length(imagelist)
    segmask = uint16(imbinarize(imagelist_r{ni,1}, T));
    imagelist_gseg{ni,1} = segmask.*imagelist_g{ni,1};
    imagelist_segratio{ni,1} = imagelist_gseg{ni,1}./imagelist_r{ni,1};
    imagelist_segratiofilt{ni,1} = imgaussfilt(imagelist_segratio{ni,1},gsfilt);
%     imagelist_gsegdecon = deconvlucy(imagelist_gseg{ni,1},gpsf);
%     imagelist_segratiofilt{ni,1} = imagelist_gseggfilt./imagelist_r{ni,1};
%     imagelist_segratiodecon{ni,1} = imagelist_gsegdecon./imagelist_r{ni,1};
end

% Test segmented images and filtering
close all;
clrtheme = 'inferno';
mini = 0; 
maxi = 4; 
subplot(131); imagesc(imagelist_g{end,1}); caxis([0 1000]); title('GCaMP'); axis equal; colormap(clrtheme);
subplot(132); imagesc(imagelist_segratio{end,1}); caxis([mini maxi]); title('Segmented ratio'); axis equal; colormap(clrtheme);
subplot(133); imagesc(imagelist_segratiofilt{end,1}); caxis([mini maxi]); title('Segmented filtered ratio'); axis equal; colormap(clrtheme);

%% Write gif file

close all;
mini = 0; 
maxi = 1.5; 
frmnum = length(imagelist); 
imstr = 1; 
imend = frmnum;
ftl = 90;
fps = frmnum/ftl; 
ff = 0.3;
clrtheme = 'inferno';

for frm = imstr:imend
    
    figure(1);
%     set(gcf,'units','normalized','outerposition',[0 0 1 1], 'color', 'k');
    set(gcf, 'color', 'w', 'pos', [10 10 500 600]);
    hold off;
%     imagedraw = imagelist{frm,1};
%     imagedraw = imrotate((imagelist{frm,1}), rotangle, 'bilinear', 'crop');
%     imagelogic = imagedraw>cutoff;
    imagedraw = imagelist_segratiofilt{frm,1};
    imagesc(imagedraw); % Have not added borders!!!!
    hold on;
    text(10, 10, [num2str((frm-imstr)/(frmnum/ftl), '%.2f') 's' ], 'Color', 'w', 'fontsize', 10);
    
    axis equal;
    colormap(clrtheme);
    c = colorbar('southoutside');
    set(c, 'xcolor', 'none', 'ycolor', 'none', 'ticklabels', '');
    cpos = get(c, 'Position');
%     cpos(1) = 0.41; cpos(3) = cpos(3)/3.5; 
    cpos(1) = 0.41; cpos(3) = cpos(3)/3;
    cpos(4) = cpos(4)/3;    
    set(c, 'Position', cpos);        
    caxis([mini maxi]);
    set(gca, 'visible', 'off');
    
    hold on
    
    % Create gif file
    drawnow;
    frame = getframe(1);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);
    if frm == imstr
      imwrite(imind, cm, [filename(1:end-4) '_regsegfiltrat.gif'], 'gif', 'LoopCount', inf, 'DelayTime', ff/fps);
    else
      imwrite(imind, cm, [filename(1:end-4) '_regsegfiltrat.gif'], 'gif', 'WriteMode', 'append', 'DelayTime', ff/fps);
    end

end

save([filename '_movpara.mat'], 'tform', 'th', 'gsfilt', 'ftl', 'ff');

%% Write registered tiff file

for i = 1:length(imagelist)
    if i==1
        imwrite(imagelist_segratiofilt{i,1}, [filename(1:end-4) '_regsegfiltrat.tif']);
    else
        imwrite(imagelist_segratiofilt{i,1}, [filename(1:end-4) '_regsegfiltrat.tif'], 'writemode', 'append');
    end
end

fprintf('segmented filtered ratio channel tiff file saved. \n');

%% Write registered tiff file
for i = 1:length(imagelist)
    if i==1
        imwrite(imagelist_r{i,1}, [filename(1:end-4) '_regred.tif']);
    else
        imwrite(imagelist_r{i,1}, [filename(1:end-4) '_regred.tif'], 'writemode', 'append');
    end
end

fprintf('red channel tiff file saved. \n');

%% Write registered tiff file
for i = 1:length(imagelist)
    if i==1
        imwrite(imagelist_gseg{i,1}, [filename(1:end-4) '_reggreen.tif']);
    else
        imwrite(imagelist_gseg{i,1}, [filename(1:end-4) '_reggreen.tif'], 'writemode', 'append');
    end
end

fprintf('green segmented channel tiff file saved. \n');