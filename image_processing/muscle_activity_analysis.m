%% Clear everything

clear; close all;

%% Load tiff files

setup_proof_reading_split_alternate;
updated = 0;

% %% Check which channel needs to move 
% % for analyzed data!!!!!!!!!!!!!!
% 
% if exist('dorsal_data', 'var')
%     close all;
%     figure;
%     subplot(121);
%     imshow(imagelist_r{1}); caxis('auto');
%     hold on
%     plot(dorsal_data{2*1-1,1}, dorsal_data{2*1,1}, 'r');
%     plot(ventral_data{2*1-1,1}, ventral_data{2*1,1}, 'b');
%     subplot(122);%     title('RFP/Left');

%     imshow(imagelist_g{1}); caxis('auto');
%     hold on
%     plot(dorsal_data{2*1-1,1}, dorsal_data{2*1,1}, 'r');
%     plot(ventral_data{2*1-1,1}, ventral_data{2*1,1}, 'b');
%     title('GCaMP/Right');
% else
%     fprintf('You have not loaded analyzed data yet. \n');
% end

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

% %% (Individual Data Struct) Check using analyzed data 
% % For those processed before calculated registration
% % Load data first!!!!!!!!!
% 
% hold on;
% plot(dorsal_data{1,1}, dorsal_data{2,1}, 'r');
% plot(ventral_data{1,1}, ventral_data{2,1}, 'b');
% fprintf('contour overlaid. \n');

%% Update the channel that has been moved

% Remove edges from registration
imagelist_moved = movingRegistered;
frmnum = size(imagelist_g,1);
for i = 1:frmnum
    img = movingRegistered{i};
    img = img + mean(img,[1 2])*uint16(img==0);
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

% Overlay two channels then top-bottom hat operations if necessary
imagelist_use = imagelist_moved;
disksize = 3;
for i = 1:frmnum
    imgoverlaid = imagelist_g{i,1}+imagelist_r{i,1};
    se = strel('disk',disksize);
    imagelist_use{i,1} = ...
        imsubtract(imadd(imgoverlaid,imtophat(imgoverlaid,se)),imbothat(imgoverlaid,se));
end
subplot(121); imagesc(imagelist_g{1,1}+imagelist_r{1,1}); axis equal;
subplot(122); imagesc(imagelist_use{1,1}); axis equal;
fprintf('channels overlaid and top-bottom hat operated. \n');

%% Remove contact area manually if necessary

% imagelist_g = movingRegistered;
% 
% % Select channel
% channel = questdlg('Use GFP or RFP channel?', 'Channel', 'GFP', 'RFP', 'Both', 'RFP');
% switch channel
%     case 'GFP'
%         imagelist_use = imagelist_g;
%     case 'RFP'
%         imagelist_use = imagelist_r;
%     case 'Both'
%         imagelist_gr = imagelist_g;
%         for i = 1:size(imagelist_gr, 1)
%             imagelist_gr{i,1} = imagelist_g{i,1}+imagelist_r{i,1};
%         end
%         imagelist_use = imagelist_gr;
% end

% Draw contact area
mm = 1.1*max(imagelist_use{1,1}, [], 'all');
figure;
subplot(121); imshow(imagelist_use{1,1}); 
title('Before removing contact area');
caxis([0 mm]);
h = drawpolygon;
fprintf('polygon is drawn. \n');

% Obtain contact area
pos = h.Position;
col_num = size(imagelist_use{1,1}, 1);
row_num = size(imagelist_use{1,1}, 2);
xq = 1:row_num; xq_rep = repmat(xq, col_num, 1); xq_rep_shape = reshape(xq_rep, [], 1);
yq = 1:col_num; yq_rep = repmat(yq, row_num, 1); yq_rep_shape = reshape(yq_rep', [], 1);
in = inpolygon(xq_rep_shape, yq_rep_shape, pos(:,1), pos(:,2));
logic_all = reshape(in, col_num, row_num);
image_in = uint16(logic_all).*imagelist_use{1,1};
logic = image_in > 0;
imagelist_updated = imagelist_use;

% Remove contact area
for i = 1:size(imagelist_use,1)
    img = imagelist_use{i,1};
    pad = mean(img,[1 2])*ones(col_num, row_num);
    imgrmvd = uint16(~(logic)).*img;
    se = strel('disk',5);
    imagelist_updated{i,1} = ...
        imsubtract(imadd(imgrmvd,imtophat(imgrmvd,se)),imbothat(imgrmvd,se));
%     pad_noise = imnoise(pad, 'gaussian');
%     imagelist_updated{i,1} = uint16(pad.*pad_noise.*double(logic)) + uint16(~(logic)).*img;
end
subplot(122); imshow(imagelist_updated{1,1}); 
title('After removing contact area');
caxis([0 mm]);

% Update images to use later
imagelist_use = imagelist_updated;
fprintf('images are updated. \n');

%% Process data

% Delineate dorsal and ventral muscles
frmtotal = size(imagelist_g,1);
prompt = {'First frame', ['Last frame (total of ' num2str(frmtotal) ' frames)']};
inputdlgtitle = 'Range for analysis';
dims = [1 35];
definput = {'1',num2str(frmtotal)};
frmrangecell = inputdlg(prompt,inputdlgtitle,dims,definput);
frmrange = str2double(frmrangecell(1)):str2double(frmrangecell(2));

if updated==1
    
    tic; close all;
    [dorsal_data, ventral_data, ...
        centerline_data, centerline_data_spline, ...
        curvdata, curvdatafiltered] = extract_centerline_vd(...
        imagelist_use(frmrange), 6, 4, 3);
    toc;
    
    fprintf('segmentation completed. \n');
    
else
    
    fprintf('please update channel after registration. \n');

end

%% Calculate dorsal and ventral muscle activities, save output

if updated==1
    
    % Generate data
    tic;
    fprintf('activity analysis kicks off \n');
    [dorsal_smd, ventral_smd, dorsal_smd_r, ventral_smd_r, ...
        dorsal_raw, ventral_raw, dorsal_raw_r, ventral_raw_r] = ...
        activity_all(imagelist_g(frmrange), imagelist_r(frmrange),...
        [1, length(frmrange)], dorsal_data, ventral_data, ...
        centerline_data_spline, curvdatafiltered);
    fprintf('activity analysis finished. \n');
    figure;
    subplot(1,2,1); imagesc(dorsal_smd./dorsal_smd_r); title('Dorsal');
    subplot(1,2,2); imagesc(ventral_smd./ventral_smd_r); title('Ventral');
    toc

    % Save data
    parts = strsplit(pathname, '\');
    masterfolder = find(cellfun(@(x) strcmp(x,'L1'), parts, 'UniformOutput', 1));
    data_path = fullfile(parts{1,1:masterfolder}, 'Alpha_Data_Raw', parts{1,masterfolder+2:end});
    warning('off'); mkdir(data_path); 
    data_path_name = fullfile(data_path, [filename(1:end-4) '.mat']);
    save(data_path_name, ...
        'centerline_data_spline', 'curvdatafiltered', ...
        'dorsal_data', 'dorsal_smd', 'dorsal_smd_r', 'dorsal_raw', 'dorsal_raw_r', ...
        'ventral_data', 'ventral_smd', 'ventral_smd_r', 'ventral_raw', 'ventral_raw_r', ...
        'tform', 'channeltomove');
    fprintf('data saved. \n');

%     system('shutdown -s');
    
else
    
    fprintf('please update channel after registration. \n');
    
end

if exist('recnumnum', 'var')
    fprintf([num2str(recnumnum) '/' num2str(recnumtotal) ' file in \n' masdata '\n' filename '\n']);
end

%% Write registered tiff file

for i = 1:length(imagelist)
    if i==1
        imwrite(imagelist{i,1}, [filename(1:end-4) '_reg.tif']);
    else
        imwrite(imagelist{i,1}, [filename(1:end-4) '_reg.tif'], 'writemode', 'append');
    end
end

fprintf('registered tiff file saved. \n');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%Plot sample images%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot ratiometric image

close all;
% page = 438;
% thresh = 250;
% sigma = 1;
% img_g = imgaussfilt(imagelist_g{page,1});
% img_r = imgaussfilt(imagelist_r{page,1});
% imagelist_logic = imagelist_g{page,1} > thresh;
% imagelist_ratio = ...
%     uint16(imagelist_logic) .* imagelist_g{page,1} ./ imagelist_r{page,1};
% imagelist_ratio_gauss = ...
%     imgaussfilt(uint16(imagelist_logic) .* img_g ./ img_r, sigma);
% % dilatedImage = imdilate(imagelist_ratio,strel('disk',10));
% % thinedImage = bwmorph(dilatedImage,'thin',10);
% 
% clr = plasma;
% clr(1,:) = [1 1 1];
% % subplot(131);
% imagesc(imagelist_ratio_gauss); colormap(clr); colorbar;
% axis equal; 
% %

% page = 438; 
% maxi = 10;
% sigma = 2; 
prompt = {'Page number:','Maximum of activity', ...
    'Rotation','Complete (1) or step-by-step (0) sample image'};
inputdlgtitle = 'Range and channel order';
dims = [1 35];
definput = {'1','6','0','1'}; % for muscle activity on Fig. 1: '438','10','2','east'
answer = inputdlg(prompt,inputdlgtitle,dims,definput);
page = str2double(answer{1,1});
maxi = str2double(answer{2,1});
rotangle = str2double(answer{3,1});
imagetype = answer{4,1};

col_num = size(imagelist_g{1,1},1);
row_num = size(imagelist_g{1,1},2);
xq_base = 1:row_num; 
yq_base = 1:col_num;
xq = kron(xq_base', ones(col_num,1));
yq = kron(ones(row_num,1), yq_base');
% xq_rep = repmat(xq_base, col_num, 1);
% xq = reshape(xq_rep, [], 1);
% yq_rep = repmat(yq_base, row_num, 1);
% yq = reshape(yq_rep', [], 1);
xv = [dorsal_data{2*page-1,1}; flipud(ventral_data{2*page-1,1})];
yv = [dorsal_data{2*page,1}; flipud(ventral_data{2*page,1})];
in = inpolygon(xq, yq, xv, yv);
% plot(xq(in), yq(in), 'or');
% 
xin = xq(in);
yin = yq(in);
img = double(imagelist_g{page,1})./double(imagelist_r{page,1});
% imgrot = imrotate(img, rotangle, 'bilinear', 'crop');
idx = sub2ind(size(img), yin, xin);
img_zeros = zeros(size(img));
img_zeros(idx) = img(idx);
% img_gauss = imgaussfilt(img_zeros, sigma);

clr = flipud(gray); 
clr(1,:) = [1 1 1];

f = figure;
switch imagetype
    case '1'
        imagesc(img_zeros);
    case '0'
        imagesc(imgrot);
end
% imagesc(img_zeros); % Choose this for a complete sample image
% imagesc(img); % Choose this for a step-by-step display of image analysis
colormap(clr); 
c = colorbar;
set(c, 'xcolor', 'k', 'ycolor', 'k', ...
    'ticklabels', '', 'units', 'normalized',...
    'direction', 'normal',...
    'location', 'north');
cpos = get(c, 'Position');
cpos(1) = 0.5; 
% cpos(2) = 0.6;
cpos(3) = 0.2*cpos(3);
cpos(4) = 0.2*cpos(4);
set(c, 'Position', cpos);
caxis([0 maxi]);
axis equal; 
set(gca, 'visible', 'off');

%% Adding Segmentation

lnwidth = 1.5;
hold on;
plot(dorsal_data{2*page-1,1}, dorsal_data{2*page,1}, ...
    'color', [1 1 1], 'linewidth', lnwidth);
plot(ventral_data{2*page-1,1}, ventral_data{2*page,1}, ...
    'color', [1 1 1], 'linewidth', lnwidth);

%% Adding skeletonization

plot(centerline_data_spline(:,2*page-1), centerline_data_spline(:,2*page), ...
    'color', [1 1 1], 'linewidth', lnwidth);

%% Adding body segments

xpt = [centerline_data_spline(:,2*page-1)'; dorsal_data{2*page-1,1}'];
ypt = [centerline_data_spline(:,2*page)'; dorsal_data{2*page,1}'];
plot(xpt, ypt, 'w', 'linewidth', lnwidth/2);
xpt = [centerline_data_spline(:,2*page-1)'; ventral_data{2*page-1,1}'];
ypt = [centerline_data_spline(:,2*page)'; ventral_data{2*page,1}'];
plot(xpt, ypt, 'w', 'linewidth', lnwidth/2);

%% Adding region of interest

recx = 76; recy = 150; recw = 20; rech = 9;
rectangle('position', [recx recy recw rech], 'edgecolor', 'k', 'linewidth', 2);

%% Save ratiometric figure

figname = [filename(1:end-4) '_' num2str(page)];
savefig(f, figname);
saveas(f, [figname '.emf'], 'meta');
fprintf('figure saved.\n');

%% Plot magnified ratiometric image

close all;
block = 74;
m = figure;
imagesc(img_zeros); colormap(clr); 
caxis([0 maxi]);
axis equal; 
set(gca, 'visible', 'off');
hold on;
plot(dorsal_data{2*page-1,1}, dorsal_data{2*page,1}, 'w', 'linewidth', lnwidth);
plot(centerline_data_spline(:,2*page-1), centerline_data_spline(:,2*page), 'w', 'linewidth', lnwidth);
plot(ventral_data{2*page-1,1}, ventral_data{2*page,1}, 'w', 'linewidth', lnwidth);

xpt = [centerline_data_spline(:,2*page-1)'; ...
    dorsal_data{2*page-1,1}'];
ypt = [centerline_data_spline(:,2*page)'; ...
    dorsal_data{2*page,1}'];
plot(xpt, ypt, 'w', 'linewidth', lnwidth);
h = fill([xpt(1,block) xpt(2,block) xpt(2,block-1) xpt(1,block-1)],...
    [ypt(1,block) ypt(2,block) ypt(2,block-1) ypt(1,block-1)], 'r', 'edgecolor', 'none');
set(h, 'facealpha', 0.4);

xpt = [centerline_data_spline(:,2*page-1)'; ...
    ventral_data{2*page-1,1}'];
ypt = [centerline_data_spline(:,2*page)'; ...
    ventral_data{2*page,1}'];
plot(xpt, ypt, 'w', 'linewidth', lnwidth);
h = fill([xpt(1,block) xpt(2,block) xpt(2,block-1) xpt(1,block-1)],...
    [ypt(1,block) ypt(2,block) ypt(2,block-1) ypt(1,block-1)], 'r', 'edgecolor', 'none');
set(h, 'facealpha', 0.4);

xlim([recx recx+recw]); ylim([recy recy+rech]);

%% Save magnified ratiometric image

savefig(m, [figname '_magnified']);
saveas(m, [figname '_magnified'], 'meta');
fprintf('figure saved.\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%Plot sample images%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%