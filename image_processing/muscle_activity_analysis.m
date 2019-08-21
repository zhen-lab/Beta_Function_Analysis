%% Load tiff files

setup_proof_reading;

%% Check which channel needs to move for analyzed data

close all;
figure;
subplot(121);
imshow(imagelist_r{1}); caxis('auto');
hold on
plot(dorsal_data{2*1-1,1}, dorsal_data{2*1,1}, 'r');
plot(ventral_data{2*1-1,1}, ventral_data{2*1,1}, 'b');
title('RFP/Left');
subplot(122);
imshow(imagelist_g{1}); caxis('auto');
hold on
plot(dorsal_data{2*1-1,1}, dorsal_data{2*1,1}, 'r');
plot(ventral_data{2*1-1,1}, ventral_data{2*1,1}, 'b');
title('GCaMP/Right');

%% Registration of GFP and RFP channels

if ~(filename==0)
    
    fprintf('tiff files loading finished. \n');
    channeltomove = questdlg('Move left or right channel?',...
        'Channel to move', 'Left', 'Right', 'Right');
    
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

%% Remove edges from registration

imagelist_g = movingRegistered;
imagelist_use = imagelist_g;

for i = 1:length(imagelist)
    img = movingRegistered{i};
    img = img + mean(img,[1 2])*uint16(img==0);
    imagelist_use{i,1} = img;
end
fprintf('edges removed. \n');

% Overlay two channels if necessary

imagelist_gr = imagelist_use;

for i = 1:length(imagelist)
    imagelist_gr{i,1} = imagelist_use{i,1}+imagelist_r{i,1};
end
imagelist_use = imagelist_gr;
fprintf('channels overlaid. \n');

%% Remove contact area manually if necessary

imagelist_g = movingRegistered;

% Select channel
channel = questdlg('Use GFP or RFP channel?', 'Channel', 'GFP', 'RFP', 'Both', 'RFP');
switch channel
    case 'GFP'
        imagelist_use = imagelist_g;
    case 'RFP'
        imagelist_use = imagelist_r;
    case 'Both'
        imagelist_gr = imagelist_g;
        for i = 1:size(imagelist_gr, 1)
            imagelist_gr{i,1} = imagelist_g{i,1}+imagelist_r{i,1};
        end
        imagelist_use = imagelist_gr;
end

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
    pad_noise = imnoise(pad, 'gaussian');
    imagelist_updated{i,1} = uint16(pad.*pad_noise.*double(logic)) + uint16(~(logic)).*img;
end
subplot(122); imshow(imagelist_updated{1,1}); 
title('After removing contact area');
caxis([0 mm]);

% Update images to use later
imagelist_use = imagelist_updated;
fprintf('images are updated. \n');

%% Process data

%imagelist_g = movingRegistered;
% Delineate dorsal and ventral muscles
tic; close all;
[dorsal_data, ventral_data, centerline_data, centerline_data_spline, curvdata, curvdatafiltered] = ...
    extract_centerline_vd(...
    imagelist_use, 6, 4, 2.5);
toc;

% Calculate dorsal and ventral muscle activities, save output

% Generate data
tic;  fprintf('activity analysis kicks off \n');
[dorsal_smd, ventral_smd, dorsal_smd_r, ventral_smd_r] = ...
    activity_all(imagelist_g, imagelist_r, range, dorsal_data, ventral_data, centerline_data_spline, curvdatafiltered);
fprintf('activity analysis finished. \n');
figure;
subplot(1,2,1); imagesc(dorsal_smd./dorsal_smd_r); title('Dorsal');
subplot(1,2,2); imagesc(ventral_smd./ventral_smd_r); title('Ventral');
toc;

% Save data
parts = strsplit(pathname, '\');
data_path = fullfile(parts{1,1:end-3}, 'Alpha_Data_Raw', parts{1,end-1});
warning('off'); mkdir(data_path); 
data_path_name = fullfile(data_path, [filename(1:end-4) '.mat']);
save(data_path_name, ...
    'centerline_data_spline', 'curvdatafiltered', 'dorsal_data', 'dorsal_smd', 'dorsal_smd_r', 'ventral_data', 'ventral_smd', 'ventral_smd_r');
fprintf('data saved. \n');

% system('shutdown -s');

%% Write registered tiff file

for i = 1:length(imagelist)
    if i==1
        imwrite(imagelist{i,1}, [filename(1:end-4) '_reg.tif']);
    else
        imwrite(imagelist{i,1}, [filename(1:end-4) '_reg.tif'], 'writemode', 'append');
    end
end
fprintf('registered tiff file saved. \n');