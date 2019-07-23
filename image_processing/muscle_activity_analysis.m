%% Pre-prcoess data

% Load tiff files
setup_proof_reading;

if ~(filename==0)
    
    fprintf('tiff files loading finished. \n');

    % Register RFP and GFP channels
    figure;

    subplot(1,2,1); imshowpair(imagelist_g{1,1}, imagelist_r{1,1});
    %
    coef = 20;
    image_registration_tform;
    subplot(1,2,2); imshowpair(movingRegistered{frmnum,1}, imagelist_r{frmnum,1});
    % pause(2);
    imagelist_g = movingRegistered;
    imagelist_use = imagelist_g;

end

% Remove edges from registration

for i = 1:length(imagelist)
    img = imagelist_g{i,1};
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

%% Process data

%imagelist_g = movingRegistered;
% Delineate dorsal and ventral muscles
tic; close all;
[dorsal_data, ventral_data, centerline_data, centerline_data_spline, curvdata, curvdatafiltered] = ...
    extract_centerline_vd(...
    imagelist_use, 6, 5, 2.5);
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
% data_path_new = fullfile(data_path, 'Alpha_Data_Raw', 'Muscle_Interneurons_Ablated');
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