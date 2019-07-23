%% Pre-prcoess data

% Select tiff files
[totalfilename, totalpathname]  = uigetfile({'*.tif'}, 'Select tiff files', 'MultiSelect', 'on');
numfiles = size(totalfilename, 2);

% Prepare cell arrays for storage
mst_imagelist_g = cell(numfiles, 1);
mst_imagelist_r = cell(numfiles, 1);
mst_range = cell(numfiles, 1);
mst_dorsal_data = cell(numfiles, 1); 
mst_ventral_data = cell(numfiles, 1); 
mst_centerline_data_spline = cell(numfiles, 1); 
mst_curvdatafiltered = cell(numfiles, 1);

%% Registration and image processing
for l = 1:numfiles
    
    % Load tiff files
    [imagelist, imagelist_g, imagelist_r, numframes, num, range] = setup_proof_reading_automated_adjusted(totalpathname, totalfilename{1,l}); % 'range' variable comes from here
    fprintf('tiff files loading finished. \n');

    % Register RFP and GFP channels
    coef = 5; image_registration_tform;
    figure; 
    subplot(1,2,1); imshowpair(imagelist_g{frmnum,1}, imagelist_r{frmnum,1});
    subplot(1,2,2); imshowpair(movingRegistered{frmnum,1}, imagelist_r{frmnum,1});
    imagelist_g = movingRegistered;
    fprintf('channel registration finished. \n');
    pause;
    close all;
    
    % Delineate dorsal and ventral muscles
%     tic; 
    [dorsal_data, ventral_data, centerline_data, centerline_data_spline, curvdata, curvdatafiltered] = ...
        extract_centerline_vd_automated_adjusted(...
        imagelist_g, totalfilename{1,l}, 5, 5);
%     toc;

    % Store information for individual movie
    mst_imagelist_g{l,1} = imagelist_g;
    mst_imagelist_r{l,1} = imagelist_r;
    mst_range{l,1} = range;
    mst_dorsal_data{l,1} = dorsal_data; 
    mst_ventral_data{l,1} = ventral_data; 
    mst_centerline_data_spline{l,1} = centerline_data_spline; 
    mst_curvdatafiltered{l,1} = curvdatafiltered;

end

% Save information for the next step
parts = strsplit(totalpathname, '\');
data_path = fullfile(parts{1,1:end-3}, 'Alpha_Data_Raw', parts{1,end-1});
warning('off'); mkdir(data_path); 
data_path_name = fullfile(data_path, ['Image_processing_' totalfilename{1,1} '-' totalfilename{1,end} '.mat']);
save(data_path_name, ...
    'mst_imagelist_g', 'mst_imagelist_r', 'mst_range', 'mst_dorsal_data', 'mst_ventral_data', 'mst_centerline_data_spline', 'mst_curvdatafiltered');
fprintf('information data saved. \n');

%% Calculate dorsal and ventral muscle activities, save output

% Select computer behavior
answer = questdlg('Turn off computer after analysis?', 'Select computer behavior', 'Yes', 'No', 'No');


% Loop for analysis
for j = 1:numfiles
    % Generate data
%     tic;  
    fprintf('activity analysis kicks off \n');
    [dorsal_smd, ventral_smd, dorsal_smd_r, ventral_smd_r] = ...
        activity_all(mst_imagelist_g{j,1}, mst_imagelist_r{j,1}, mst_range{j,1}, mst_dorsal_data{j,1}, mst_ventral_data{j,1}, mst_centerline_data_spline{j,1}, mst_curvdatafiltered{j,1});
    fprintf('activity analysis finished. \n');
%     figure;
%     subplot(1,2,1); imagesc(dorsal_smd./dorsal_smd_r); title('Dorsal');
%     subplot(1,2,2); imagesc(ventral_smd./ventral_smd_r); title('Ventral');
%     toc;

    % Save data
    parts = strsplit(totalpathname, '\');
    data_path = fullfile(parts{1,1:end-3}, 'Alpha_Data_Raw', parts{1,end-1});
    warning('off'); mkdir(data_path); 
    data_path_name = fullfile(data_path, [totalfilename{1,j}(1:end-4) '.mat']);
    save(data_path_name, ...
        'centerline_data_spline', 'curvdatafiltered', 'dorsal_data', 'dorsal_smd', 'dorsal_smd_r', 'ventral_data', 'ventral_smd', 'ventral_smd_r');
    % data_path_new = fullfile(data_path, 'Alpha_Data_Raw', 'Muscle_Interneurons_Ablated');
    fprintf([num2str(j) '/'  num2str(numfiles) ' activity data saved. \n']);

end

fprintf('all data saved. \n');

% Determine computer behavior
if strcmp(answer, 'Yes')==1
    system('shutdown -s');
end

%% Write registered tiff file

for l = 1:length(imagelist)
    if l==1
        imwrite(imagelist{l,1}, [filename(1:end-4) '_reg.tif']);
    else
        imwrite(imagelist{l,1}, [filename(1:end-4) '_reg.tif'], 'writemode', 'append');
    end
end