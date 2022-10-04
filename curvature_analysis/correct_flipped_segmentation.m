%% Clear everything

clear; close all;

%% Find the start of incorrect frames
% Navigate to the folder containing analyzed data
% Load the data and check approximately where the flip occurs

imagesc(curvdatafiltered);

%% Estimate the frame immediately before the flip

frame_correct_est = 300;
nn = 15;
close all; 
figure;
for i = 1:nn
    
    subplot(1,nn,i); title(num2str(frame_correct_est+i));
%     imagesc(imagelist_g{frame_correct_est+i,1});
    xd = dorsal_data{(frame_correct_est+i)*2-1,1};
    yd = dorsal_data{(frame_correct_est+i)*2,1};
    xv = ventral_data{(frame_correct_est+i)*2-1,1};
    yv = ventral_data{(frame_correct_est+i)*2,1};
    hold on;
    plot(xd, yd, 'r'); plot(xv, yv, 'b');
    plot(xd(1), yd(1), 'om');
    set(gca, 'ydir', 'reverse', 'xticklabel', [], 'yticklabel', [], 'ticklength', [0 0]);

end

%% Correct segmentations

frame_correct_ini = 312;
frame_correct_fin = 485;

dorsal_data_new = dorsal_data;
ventral_data_new = ventral_data;
dorsal_data_new(frame_correct_ini*2-1:frame_correct_fin*2-1) = ventral_data(frame_correct_ini*2-1:frame_correct_fin*2-1);
ventral_data_new(frame_correct_ini*2-1:frame_correct_fin*2-1) = dorsal_data(frame_correct_ini*2-1:frame_correct_fin*2-1);
dorsal_data_new_flip = dorsal_data_new;
ventral_data_new_flip = ventral_data_new;
dorsal_data_new_flip(frame_correct_ini*2-1:frame_correct_fin*2-1) = cellfun(@flipud, dorsal_data_new(frame_correct_ini*2-1:frame_correct_fin*2-1), 'UniformOutput', false);
ventral_data_new_flip(frame_correct_ini*2-1:frame_correct_fin*2-1) = cellfun(@flipud, ventral_data_new(frame_correct_ini*2-1:frame_correct_fin*2-1), 'UniformOutput', false);

centerline_data_spline_flip = centerline_data_spline;
centerline_data_spline_flip(:, frame_correct_ini*2-1:frame_correct_fin*2-1) = flipud(centerline_data_spline(:, frame_correct_ini*2-1:frame_correct_fin*2-1));

curvdatafiltered_flip = curvdatafiltered;
curvdatafiltered_flip(:, frame_correct_ini:frame_correct_fin) = -flipud(curvdatafiltered(:, frame_correct_ini:frame_correct_fin));

nn = 20;
figure; 
for i = 1:nn
    
    subplot(1,nn,i); hold on;
    title(num2str(frame_correct_est+i)); 
    xd = dorsal_data_new_flip{(frame_correct_est+i)*2-1,1};
    yd = dorsal_data_new_flip{(frame_correct_est+i)*2,1};
    xv = ventral_data_new_flip{(frame_correct_est+i)*2-1,1};
    yv = ventral_data_new_flip{(frame_correct_est+i)*2,1};
    xc = centerline_data_spline_flip(:,(frame_correct_est+i)*2-1);
    yc = centerline_data_spline_flip(:,(frame_correct_est+i)*2);
    plot(xd, yd, 'r'); plot(xv, yv, 'b'); plot(xc, yc, 'y');
    plot(xd(1), yd(1), 'om'); plot(xc(1),yc(1), 'oy', 'markersize', 10);
    set(gca, 'ydir', 'reverse', 'xticklabel', [], 'yticklabel', [], 'ticklength', [0 0]);

end

%% Import the corresponding recording for correction
% Nagivate to the folder containing recordings
setup_proof_reading_split_alternate;

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

% if ~(filename==0)
%     
%     fprintf('tiff files loading finished. \n');
% 
%     % Register RFP and GFP channels
%     figure;
%     subplot(1,2,1); imshowpair(imagelist_g{1,1}, imagelist_r{1,1});
%     image_registration_tform;
%     subplot(1,2,2); imshowpair(movingRegistered{frmnum,1}, imagelist_r{frmnum,1});
%     
% end

%% Caculate activity correctly
% Stay in the folder containing recordings

imagelist_g = movingRegistered;
range = [1 size(imagelist_g,1)];
[dorsal_smd, ventral_smd, dorsal_smd_r, ventral_smd_r] ...
    = activity_all(imagelist_g, imagelist_r, range, dorsal_data_new_flip, ventral_data_new_flip, centerline_data_spline_flip, curvdatafiltered_flip);

% Update data
dorsal_data = dorsal_data_new_flip;
ventral_data = ventral_data_new_flip;
centerline_data_spline = centerline_data_spline_flip;
curvdatafiltered = curvdatafiltered_flip;

% Save data
parts = strsplit(pathname, '\');
data_path = fullfile(parts{1,1:end-3}, 'Alpha_Data_Raw', parts{1,end-1});
warning('off'); mkdir(data_path); 
data_path_name = fullfile(data_path, [filename(1:end-4) '.mat']);
% save(data_path_name, ...
%     'centerline_data_spline', 'curvdatafiltered', 'dorsal_data', 'dorsal_smd', 'dorsal_smd_r', 'ventral_data', 'ventral_smd', 'ventral_smd_r');
save(data_path_name, ...
    'centerline_data_spline', 'curvdatafiltered', ...
    'dorsal_data', 'dorsal_smd', 'dorsal_smd_r', 'dorsal_raw', 'dorsal_raw_r', ...
    'ventral_data', 'ventral_smd', 'ventral_smd_r', 'ventral_raw', 'ventral_raw_r', ...
    'tform', 'channeltomove');
% data_path_new = fullfile(data_path, 'Alpha_Data_Raw', 'Muscle_Interneurons_Ablated');
fprintf('data saved. \n');
