%% Prcoess data and registrate channels

% Load tiff files
setup_proof_reading;
fprintf('tiff files loading finished. \n');
imshowpair(imagelist_g{1,1}, imagelist_r{1,1});
coef = 20;

% Registration
answer = questdlg('Need registration?', 'Image processing');
if strcmp(answer,'Yes')==1
    close all;
    figure; 
    subplot(1,2,1); imshowpair(imagelist_g{1,1}, imagelist_r{1,1});
    image_registration_tform; 
    subplot(1,2,2); imshowpair(movingRegistered{1,1}, imagelist_r{1,1});
end

%% Track neuron and generate signal

proof_reading(imagelist, [], filename, ...
    istart, iend, 1);

%% Save data

parts = strsplit(pathname, '\');
data_path = fullfile(parts{1,1:end-3}, 'Alpha_Data_Raw', parts{1,end-1});
warning('off'); mkdir(data_path); 
data_path_name = fullfile(data_path, [filename(1:end-4) '.mat']);
save(data_path_name, ...
    'signal', 'signal_mirror', 'ratio', 'neuron_position_data', 'dual_position_data');
fprintf('data saved. \n');

%% Write registered tiff file

for i = 1:length(imagelst)
    if i==1
        imwrite(imagelist{i,1}, [filename(1:end-4) '_reg.tif']);
    else
        imwrite(imagelist{i,1}, [filename(1:end-4) '_reg.tif'], 'writemode', 'append');
    end
end