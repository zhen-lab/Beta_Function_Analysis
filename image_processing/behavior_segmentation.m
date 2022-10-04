%% Image loading

setup_proof_reading;

%% Image segmentation

[dorsal_data, ventral_data, centerline_data, ...
    centerline_data_spline, curvdata, curvdatafiltered] = ...
    extract_centerline_behavior(imagelist, 0, ...
    4, 4, 1, 0.565, 0.565);

% Save data

parts = strsplit(pathname, '\');
data_path = fullfile(parts{1,1:end-3}, 'Alpha_Data_Raw', parts{1,end-1});
warning('off'); mkdir(data_path); 
data_path_name = fullfile(data_path, [filename(1:end-4) '.mat']);
save(data_path_name, ...
    'centerline_data_spline', 'curvdatafiltered', 'dorsal_data', 'ventral_data');
fprintf('data saved. \n');

% system('shutdown -s');