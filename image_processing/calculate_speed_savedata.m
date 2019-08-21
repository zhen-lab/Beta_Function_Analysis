function [] = calculate_speed_savedata(magnification, totaltime, totalframes)

% Caculate velocity
[filename, pathname, vel_anterior_sign, vel_posterior_sign, ...
    vel_anterior_sign_smd, vel_posterior_sign_smd, vel_ap_sign_smd_mean] = calculate_speed(magnification, totaltime, totalframes);

% Save data
parts = strsplit(pathname, '\');
data_path = fullfile(parts{1,1:end-3}, 'Alpha_Data_Raw', parts{1,end-1});
warning('off'); mkdir(data_path); 
data_path_name = fullfile(data_path, [filename(1:end-4) '_velocity.mat']);
save(data_path_name, ...
    'vel_anterior_sign', 'vel_posterior_sign', 'vel_anterior_sign_smd', 'vel_posterior_sign_smd', 'vel_ap_sign_smd_mean');
% data_path_new = fullfile(data_path, 'Alpha_Data_Raw', 'Muscle_Interneurons_Ablated');
fprintf('data saved. \n');

end