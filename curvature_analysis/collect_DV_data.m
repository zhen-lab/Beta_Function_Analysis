function [dorsal_ratio, ventral_ratio, dorsal, dorsal_gfp, ventral, ventral_gfp] = collect_DV_data(folder, frames, formats)

f = dir(fullfile(folder, '*.mat'));
numsamples = numel(f);
dorsal = NaN(frames, numsamples);
dorsal_gfp = NaN(frames, numsamples);
ventral = NaN(frames, numsamples);
ventral_gfp = NaN(frames, numsamples);

namespeciallist = uigetfile('.mat', 'Select contact-area-removed files', 'MultiSelect', 'on');
namestore = cell(1,1);
% headneck = 36;

for idx = 1 : numsamples
    
    namestore{1,1} = f(idx).name;
    ismb = cellfun(@(x) ismember(x, namespeciallist), namestore, 'UniformOutput', 0);
    pathname = fullfile(folder, f(idx).name);
    
    % Determine if the current file was analyzed with contact-area removal 
    if ismb{1,1}
        headneck = 51; % With contact-area removed
    else
        headneck = 36; % Without contact-area removed
    end
    
    if ~isempty(regexp(pathname, '\w*(?=concat)', 'match')) % In case of concatenated data
        
        load(pathname, 'total_dorsal', 'total_dorsal_r', 'total_ventral', 'total_ventral_r');
        frame_length = min(size(total_dorsal,2), frames); % in case there are recordings with shorter time window
        
        r_d = total_dorsal(headneck:end, 1:frame_length)./total_dorsal_r(headneck:end, 1:frame_length);
        finite_d = double(isfinite(r_d)); 
        finite_sum_d = sum(finite_d);
        r_d_sum = nansum(finite_d.*r_d);
        r_d_mean_finite = r_d_sum./finite_sum_d;
        
        r_v = total_ventral(headneck:end, 1:frame_length)./total_ventral_r(headneck:end, 1:frame_length);
        finite_v = double(isfinite(r_v)); 
        finite_sum_v = sum(finite_v);
        r_v_sum = nansum(finite_v.*r_v);
        r_v_mean_finite = r_v_sum./finite_sum_v;
        
        dorsal(1:frame_length, idx) = r_d_mean_finite'; % removes any inf numbers
        dorsal_gfp(1:frame_length, idx) = mean(total_dorsal(headneck:end, 1:frame_length))';
        
        ventral(1:frame_length, idx) = r_v_mean_finite'; % removes any inf numbers
        ventral_gfp(1:frame_length, idx) = mean(total_ventral(headneck:end, 1:frame_length))';
    
    else % In case of original data
        
        switch formats
            
            case 'New'
                
                load(pathname, 'dorsal_smd', 'dorsal_smd_r', 'ventral_smd', 'ventral_smd_r');
                frame_length = min(size(dorsal_smd,2), frames); % in case there are recordings with shorter time window
        
                r_d = dorsal_smd(headneck:end, 1:frame_length)./dorsal_smd_r(headneck:end, 1:frame_length);
                finite_d = double(isfinite(r_d)); 
                finite_sum_d = sum(finite_d);
                r_d_sum = nansum(finite_d.*r_d);
                r_d_mean_finite = r_d_sum./finite_sum_d;

                r_v = ventral_smd(headneck:end, 1:frame_length)./ventral_smd_r(headneck:end, 1:frame_length);        
                finite_v = double(isfinite(r_v)); 
                finite_sum_v = sum(finite_v);
                r_v_sum = nansum(finite_v.*r_v);
                r_v_mean_finite = r_v_sum./finite_sum_v;

                dorsal(1:frame_length, idx) = r_d_mean_finite'; % removes any inf numbers
                dorsal_gfp(1:frame_length, idx) = mean(dorsal_smd(headneck:end, 1:frame_length))';       

                ventral(1:frame_length, idx) = r_v_mean_finite'; % removes any inf numbers
                ventral_gfp(1:frame_length, idx) = mean(ventral_smd(headneck:end, 1:frame_length))';
            
            case 'Old'
                
                load(pathname, 'dorsal_smd', 'ventral_smd');
                frame_length = min(size(dorsal_smd,2), frames); % in case there are recordings with shorter time window
        
                r_d = dorsal_smd(headneck:end, 1:frame_length);
                finite_d = double(isfinite(r_d)); 
                finite_sum_d = sum(finite_d);
                r_d_sum = nansum(finite_d.*r_d);
                r_d_mean_finite = r_d_sum./finite_sum_d;

                r_v = ventral_smd(headneck:end, 1:frame_length);        
                finite_v = double(isfinite(r_v)); 
                finite_sum_v = sum(finite_v);
                r_v_sum = nansum(finite_v.*r_v);
                r_v_mean_finite = r_v_sum./finite_sum_v;

                dorsal(1:frame_length, idx) = r_d_mean_finite'; % removes any inf numbers
                dorsal_gfp(1:frame_length, idx) = mean(dorsal_smd(headneck:end, 1:frame_length))';       

                ventral(1:frame_length, idx) = r_v_mean_finite'; % removes any inf numbers
                ventral_gfp(1:frame_length, idx) = mean(ventral_smd(headneck:end, 1:frame_length))';
                
        end
        
    end
    
end

dorsal_start = repmat(dorsal(1, :), frames, 1);
ventral_start = repmat(ventral(1, :), frames, 1);
dorsal_ratio = (dorsal./dorsal_start)';
ventral_ratio = (ventral./ventral_start)';

dorsal = dorsal'; dorsal_gfp = dorsal_gfp';
ventral = ventral'; ventral_gfp = ventral_gfp';

end