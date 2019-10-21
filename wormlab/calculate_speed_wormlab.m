function [vel_anterior_sign, vel_posterior_sign, vel_ap_sign_smd_mean] = calculate_speed_wormlab()

%%%%%%%%%%%%%%%%%%%%
% Set up constant variables
windowlength_after = 9;
binning = 4;
ysize = 256;
    
% Obtain stage position data by selecting .tif file
[filename, pathname]  = uigetfile('*.tif', 'Select One TIFF file');

if str2double(filename)==0

    fprintf('User canceled. \n');

else

    % Use prompt to define variables
    prompt = {'Magnification',...
        'Total time','Total frames',...
        'Start frame','End frame'};
%                 'Minimum consecutive frames'};
    dlgtitle = ['Input for speed calculation for ' filename];
    dims = [1 35];
    definput = {'10','75','2000','1','5000'};
    answer = inputdlg(prompt,dlgtitle,dims,definput);

    if isempty(answer)

        fprintf('User canceled. \n');

    else

        timeperframe = str2double(answer{2,1})/str2double(answer{3,1});
        magnification = str2double(answer{1,1});
        istart = str2double(answer{4,1});
        iend = str2double(answer{5,1});
        scale = 6.3/magnification*binning;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Read stage positions
        pos_stage = get_stage_position;
        pos_stage_size = size(pos_stage,1);
        pos_stage_cropped = pos_stage(istart:iend, 1:2);
        fprintf('Stage position loaded. \n');  
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Read head and midpoint positions
        position_file = uigetfile('.csv', 'Select head file');

        if str2double(position_file)==0

            fprintf('User canceled. \n');

        else

            T = readtable(position_file, 'HeaderLines', 5);
            pos_anterior_all = table2array(T);

            position_file = uigetfile('.csv', 'Select midpoint file');

            if str2double(position_file)==0

                fprintf('User canceled. \n');

            else

                T = readtable(position_file, 'HeaderLines', 5);
                pos_posterior_all = table2array(T);    
                headmidcolstart = 3;
                headmidcolend = 4;
                coef = 5;
                
                % Create vector for the intact frames
                frm_all = 1:pos_stage_size;
                
                % Anterior position
                % Pad missing frames with NaN         
                stage_obj_overlap_anterior = intersect(frm_all, pos_anterior_all(:,1));
                pos_anterior_pos_padded = nan(length(frm_all), 2);
                pos_anterior_pos_padded(stage_obj_overlap_anterior, :) = pos_anterior_all(:,headmidcolstart:headmidcolend);
                % Crop by user defined frames
                pos_anterior_bfcvs = pos_anterior_pos_padded(istart:iend, :);
                % Convert to pixels (WormLab magnifies units by a factor of 5, and reverses y axis)
                pos_anterior = [pos_anterior_bfcvs(:,1)/coef ysize-pos_anterior_bfcvs(:,2)/coef];
                
                % Posterior position
                % Pad missing frames with NaN 
                stage_obj_overlap_posterior= intersect(frm_all, pos_posterior_all(:,1));
                pos_posterior_pos_padded = nan(length(frm_all), 2);
                pos_posterior_pos_padded(stage_obj_overlap_posterior, :) = pos_posterior_all(:,headmidcolstart:headmidcolend);
                % Crop by user defined frames
                pos_posterior_bfcvs = pos_posterior_pos_padded(istart:iend, :);
                % Convert to pixels (WormLab magnifies units by a factor of 5, and reverses y axis)
                pos_posterior = [pos_posterior_bfcvs(:,1)/coef ysize-pos_posterior_bfcvs(:,2)/coef];
                
                fprintf('Head and center positions padded and loaded. \n');
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Smooth neurons positions, and substract stage positions by a scale
                pos_anterior_adj = pos_anterior * scale - pos_stage_cropped;
                pos_posterior_adj = pos_posterior * scale - pos_stage_cropped;

                % Reference vetor AP
                vec_ap = pos_anterior - pos_posterior; 
                vec_ap = vec_ap(1:end-1, :); % To match sizes with velocity
                % Velocity vectors A'A and P'P
                vec_aa = diff(pos_anterior_adj);
                vec_pp = diff(pos_posterior_adj);
                % Dot product between A'A and AP, P'P and AP to determine direction
                sign_aa_ap = sign(dot(vec_aa, vec_ap, 2));
                sign_pp_ap = sign(dot(vec_pp, vec_ap, 2));
                % Absolute values for velocity
                vel_anterior_abs = sqrt(nansum(vec_aa.^2, 2))/timeperframe;
                vel_posterior_abs = sqrt(nansum(vec_pp.^2, 2))/timeperframe;
                % Signed velocity
                vel_anterior_sign = sign_aa_ap .* vel_anterior_abs;
                vel_posterior_sign = sign_pp_ap .* vel_posterior_abs;
                % Smoothened signed velocity
                vel_ap_sign_mean = mean([vel_anterior_sign vel_posterior_sign], 2);
                % Smoothened signed velocity
                vel_anterior_sign_smd = smoothdata(vel_anterior_sign, 'movmedian', windowlength_after);
                vel_posterior_sign_smd = smoothdata(vel_posterior_sign, 'movmedian', windowlength_after);
                vel_ap_sign_smd_mean = smoothdata(mean([vel_anterior_sign_smd vel_posterior_sign_smd], 2), 'movmedian', windowlength_after);

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Plot figures
                figure(1);
                hold on; title('Velocity for behavior wout smoothening');
                plot(vel_anterior_sign, 'r');
                plot(vel_posterior_sign, 'b');
                plot(vel_ap_sign_mean, 'm');
                plot([1 length(vel_anterior_sign)], [0 0], 'k');
                hold off;

                s = figure(2);
                hold on; title('Velocity for behavior with smoothening');                  
                plot(vel_ap_sign_smd_mean, 'k');
                plot([1 length(vel_anterior_sign)], [0 0], 'k');
                hold off;

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% Save data
                parts = strsplit(pathname, '/');
                data_path = fullfile(parts{1,1:end-1}, 'Alpha_Data_Raw', parts{1,end});
                warning('off'); mkdir(data_path); 
                data_path_name = fullfile(data_path, [filename(1:end-4) '_velocity.mat']);
                save(data_path_name, ...
                    'vel_anterior_sign', 'vel_posterior_sign',...
                    'vel_ap_sign_mean', 'vel_ap_sign_smd_mean');
                fprintf('data saved. \n');

                % Save figures
                data_path = fullfile(parts{1,1:end-1}, 'Alpha_Data_Plot', parts{1,end});
                warning('off'); mkdir(data_path); 
                savefig(s, [data_path '\' filename(1:end-4) '_velocity.fig']);
                fprintf('figures saved. \n');
                fprintf([filename ' analysis finished. \n']);

            end

        end

    end

end
    
end

% N.B to smooth graphs, try: hampel(vel_posterior_sign, 13);
