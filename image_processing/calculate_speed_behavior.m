function [vel_anterior_sign, vel_posterior_sign, dorsal_data, ventral_data, centerline_data_spline] = calculate_speed_behavior(filename, pathname, data)

%%%%%%%%%%%%%%%%%%%%
% Set up constant variables
binning = 4;
        
% Obtain stage position data by selecting .tif file
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
    definput = {'4','75','2000','1','2000'};
    answer = inputdlg(prompt,dlgtitle,dims,definput);

    if isempty(answer)

        fprintf('User canceled. \n');

    else

        timeperframe = str2double(answer{2,1})/str2double(answer{3,1});
        magnification = str2double(answer{1,1});
        istart = str2double(answer{4,1});
        iend = str2double(answer{5,1});
%                 min_consecutive = str2double(answer{6,1});

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Determine the scale by magnification
%                 switch magnification  
%                     case 10
%                         % 10x obj um/pix with  4x binning
%                         scale = 0.63 * binning;        
%                     case 16    
%                         % 16x obj um/pix with  4x binning
%                         scale = 0.39375 * binning;        
%                     case 20
%                         % 20X obj um/pix with  4x binning
%                         scale = 0.315 * binning;        
%                     case 40
%                         % 40x obj um/pix with  4x binning
%                         scale = 0.1575 * binning;       
%                     case 63
%                         % 63x obj um/pix with  4x binning (0.1 um/pix?)
%                         scale = 0.1 * binning;     
%                     otherwise
%                         warning('No such magnification available. Only type in 10x, 16x, 20x, 40x, or 63x. \n');
%                 end
        scale = 6.3/magnification*binning;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$$$
        % Read stage positions
        xyz = read_motorstage_position_behavior(data);
        % xyz_col is a cell array
        % Not sure if cellfun is faster here
        pos_stage = zeros(size(xyz, 1), 2);
        for j = 1:size(xyz,1)
            for i = 1:2
                pos_stage(j,i) = str2double(xyz{j,i});
            end
        end
        pos_stage_cropped = pos_stage(istart:iend, :);
        fprintf('Stage position loaded. \n');        

        [position_file, ~] = uigetfile('*.mat', ['Select position data for ' filename]);

        if str2double(position_file)==0

            fprintf('User canceled. \n');

        else

            load(position_file, 'centerline_data_spline', 'dorsal_data', 'ventral_data');
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             % Revert the position data if the segmentation used smaller
%             % area (1/4 to 3/4 of the size) for analysis
%             if smallarea == 1
%                 imgsize = fliplr(size(data{1,1}{1,1})); % size 1 is for y axis, size 2 is for x axis
%                 clsize = size(centerline_data_spline);
%                 displacement = kron(ones(clsize(1), clsize(2)/2), imgsize/4);
%                 centerline_data_spline  = centerline_data_spline + displacement;
%                 dorsal_data(1:2:end) = cellfun(@(x) x+imgsize(1)/4, dorsal_data(1:2:end), 'UniformOutput', false);
%                 dorsal_data(2:2:end) = cellfun(@(x) x+imgsize(2)/4, dorsal_data(2:2:end), 'UniformOutput', false);
%                 ventral_data(1:2:end) = cellfun(@(x) x+imgsize(1)/4, ventral_data(1:2:end), 'UniformOutput', false);
%                 ventral_data(2:2:end) = cellfun(@(x) x+imgsize(2)/4, ventral_data(2:2:end), 'UniformOutput', false);
%             end
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Head position
            pos_anterior_all = [centerline_data_spline(1, 1:2:end)' centerline_data_spline(1, 2:2:end)'];
            pos_anterior = pos_anterior_all(istart:iend, :);
            % Center position
            pos_posterior_all = [centerline_data_spline(end/2, 1:2:end)' centerline_data_spline(end/2, 2:2:end)'];
            pos_posterior = pos_posterior_all(istart:iend, :);
            fprintf('Head and center positions loaded. \n');
           
            
            % Smooth neurons positions, and substract stage positions by a scale
%                     pos_anterior_smd = smoothdata(pos_anterior, 'movmedian', windowlength_before);
%                     pos_posterior_smd = smoothdata(pos_posterior, 'movmedian', windowlength_before);
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
            vel_anterior_abs = sqrt(sum(vec_aa.^2, 2))/timeperframe;
            vel_posterior_abs = sqrt(sum(vec_pp.^2, 2))/timeperframe;
            % Signed velocity
            vel_anterior_sign = sign_aa_ap .* vel_anterior_abs;
            vel_posterior_sign = sign_pp_ap .* vel_posterior_abs;
            % Smoothened signed velocity
%                     vel_anterior_sign_smd = smoothdata(vel_anterior_sign, 'movmedian', windowlength_after);
%                     vel_posterior_sign_smd = smoothdata(vel_posterior_sign, 'movmedian', windowlength_after);
            vel_ap_sign_mean = mean([vel_anterior_sign vel_posterior_sign], 2);

            s = figure;
            hold on; title('Velocity for behavior wout smoothening')
            plot(vel_anterior_sign, 'r');
            plot(vel_posterior_sign, 'b');
            plot(vel_ap_sign_mean, 'm');
            plot(1:length(vel_anterior_sign), zeros(length(vel_anterior_sign),1), 'k');
            hold off;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Save data
            parts = strsplit(pathname, '\');
            data_path = fullfile(parts{1,1:end-3}, 'Alpha_Data_Raw', parts{1,end-1});
            warning('off'); mkdir(data_path); 
            data_path_name = fullfile(data_path, [filename(1:end-4) '_velocity.mat']);
            save(data_path_name, ...
                'vel_anterior_sign', 'vel_posterior_sign',...
                'vel_ap_sign_mean');
            fprintf('data saved. \n');

            % Save figures
            data_path = fullfile(parts{1,1:end-3}, 'Alpha_Data_Plot', parts{1,end-1});
            warning('off'); mkdir(data_path); 
            savefig(s, [data_path '\' filename(1:end-4) '_velocity.fig']);
            fprintf('figures saved. \n');
            fprintf([filename ' analysis finished. \n']);

        end

    end

end

end
