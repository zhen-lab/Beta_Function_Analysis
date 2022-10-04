function [signal_directions_separated_forward, signal_directions_separated_backward, vel_anterior_sign, vel_posterior_sign] = calculate_speed()

%%%%%%%%%%%%%%%%%%%%
% Set up constant variables
windowlength_before = 3;
windowlength_after = 9;
binning = 4;

% Obtain stage position data by selecting .tif file
[filename, pathname]  = uigetfile('*.tif', 'Select One TIFF file');

if str2double(filename)==0
    
    fprintf('User canceled. \n');
    
else
    
    % Use prompt to define variables
    prompt = {'Magnification','Total time','Total frames','Start frame','End frame','Minimum consecutive frames'};
    dlgtitle = ['Input for speed calculation for ' filename];
    dims = [1 35];
    definput = {'63','180','1800','1','1800','2'};
    answer = inputdlg(prompt,dlgtitle,dims,definput);

    if isempty(answer)

        fprintf('User canceled. \n');

    else

        timeperframe = str2double(answer{2,1})/str2double(answer{3,1});
        magnification = str2double(answer{1,1});
        istart = str2double(answer{4,1});
        iend = str2double(answer{5,1});
        min_consecutive = str2double(answer{6,1});

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Determine the scale by magnification
        switch magnification  
            case 10
                % 10x obj um/pix with  4x binning
                scale = 0.63 * binning;        
            case 16    
                % 16x obj um/pix with  4x binning
                scale = 0.39375 * binning;        
            case 20
                % 20X obj um/pix with  4x binning
                scale = 0.315 * binning;        
            case 40
                % 40x obj um/pix with  4x binning
                scale = 0.1575 * binning;       
            case 63
                % 63x obj um/pix with  4x binning (0.1 um/pix?)
                scale = 0.1 * binning;     
            otherwise
                warning('No such magnification available. Only type in 10x, 16x, 20x, 40x, or 63x. \n');
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%$$$
        % Read stage positions
        xyz = read_motorstage_position(filename);
        % xyz_col is a cell array
        % Not sure if there is a better way to convert????
        pos_stage = zeros(size(xyz, 1), 2);
        for j = 1:size(xyz,1)
            for i = 1:2
                pos_stage(j,i) = str2double(xyz{j,i});
            end
        end
        pos_stage_cropped = pos_stage(istart:iend, :);
        fprintf('Stage position loaded. \n');

        % Obtain neurons positions data by selecting two .mat files
        [position_file, ~] = uigetfile('*.mat', ['Select anterior and posterior neurons for ' filename], 'MultiSelect', 'on');
        
        if str2double(position_file)==0
            
            fprintf('User canceled. \n');
            
        else
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Caculate velocity
            if size(position_file, 2) == 2

                % Read neurons positions
                load(position_file{1});
                pos_anterior = neuron_position_data;
                signal_anterior = signal{1,1}./signal_mirror{1,1};
                signal_anterior = signal_anterior(2:end); % To match the number of data points in velocity
                ymax_1 = max(signal_anterior);
                load(position_file{2});
                pos_posterior = neuron_position_data;
                signal_posterior = signal{1,1}./signal_mirror{1,1};
                signal_posterior = signal_posterior(2:end); % To match the number of data points in velocity
                ymax_2 = max(signal_posterior);
                ymax = max([ymax_1 ymax_2]);
                fprintf('Neurons positions loaded. \n');

                % Smooth neurons positions, and substract stage positions by a scale
                pos_anterior_smd = smoothdata(pos_anterior, 'movmedian', windowlength_before);
                pos_posterior_smd = smoothdata(pos_posterior, 'movmedian', windowlength_before);
                pos_anterior_adj = pos_anterior_smd * scale - pos_stage_cropped;
                pos_posterior_adj = pos_posterior_smd * scale - pos_stage_cropped;

                % Reference vetor AP
                vec_ap = pos_anterior_smd - pos_posterior_smd; 
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
                vel_anterior_sign_smd = smoothdata(vel_anterior_sign, 'movmedian', windowlength_after);
                vel_posterior_sign_smd = smoothdata(vel_posterior_sign, 'movmedian', windowlength_after);
                vel_ap_sign_smd_mean = smoothdata(mean([vel_anterior_sign_smd vel_posterior_sign_smd], 2), 'movmedian', windowlength_after);

                figure;
                subplot(121); hold on; title('Before smoothening')
                    plot(vel_anterior_sign, 'r');
                    plot(vel_posterior_sign, 'b');
                    plot(1:length(vel_anterior_sign), zeros(length(vel_anterior_sign),1), 'k');
                    yl = ylim;
                    hold off;
                subplot(122); hold on; title ('After smoothening')
                    plot(vel_anterior_sign_smd, 'r');
                    plot(vel_posterior_sign_smd, 'b');
                    plot(vel_ap_sign_smd_mean, 'm');
                    plot(1:length(vel_anterior_sign), zeros(length(vel_anterior_sign),1), 'k');
                    ylim(yl);
                    hold off;

            else
                warning('Please select two neurons. \n');

            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Obtain blocks of forward and backward movements
            forward = double(vel_ap_sign_smd_mean > 0); 
            forward_buffered = [0; double(forward); 0]; % To enclose pulses on both ends
            backward = double(vel_ap_sign_smd_mean < 0); 
            backward_buffered = [0; double(backward); 0]; % To enclose pulses on both ends
            [W_forward,INITCROSS_forward,FINALCROSS_forward] = pulsewidth(double(forward_buffered));
            [W_backward,INITCROSS_backward,FINALCROSS_backward] = pulsewidth(double(backward_buffered));

            % Remove blocks with width no greater than the threshold for forward
            forward_pruned = forward;
            forward_pruned((ceil(INITCROSS_forward(W_forward<=min_consecutive))-1)...
                :(floor(FINALCROSS_forward(W_forward<=min_consecutive))-1)) = 0;
            INITCROSS_forward(W_forward<=min_consecutive) = [];
            FINALCROSS_forward(W_forward<=min_consecutive) = [];
            W_forward(W_forward<=min_consecutive) = [];
            if isequal(size(W_forward), size(INITCROSS_forward), size(FINALCROSS_forward))
                fprintf(['Forward Events #: ' num2str(size(W_forward, 1)) '\n']);
            else
                fprintf(['Forward Unequal: ' num2str(size(W_forward)) ' ' num2str(size(INITCROSS_forward)) ' ' num2str(size(FINALCROSS_forward)) '\n']);
            end
            % Remove blocks with width no greater than the threshold for backward
            backward_pruned = backward;
            backward_pruned((ceil(INITCROSS_backward(W_backward<=min_consecutive))-1)...
                :(floor(FINALCROSS_backward(W_backward<=min_consecutive))-1)) = 0;
            INITCROSS_backward(W_backward<=min_consecutive) = [];
            FINALCROSS_backward(W_backward<=min_consecutive) = [];
            W_backward(W_backward<=min_consecutive) = [];
            if isequal(size(W_backward), size(INITCROSS_backward), size(FINALCROSS_backward))
                fprintf(['Backward Events #: ' num2str(size(W_backward, 1)) '\n']);
            else
                fprintf(['Forward Unequal: ' num2str(size(W_backward)) ' ' num2str(size(INITCROSS_backward)) ' ' num2str(size(FINALCROSS_backward)) '\n']);
            end
            blocks_forward = ~isempty(W_forward)*size(W_forward,1);
            blocks_backward = ~isempty(W_backward)*size(W_backward,1);
            % blocks_total = blocks_forward + blocks_backward;

            % Store signals by forward and backward movements separately
            % Forward
            if ~(blocks_forward==0)
                signal_directions_separated_forward = cell(blocks_forward,2);
                for m = 1:blocks_forward
                    signal_directions_separated_forward{m,1} = signal_anterior(ceil(INITCROSS_forward(m))-1:floor(FINALCROSS_forward(m))-1);
                    signal_directions_separated_forward{m,2} = signal_posterior(ceil(INITCROSS_forward(m))-1:floor(FINALCROSS_forward(m))-1);
                end
            else
                signal_directions_separated_forward = {};
            end
            % Backward
            if ~(blocks_backward==0)
                signal_directions_separated_backward = cell(blocks_backward,2);
                for m = 1:blocks_backward
                    signal_directions_separated_backward{m,1} = signal_anterior(ceil(INITCROSS_backward(m))-1:floor(FINALCROSS_backward(m))-1);
                    signal_directions_separated_backward{m,2} = signal_posterior(ceil(INITCROSS_backward(m))-1:floor(FINALCROSS_backward(m))-1);
                end
            else
                signal_directions_separated_backward = {};
            end
            % Cell array for storage of pooled signal data
            % First two store forward signals, second two store backward signals
            % signal_directions = cell(4,1);
            % Neurons signals during forward
            % signal_anterior_forward = signal_anterior .* forward_pruned;
            % signal_directions{1,1} = signal_anterior_forward(signal_anterior_forward>0);
            % signal_posterior_forward = signal_posterior .* forward_pruned;
            % signal_directions{2,1} = signal_posterior_forward(signal_posterior_forward>0);
            % Neurons signals during backward
            % signal_anterior_backward = signal_anterior .* backward_pruned;
            % signal_directions{3,1} = signal_anterior_backward(signal_anterior_backward>0);
            % signal_posterior_backward = signal_posterior .* backward_pruned;
            % signal_directions{4,1} = signal_posterior_backward(signal_posterior_backward>0);

            % Figure for velocity
            close all;
            s(1) = figure;
            hold on;
            plot(vel_ap_sign_smd_mean, 'k', 'linewidth', 2);
            plot(1:length(vel_ap_sign_smd_mean), zeros(length(vel_ap_sign_smd_mean),1), 'k', 'linewidth', 1);
            axis tight; set(gca, 'xticklabel', []); hold off;
            % Figure highlighted with forward locomotion
            clrmap = lbmap(10);
            s(2) = figure;
            if ~(blocks_forward==0)
                for blocks = 1:length(W_forward)
                hold on;
                st = ceil(INITCROSS_forward(blocks))-1; en = floor(FINALCROSS_forward(blocks))-1;
                fill([st en en st], 1.5*[0 0 ymax ymax], 0.8*[1 1 1], 'edgecolor', 'none');
                end
                plot(signal_anterior, 'color', clrmap(2,:), 'linewidth', 2);
                plot(signal_posterior, 'color', clrmap(end-1,:), 'linewidth', 2);
                axis tight; set(gca, 'xticklabel', []); hold off;
            else
                plot(signal_anterior, 'color', clrmap(2,:), 'linewidth', 2);
                plot(signal_posterior, 'color', clrmap(end-1,:), 'linewidth', 2);
                axis tight; set(gca, 'xticklabel', []); hold off;
            end
            % Figure highlighted with backward locomotion
            s(3) = figure;
            if ~(blocks_backward==0)
                for blocks = 1:length(W_backward)
                hold on;
                st = ceil(INITCROSS_backward(blocks))-1; en = floor(FINALCROSS_backward(blocks))-1;
                fill([st en en st], 1.5*[0 0 ymax ymax], 0.8*[1 1 1], 'edgecolor', 'none');
                end
                plot(signal_anterior, 'color', clrmap(2,:), 'linewidth', 2);
                plot(signal_posterior, 'color', clrmap(end-1,:), 'linewidth', 2);
                axis tight; set(gca, 'xticklabel', []); hold off;
            else
                plot(signal_anterior, 'color', clrmap(2,:), 'linewidth', 2);
                plot(signal_posterior, 'color', clrmap(end-1,:), 'linewidth', 2);
                axis tight; set(gca, 'xticklabel', []); hold off;
            end
            % Figure aligned with initial directional change
            s(4) = figure;
            if ~(blocks_forward==0)
                subplot(121); hold on;
                for j = 1:blocks_forward
                    a = signal_directions_separated_forward{j,1};
                    b = signal_directions_separated_forward{j,2};
                    plot(a/a(1), 'color', clrmap(2,:), 'linewidth', 2);
                    plot(b/b(1), 'color', clrmap(end-1,:), 'linewidth', 2);
                    set(gca, 'xticklabel', []);
                end
                hold off; 
            else
                set(gca, 'xticklabel', []);
            end
            if ~(blocks_backward==0)
                subplot(122); hold on
                for j = 1:blocks_backward
                    a = signal_directions_separated_backward{j,1};
                    b = signal_directions_separated_backward{j,2};
                    plot(a/a(1), 'color', clrmap(2,:), 'linewidth', 2);
                    plot(b/b(1), 'color', clrmap(end-1,:), 'linewidth', 2);
                    set(gca, 'xticklabel', []);
                end
                hold off; 
            else
                set(gca, 'xticklabel', []);
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Save data
            parts = strsplit(pathname, '\');
            data_path = fullfile(parts{1,1:end-3}, 'Alpha_Data_Raw', parts{1,end-1});
            warning('off'); mkdir(data_path); 
            data_path_name = fullfile(data_path, [filename(1:end-4) '_velocity.mat']);
            save(data_path_name, ...
                'signal_directions_separated_forward', 'signal_directions_separated_backward',...
                'vel_anterior_sign', 'vel_posterior_sign', ...
                'vel_anterior_sign_smd', 'vel_posterior_sign_smd', 'vel_ap_sign_smd_mean');
            % data_path_new = fullfile(data_path, 'Alpha_Data_Raw', 'Muscle_Interneurons_Ablated');
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