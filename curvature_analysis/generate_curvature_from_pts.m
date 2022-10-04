function [curvdata, curvdatafiltered] = generate_curvature_from_pts()

    numpts = 3; 
    imagewidth = 336;
    % fps = 20; 
    % spline_p = 0.001;
    
    % Select anterior point %%%%%%%%%%%%%%%%%%%%%%%%%%
    pt_a = uigetfile('.mat', 'Select ANTERIOR point'); 
    load(pt_a, 'neuron_position_data');
    if pt_a~=0
        if median(neuron_position_data(:,1))>imagewidth/2
            neuron_position_data(:,1) = neuron_position_data(:,1) - imagewidth/2;
        end
        numframes = size(neuron_position_data, 1);
        np = zeros(numframes, numpts*2);
        np(:,1:2) = neuron_position_data; % Pad positions of the first point
        curvdata = zeros(size(np, 2)/2-2, numframes);
        % lineAllFrame= zeros(numcurvpts+2, 2*numframes);
        fprintf(['Anterior point file name: ' pt_a '. \n']);
    else
        fprintf('user cancelled selection.\n');
    end
    pref = regexp(pt_a, '[A-Za-z]+\d*(?=[._])', 'match');
    % Select middel and posterior point %%%%%%%%%%%%%%%%%%%%%%%%%%
    if isempty(pref)
        fprintf('failed to find files with prefix word#.\n');
    else
        pt_m = uigetfile('.mat', 'Select MIDDLE point', [pref{1} '*.mat']); 
        load(pt_m, 'neuron_position_data', 'signal', 'signal_mirror');
        if median(neuron_position_data(:,1))>imagewidth/2
            neuron_position_data(:,1) = neuron_position_data(:,1) - imagewidth/2;
        end
        np(:,3:4) = neuron_position_data; r = signal{1,1}./signal_mirror{1,1};
        fprintf(['Middle point file name: ' pt_m '. \n']);

        pt_p = uigetfile('.mat', 'Select POSTERIOR point', [pref{1} '*.mat']); 
        load(pt_p, 'neuron_position_data');
        if median(neuron_position_data(:,1))>imagewidth/2
            neuron_position_data(:,1) = neuron_position_data(:,1) - imagewidth/2;
        end
        np(:,5:6) = neuron_position_data;
        fprintf(['Posterior point file name: ' pt_p '. \n']);
    end

    % Draw activity curvature plot %%%%%%%%%%%%%%%%%%%%%%%%%%
%     if ~isequal(pt_a(1:end-7),pt_m(1:end-7),pt_p(1:end-7))
%         warning('Files are not from the same recording. \n');
%     else
    fprintf('Drawing started.\n');
    fname = pt_m(1:end-4);
    for i = 1:numframes

        line = np(i, :);
        line = reshape(line, 2, [])';

    %     [splineLine, cv2, sp_curv_num] = spline_line(line, spline_p, numcurvpts);
    %     cv2i = interp1(sp_curv_num+.00001*(0:length(sp_curv_num)-1), cv2, (0:(sp_curv_num(end)-1)/(numcurvpts+1):(sp_curv_num(end)-1)));    
    %     lineAllFrame(:, 2*i-1:2*i) = cv2i;

        df2 = diff(line, 1, 1);
        atdf2 = unwrap(atan2(-df2(:, 2), df2(:, 1))); % atan2(Y, X)
    %     curv = diff(atdf2, 1);
    %     curvdata(:, i) = atdf2';
        curv = unwrap(diff(atdf2, 1)); 
        curvdata(:, i) = curv;

    %     hold off;
    %     figure(1); imagesc(imagelist{i, 1}); axis equal; axis off; hold on;
    %     plot(line(:, 1), line(:, 2), 'w', 'LineWidth', 3); hold on;
    %     plot(splineLine(:, 2), splineLine(:, 1), 'g', 'LineWidth', 1.5); hold on;
    %     plot(cv2i(:, 2), cv2i(:, 1), 'r', 'LineWidth', 0.8); hold on;
    %     text(10, 10, ['Frame #' num2str(i)], 'Color', 'w');

    end

    %%%%%Filter curvature by averaging%%%%%%%

    timefilter = 2;
    bodyfilter = 10;

    h = fspecial('average', [timefilter bodyfilter]);
    curvdatafiltered = imfilter(curvdata*100,  h , 'replicate');
    curvdatafiltered = smooth(curvdatafiltered);
    curvdatafiltered = reshape(curvdatafiltered, [], numframes);
    % maxcurv = max(max(curvdatafiltered));
    % mincurv = min(min(curvdatafiltered));

    % figure(2); imagesc(curvdatafiltered); colorbar; caxis([-100 100]);
    % 
    % title('Cuvature Diagram');
    % 
    % set(gca, 'YDir', 'normal'); % To display the head-tail direction from bottom to top
    % set(gca,'YTICK',[1 20 40 60 80 100]);
    % set(gca,'YTICKLABEL',[0 0.2 0.4 0.6 0.8 1]);
    % 
    % set(gca,'XTICK',1:10*fps:numframes);
    % x_tick=get(gca,'XTICK');
    % set(gca,'XTICKLABEL',(x_tick-1)/fps);
    % 
    % ylabel('Fractional distance along the centerline/head=0,tail=1');
    % xlabel('Time/s');
    % plotyy(1:length(r), r, 1:length(r), curvdatafiltered);
    figure; 

    subplot(211); title('Positive curvdata');
    yyaxis left; plot(r, 'color', [1 0.5 0]); ylabel('activity');
    yyaxis right; plot(curvdatafiltered/100, 'color', [0 0.5 1]); ylabel('curvature');
    xl = xlim;
    hold on; plot([xl(1) xl(end)], [0 0], 'k'); axis tight;

    subplot(212); title('Negative curvdata');
    yyaxis left; plot(r, 'color', [1 0.5 0]); ylabel('activity');
    yyaxis right; plot(-curvdatafiltered/100, 'color', [0 0.5 1]); ylabel('curvature');
    xl = xlim;
    hold on; plot([xl(1) xl(end)], [0 0], 'k'); axis tight;

    posneg = questdlg('Postive or negative curvdata?',...
        '+ or -', 'Positive', 'Negative', 'Positive');
    switch posneg
        case 'Positive'
            pn = 1;
        case 'Negative'
            pn = -1;
    end

    curvdata = pn * curvdata';
    curvdatafiltered = pn * curvdatafiltered';

    save([fname '_curvature_curated.mat'], 'curvdata', 'curvdatafiltered');
    fprintf('curvdata curated and saved. \n');

%     end

end