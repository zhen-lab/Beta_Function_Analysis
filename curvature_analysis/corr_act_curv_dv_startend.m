% frmstart = 1;
% frmend = 60;
% headbody = 35;
% thresh = 0.1;
% 
% dorsal_ratio = dorsal_smd./dorsal_smd_r;
% ventral_ratio = ventral_smd./ventral_smd_r;
% dorsal_start = dorsal_ratio(headbody:end,frmstart);
% dorsal_end = dorsal_ratio(headbody:end,frmend);
% ventral_start = ventral_ratio(headbody:end,frmstart);
% ventral_end = ventral_ratio(headbody:end,frmend);
% curv_start = curvdatafiltered(headbody:end-1,frmstart);
% curv_end = curvdatafiltered(headbody:end-1,frmend);
% 
% corr_dorsal_start = corr(dorsal_start, curv_start);
% corr_ventral_start = corr(ventral_start, -curv_start);
% corr_dorsal_end = corr(dorsal_end, curv_end);
% corr_ventral_end = corr(ventral_end, -curv_end);
% 
% corr_total = [corr_dorsal_start corr_dorsal_end; ...
%     corr_ventral_start corr_ventral_end];

%% Determine the curvature direction

filename = uigetfile('.mat');

if isempty(filename)
    
    fprintf('No file is chosen. \n');

else
    
    load(filename);

    close all;
    curvdataBody = curvdatafiltered;

    figure;
    subplot(141); hold on;
    plot(dorsal_data{1,1}, dorsal_data{2,1}, 'r');
    plot(ventral_data{1,1}, ventral_data{2,1}, 'b');
    plot(centerline_data_spline(:,1), centerline_data_spline(:,2), 'k');
    plot(centerline_data_spline(1,1), centerline_data_spline(1,2), 'og');
    plot(centerline_data_spline(end,1), centerline_data_spline(end,2), 'oy');
    set(gca, 'ydir', 'reverse');
    subplot(142); imagesc(curvdataBody); title('Curvature');
    subplot(143); imagesc(dorsal_smd./dorsal_smd_r); title('Dorsal activity');
    subplot(144); imagesc(ventral_smd./ventral_smd_r); title('Ventral activity');

    answer = questdlg('Flip curvature data?', 'Sign of curvature', ...
        'No', 'Yes', 'No');

    if isequal(answer,'Yes')
        curvdataBody = -curvdatafiltered;
        fprintf('Curvature data is flipped. \n');
        subplot(142); imagesc(curvdataBody); title('Curvature after flipping');
    else
        fprintf('Curvature data remains unchanged. \n');
    end
    
end

%% Calculate correlations between D/V curvature and activity

thresh = 0.1;
body = 35;
istart = 1;
iend = 60;
range = [istart];
% iend = 1460;
% range = [1:780 970:iend];

curvdataBodyAdj = curvdataBody(1:end-1, range);
shapeD = curvdataBodyAdj > thresh;
shapeV = curvdataBodyAdj < -thresh;

curvdatafilteredCutD = shapeD .* curvdataBodyAdj;
curvdatafilteredCutV = shapeV .* curvdataBodyAdj;
activityCutD = shapeD .* (dorsal_smd(:,range)./dorsal_smd_r(:,range));
activityCutV = shapeV .* (ventral_smd(:,range)./ventral_smd_r(:,range));

bodycurvD = curvdatafilteredCutD(body:end,:);
bodycurvV = curvdatafilteredCutV(body:end,:);
bodyactD = activityCutD(body:end,:);
bodyactV = activityCutV(body:end,:);

subbodycurvD = bodycurvD(bodycurvD>0);
subbodycurvV = -bodycurvV(bodycurvV<0);
subbodyactD = bodyactD(bodyactD>0);
subbodyactV = bodyactV(bodyactV>0);

% Remove any inf numbers
% Only in activity vectors can inf numbers be found
DA = subbodyactD(~isinf(subbodyactD));
DC = subbodycurvD(~isinf(subbodyactD));
VA = subbodyactV(~isinf(subbodyactV));
VC = subbodycurvV(~isinf(subbodyactV));

corrD = corr(DC, DA);
corrV = corr(VC, VA);
corrDV = [corrD corrV];

fprintf('Correlations have been calculated. \n');

% Save data

data_path_name = fullfile(pwd, ['new_act-cur-corr_' filename(1:end-4) '.mat']);
save(data_path_name, 'corrDV', 'curvdataBody');
fprintf('data saved. \n');

%% Collect all correlation data

f = dir(fullfile(pwd, 'new_act-cur-corr*.mat'));
numsamples = numel(f);
dorsal = zeros(numsamples,1);
ventral = zeros(numsamples, 1);

for idx = 1 : numsamples
    name = fullfile(pwd, f(idx).name);
    load(name, 'corrDV');
    dorsal(idx) = corrDV(1);
    ventral(idx) = corrDV(2);
end

dv = [dorsal ventral];

pathname = strsplit(pwd, '\');
fname = pathname{end};
save([fname '_correlation_activity_curvture_all'], 'dv');