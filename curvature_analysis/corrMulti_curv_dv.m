%% Determine the curvature direction

filename = uigetfile('.mat');
load(filename);

close all;
curvdataBody = curvdatafiltered;

figure;
subplot(131); imagesc(curvdataBody); title('Curvature');
subplot(132); imagesc(dorsal_smd); title('Dorsal activity');
subplot(133); imagesc(ventral_smd); title('Ventral activity');

answer = questdlg('Flip curvature data?', 'Sign of curvature', ...
    'No', 'Yes', 'No');

if isequal(answer,'Yes')
    curvdataBody = -curvdatafiltered;
    fprintf('Curvature data is flipped. \n');
    subplot(131); imagesc(curvdataBody); title('Curvature after flipping');
else
    fprintf('Curvature data remains unchanged. \n');
end

%% Calculate correlations between D/V curvature and activity

close all;

thresh = 1.5;
threshIn = 50;
body = 35;

% curvdataBody = curvdatafiltered;
dorsalBody = dorsal_smd;
ventralBody = ventral_smd;

curvdatafilteredCutD = (curvdataBody>0) .* curvdataBody;
curvdatafilteredCutV = (curvdataBody<0) .* curvdataBody;

cD = contour(curvdatafilteredCutD, [thresh thresh]);
cV = contour(curvdatafilteredCutV, -[thresh thresh]);

sD = getcontourlines(cD);
sV = getcontourlines(cV);

curvdataBody = curvdataBody(1:end-1, :);
col_num = size(curvdataBody, 1);
row_num = size(curvdataBody, 2);
xq = 1:row_num; xq_rep = repmat(xq, col_num, 1); xq_rep_shape = reshape(xq_rep, [], 1);
yq = 1:col_num; yq_rep = repmat(yq, row_num, 1); yq_rep_shape = reshape(yq_rep', [], 1);

j = 1;

for i = 1:length(sD)
    
    line = [sD(i).x; sD(i).y];
    in = inpolygon(xq_rep_shape, yq_rep_shape, line(1, :), line(2, :));

    logi_all = reshape(in, col_num, row_num);
    logi = logi_all(body:end, :);

    if sum(sum(logi)) >= threshIn

        dorsalIn = logi.*dorsalBody(body:end, :);
        curvIn = logi.*curvdataBody(body:end, :);
        figure(j);
        subplot(2,1,1); imagesc(dorsalIn); title('Dorsal activity');
        subplot(2,1,2); imagesc(curvIn); title('Dorsal curvature');

        dorsalInPos = dorsalIn(logi==1);
%         dorsalInFor = dorsalIn(:, fwd); dorsalInForPos = dorsalInFor(logi(:, fwd)==1);
%         dorsalInBak = dorsalIn(:, bwd); dorsalInBakPos = dorsalInBak(logi(:, bwd)==1);

        curvInPos = curvIn(logi==1);
%         curvInFor = curvIn(:, fwd); curvInForPos = curvInFor(logi(:, fwd)==1);
%         curvInBak = curvIn(:, bwd); curvInBakPos = curvInBak(logi(:, bwd)==1);

        ccAllD(j, :) = corr(dorsalInPos, curvInPos);
%         fd = corrcoef(dorsalInForPos, curvInForPos); if numel(fd) == 1, ccForD(j, :) = 1; else ccForD(j, :) = fd(2); end;
%         bd = corrcoef(dorsalInBakPos, curvInBakPos); if numel(bd) == 1, ccBakD(j, :) = 1; else ccBakD(j, :) = bd(2); end;

        j = j+1;

    end
   
end

beep; pause; 
close all;

j = 1;

for i = 1:length(sV)
    
    line = [sV(i).x; sV(i).y];
    in = inpolygon(xq_rep_shape, yq_rep_shape, line(1, :), line(2, :));

    logi_all = reshape(in, col_num, row_num);
    logi = logi_all(body:end, :);

    if sum(sum(logi)) >= threshIn

        ventralIn = logi.*ventralBody(body:end, :);
        curvIn = -logi.*curvdataBody(body:end, :);
        figure(j);
        subplot(2,1,1); imagesc(ventralIn); title('Ventral activity');
        subplot(2,1,2); imagesc(curvIn); title('Ventral curvature');

        ventralInPos = ventralIn(logi==1);
%         ventralInFor = ventralIn(:, fwd); ventralInForPos = ventralInFor(logi(:, fwd)==1);
%         ventralInBak = ventralIn(:, bwd); ventralInBakPos = ventralInBak(logi(:, bwd)==1);

        curvInPos = curvIn(logi==1);
%         curvInFor = curvIn(:, fwd); curvInForPos = curvInFor(logi(:, fwd)==1);
%         curvInBak = curvIn(:, bwd); curvInBakPos = curvInBak(logi(:, bwd)==1);

        ccAllV(j, :) = corr(ventralInPos, curvInPos);
%         fv = corrcoef(ventralInForPos, curvInForPos); if numel(fv) == 1, ccForV(j, :) = 1; else ccForV(j, :) = fv(2); end;
%         bv = corrcoef(ventralInBakPos, curvInBakPos); if numel(bv) == 1, ccBakV(j, :) = 1; else ccBakV(j, :) = bv(2); end;

        j = j+1;

    end

end

beep;

fprintf('Correlations have been calculated. \n');

% Save data

data_path_name = fullfile(pwd, ['act-cur-corr_' filename(1:end-4) '.mat']);
save(data_path_name, 'ccAllD', 'ccAllV', 'curvdataBody');
fprintf('data saved. \n');