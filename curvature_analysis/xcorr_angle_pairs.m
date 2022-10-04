function [curvdata, curvdatafiltered, lagsall, rall] = xcorr_angle_pairs(np, gfp, rfp, sig)

[curvdata, curvdatafiltered] = curvature_pts(np);
curvdata = sig * curvdata;
rall = zeros(size(curvdata, 2), 2*size(curvdata,1)-1);
lagsall = zeros(size(curvdata, 2), 1);
signal = smoothdata(gfp(:, 3)./rfp(:, 3), 'rloess');

for i = 1:size(curvdata, 2)

    [r,lags] = xcorr(normalize_signal(signal),...
        normalize_signal(curvdata(:, i))); 
    rall(i, :) = r';
    lagsall(i) = lags(r == max(r));

end

rall = rall';
disp(lagsall);  

% Check lags
figure; 

    subplot(2,1,1); hold on;
    plot(normalize_signal(signal), 'linewidth', 4);
    plot(normalize_signal(curvdata));
    legend('Ratio', 'N-1', 'N', 'N+1');
    set(gca, 'xlim', [1+min(lagsall) length(gfp)+max(lagsall)]);
    
    subplot(2,1,2); hold on;
    plot(normalize_signal(signal), 'linewidth', 4);
    for i = 1:size(curvdata, 2)
       plot(1+lagsall(i):length(gfp)+lagsall(i), normalize_signal(curvdata(:, i)));
    end
    set(gca, 'xlim', [1+min(lagsall) length(gfp)+max(lagsall)]);
    hold off;

end