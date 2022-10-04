function [curvdata, curvdatafiltered] = centerline2curv()

    filename = uigetfile('*.mat', 'Select a curated file for curvature', ...
        '*_velocity_curated.mat');
    load(filename, 'centerline_corrected');
    curvdata = zeros(size(centerline_corrected,1)-2, size(centerline_corrected,2)/2); 
    for i = 1:size(centerline_corrected,2)/2
        df2 = diff(centerline_corrected(:,2*i-1:2*i), 1, 1);
        atdf2 = unwrap(atan2(-df2(:,2), df2(:,1)));
        % Collect curvature information for midline
        curv = unwrap(diff(atdf2, 1));
        curvdata(:, i) = curv;
    end

    timefilter = 2; bodyfilter = 10;
    h = fspecial('average', [timefilter bodyfilter]);
    curvdatafiltered = imfilter(curvdata*100,  h , 'replicate');
    curvdatafiltered = smoothdata(curvdatafiltered, 1);
    imagesc(curvdatafiltered);
    
    save([filename(1:end-4) '_curvature_curated.mat'], ...
    'curvdata', 'curvdatafiltered');

end