function [curvdata, curvdatafiltered] = curvature_pts(np)

% spline_p = 0.001;
fps = 20;
numframes = size(np, 1);
curvdata = zeros(size(np, 2)/2-2, numframes);
% lineAllFrame= zeros(numcurvpts+2, 2*numframes);

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

figure(2); imagesc(curvdatafiltered); colorbar; caxis([-100 100]);

title('Cuvature Diagram');

set(gca, 'YDir', 'normal'); % To display the head-tail direction from bottom to top
set(gca,'YTICK',[1 20 40 60 80 100]);
set(gca,'YTICKLABEL',[0 0.2 0.4 0.6 0.8 1]);

set(gca,'XTICK',1:10*fps:numframes);
x_tick=get(gca,'XTICK');
set(gca,'XTICKLABEL',(x_tick-1)/fps);

ylabel('Fractional distance along the centerline/head=0,tail=1');
xlabel('Time/s');


curvdata = curvdata';
curvdatafiltered = curvdatafiltered';

end