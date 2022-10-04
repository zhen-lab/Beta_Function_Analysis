function [curvdatafiltered, nan_frame] = midline2curv(filename)

filename_string = char(filename);
alldata = xlsread(filename_string);

vel = alldata(:, 1); % Velocity information stored in col 1
nan_frame = isnan(vel); % Get the NaN based on velocity
pos = alldata(:, 39:104); % Position information stored from col 39 to 104

spline_p = 0.002;
numcurvpts = size(pos,2)/2;
framenum = size(pos, 1);
% fps = 20;
curv_data = zeros(numcurvpts, framenum);

for i = 1:framenum
        
        pos_ind = reshape(pos(i,:), 2, [])';
        [~, cv2, sp_curv_num] = spline_line(pos_ind, spline_p, numcurvpts);

        cv2i = interp1(sp_curv_num+.00001*(0:length(sp_curv_num)-1), cv2, (0:(sp_curv_num(end)-1)/(numcurvpts+1):(sp_curv_num(end)-1)));
        df2 = diff(cv2i, 1, 1);
        atdf2 = unwrap(atan2(-df2(:,2), df2(:,1)));

        % Collect curvature information for midline
        curv = unwrap(diff(atdf2, 1));
        curv_data(:, i) = curv;
    
end

curv_data = reshape(smooth(curv_data), [], framenum);

timefilter = 2;
bodyfilter = 10;
h = fspecial('average', [timefilter bodyfilter]);

curvdatafiltered = imfilter(curv_data*10,  h , 'replicate');
curvdatafiltered = smooth(curvdatafiltered);
curvdatafiltered = reshape(curvdatafiltered, [], framenum);

figure('units','normalized','outerposition',[0 0 1 1]);
curvdatafiltered(:, nan_frame) = NaN;
[nr, nc] = size(curvdatafiltered);
pcolor([curvdatafiltered nan(nr,1); nan(1,nc+1)]);
shading flat;

colormap jet;
c = colorbar;
caxis([-3 3]); set(c, 'ytick', []);
% ax = gca;
% axpos = get(ax, 'Position');
% cpos = get(c, 'Position');
% cpos(3) = 0.1;
% set(c, 'Position', cpos);
% set(ax, 'Position', axpos);

% title('Cuvature Diagram');
% 
% set(gca, 'YDir', 'normal'); % To display the head-tail direction from bottom to top
% set(gca,'YTICK',[1 11 22 33]);
% set(gca,'YTICKLABEL',[0 0.33 0.67 1]);
% 
% set(gca,'XTICK',1:10*fps:framenum);
% x_tick=get(gca,'XTICK');
% set(gca,'XTICKLABEL',(x_tick-1)/fps);
% 
% ylabel('Head = 0,Tail = 1');
% xlabel('Time/s');

% set(gca, 'visible', 'off');
set(gca, 'xticklabel', [], 'yticklabel', [], 'layer', 'top', 'ytick', [], 'ticklength', [0 0]);
daspect([10,1, 1]);

end