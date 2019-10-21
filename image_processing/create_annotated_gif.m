function [] = create_annotated_gif(filename, imagelist, numframes, ornt, clr)

set(gcf, 'color', [1 1 1]);
filename = [filename(1:end-4) '.gif'];
pos = [size(imagelist{1,1},2)/2-25 10 10 10; size(imagelist{1,1},2)/2-5 10 10 10; size(imagelist{1,1},2)/2+15 10 10 10];
% clr = [0 0.2 1; 0.5 0.5 0.5; 0.9 0.2 0];
sig = [1 0 -1];
accum = zeros(3, 1);
ord = 1:3;
str = {'D', 'N', 'V'};
spf = 68/1800;

vec = 1:24*16;
% vec = reshape(vec, 21, 24)';
vec_left = vec(:, 1:16*16); 
% vec_left = reshape(vec_left, [], 1);
vec_right = vec(:, 16*16+1:end); 
% vec_right = reshape(vec_right, [], 1);

for n = 1:numframes
    hold off;
    figure(1);
%%% Draw annotated recording %%%
    subplot(24, 16, vec_left);
    hold on;
    background = imopen(imagelist{n, 1},strel('disk',15));
    I = imagelist{n, 1}-background;
    imagesc(I);    
%     imagesc(imagelist{n, 1}); 
    colormap gray; axis equal;
    set(gca, 'visible', 'off');
    % Draw three circles
    for i = 1:size(pos, 1)
        rectangle('position', pos(i, :), 'curvature', [1 1], 'facecolor', 'none', 'edgecolor', [1 1 1]);
        text(pos(i, 1), pos(i, 2)+20, str(i), 'color', clr(i, :), 'fontweight', 'bold');
        accum(i) = length(find(ornt(1:n) == sig(i))); % Calculate currently accumulated number of frames for each orientation
    end
    % Color the circle according to orientation
    num = ord(sig == ornt(n));
    rectangle('Position', pos(num, :), 'Curvature', [1 1], 'facecolor', clr(num, :),'edgecolor', [1 1 1]);
    text(size(imagelist{1,1}, 2)-50, size(imagelist{1,1}, 1)-10, [num2str(n*spf, '%.2f') 's'], 'color', 'w', 'fontweight', 'bold');
    set(gca, 'ydir', 'normal');
    hold off;
    
%%% Draw distribution of orientation %%%
    subplot(24, 16, vec_right);
    b = bar([accum(1) 0 0; 0 accum(2) 0; 0 0 accum(3)]/sum(accum));
    for i = 1:3
        b(i).FaceColor = clr(i, :);
        b(i).EdgeColor = [1 1 1];
    end
    axis normal;
    set(gca, 'visible', 'off', 'xlim', [-1 5], 'ylim', [0 1], 'layer', 'top', 'ycolor', [1 1 1]);
    
%%% Create GIF file %%%
    drawnow;
      frame = getframe(1);
      im = frame2im(frame);
      [imind,cm] = rgb2ind(im,256);
      if n == 1;
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf, 'delaytime', spf);
      else
          imwrite(imind,cm,filename,'gif','WriteMode','append', 'delaytime', spf);
      end
    
end

end