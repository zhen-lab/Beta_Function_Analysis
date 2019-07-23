function [] = create_gif(filename, img_stack, numframe, lighton)

close all;
bgcolor = [236 240 241]/255;
figure(1)
set(gcf, 'color', bgcolor);
filename = [filename(1:end-4) '.gif'];
spf = 68/1800;

for n = 1:numframe
    imagesc(img_stack(:,:,n)); colormap gray; caxis([0 255]);
    set(gca, 'visible', 'off');
    text(230, 230, [num2str(n*spf, '%0.1f') 'sec'], 'fontsize', 18, 'color', 'k');
    
    if n>=lighton(1) && n<=lighton(2)
%         text(230, 230, 'Light On', 'fontsize', 18, 'color', 'g');
        colormap summer;
    end;
    
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