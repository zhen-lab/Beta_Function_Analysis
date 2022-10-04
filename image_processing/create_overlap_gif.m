function [] = create_overlap_gif(filename, imagelist, lighton)

figure(1)
set(gcf, 'color', [1 1 1]);
filename = [filename(1:end-4) '_overlap.gif'];
spf = 68/1800;

for i = 1:lighton(2)-lighton(1)
    set(gca, 'visible', 'off');
    
    imshowpair(imagelist{mod(i, lighton(1))+1,1}, imagelist{i+lighton(1)-1, 1});    
       
    drawnow;
      frame = getframe(1);
      im = frame2im(frame);
      [imind,cm] = rgb2ind(im,256);
      if i == 1
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf, 'delaytime', spf);
      else
          imwrite(imind,cm,filename,'gif','WriteMode','append', 'delaytime', spf);
      end
      
end

end