filename = 'AVA_Activation_Annotated.gif';
h = figure; fps = 5; fast = 4;

for i = 1:100
    
    I = imread(['img-' num2str(i-1, '%03.0f') '.png']);
    image(I); axis image;
    set(gca, 'visible', 'off'); set(gcf, 'color', 'w');
    hold on; 
    text(15, 15, [num2str(i/fps) 'sec'], 'color', 'w'); drawnow; 
    hold off;
    
    frame = getframe(h); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 
    % Write to the GIF File 
    if i == 1 
      imwrite(imind,cm,filename,'gif', 'Loopcount',inf, 'delaytime', 1/fps/fast); 
    else 
      imwrite(imind,cm,filename,'gif','WriteMode','append', 'delaytime', 1/fps/fast); 
    end 
end