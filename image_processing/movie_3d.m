close all;

spf = 1/10;
frm_start = 1;
frm_end = 100;
bgcolor = [1 1 1];

% Camera travel path
flynum = 50;
xp = linspace(110,113,flynum); yp = linspace(80,103,flynum); zp = linspace(0,10,flynum);
xt = linspace(115,116,flynum); yt = linspace(109,115,flynum);
z_center = 10;

for n = frm_start:frm_end
    
    % Set 3D display using isosurface
    set(gcf, 'color', bgcolor);
    img_stack_sub = img_stack(end/2:end,:,(20*(n-1)+1):(20*(n-1)+20));
    I = img_stack_sub; thresh = 500; 
    p = patch(isosurface(I,thresh));
    cvertices = round(p.Vertices);
    cmpbase = I(sub2ind(size(I), cvertices(:,1), cvertices(:,2), cvertices(:,3)));
    cmp = cmpbase/max(cmpbase);
    isonormals(I,p)
    set(p,'FaceVertexCData',cmp,'FaceColor','interp','EdgeColor','none');
    view(3); axis tight
    xlim([80 170]); ylim([0 250]); zlim([0 20]);
    daspect([1 1 0.5])
    view(45,-45);
%     zoom(1.5)
    set(gca, 'visible', 'off', 'zdir', 'reverse');
    colormap hsv
    alpha(p, 0.6)
%     hlight = camlight('headlight');
%     p.AmbientStrength = 1;
%     camlight
%     lighting gouraud
    camproj perspective; % Perspective effect

    % Creat GIF
    % Camera travel along set path
    if n == frm_start;
        for j = 1:flynum
        campos([xp(j),yp(j),z_center]);
        camtarget([xt(j),yt(j),z_center]);
        drawnow;
        frame = getframe(1);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
            if j == 1
              imwrite(imind,cm,[filename '.gif'],'gif','Loopcount',inf,'delaytime',spf);
            else
              imwrite(imind,cm,[filename '.gif'],'gif','WriteMode','append','delaytime',spf);
            end
        end
    else
    % View from the set point after travel
%         j = flynum;
        campos([xp(j),yp(j),zp(j)]);
        camtarget([xt(j),yt(j),z_center]);
        drawnow;
        frame = getframe(1);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        imwrite(imind,cm,[filename '.gif'],'gif','WriteMode','append', 'delaytime', spf);    
%     zoom(1/1.5);
    end
    delete(p);
end