%% Pre-prcoess data

% Load tiff files
setup_proof_reading;
fprintf('tiff files loading finished. \n');

%% Register RFP and GFP channels
figure; 
subplot(1,2,1); imshowpair(imagelist_g{1,1}, imagelist_r{1,1});
image_registration_tform;
imagelist_g = movingRegistered;
subplot(1,2,2); imshowpair(imagelist_g{1,1}, imagelist_r{1,1});
%%
close all;

% Prepare filters
bump_finder = bump_detector([1,1], [3,3], [51 51]);
normalization = 0.7 * sum_all(bump_finder(bump_finder>0)*255);
dc = @(x) deconvlucy(x, fspecial('gaussian', 5, 2), 1);
flt = @(x, h) imfilter(double(dc(x)), h, 'circular');
flt_and_normalize = @(x) uint8(flt(x, bump_finder)/normalization*256-5);

% Misc variables
fps = 20; dt = 0.8/fps; imagerange = 1:700;
imagerangex = size(imagelist_g{1,1},2); imagerangey = 700; imagerangez = 300;
edgelength = 2*length(dorsal_data{1,1});
textpos = [10 400 100];
colormap hot;
cmp = colormap;
cmp(1,:) = [0 0 0];
al = -54; ez = 64; coefsteps = 25;

% Loop for GIF
for frm = 1:size(imagelist,1)

    flt_g = flt_and_normalize(imagelist_g{frm,1}(imagerange,:));
    flt_r = flt_and_normalize(imagelist_r{frm,1}(imagerange,:));
    imagelogic = flt_r > mean2(flt_r);
    dorsalx = dorsal_data{2*frm-1,1}; dorsaly = dorsal_data{2*frm,1};
    ventralx = flipud(ventral_data{2*frm-1,1}); ventraly = flipud(ventral_data{2*frm,1});

    if frm == 1 % Draw initial rotation and contour separation

        for coef = linspace(1,100,coefsteps)

            hold off;
            imagesc(double(flt_g).*double(imagelogic)./double(flt_r)); caxis([0 2.5]); 
            colorbar; colormap(cmp);

            hold on;
            plot3([dorsalx;ventralx], [dorsaly;ventraly], coef*ones(edgelength), 'color', 'w', 'linewidth', 1);

            view((coef-1)/99*al, 90-(coef-1)/99*(90-ez)); axis equal; box off;
            set(gca, 'visible', 'off'); xlim([1 imagerangex]); ylim([1 imagerangey]); zlim([0 imagerangez]);
            set(gcf,'units','normalized','outerposition',[0 0 1 1], 'color', 'k');
            zoom(1.8);

            drawnow;
            frame = getframe(1);
            im = frame2im(frame);
            [imind, cm] = rgb2ind(im, 256);

            if coef == 1 % Skip the first time point for colorbar location discrenpancy
              imwrite(imind, cm, [filename '.gif'], 'gif', 'LoopCount', 1, 'DelayTime', dt);
            else
              imwrite(imind, cm, [filename '.gif'], 'gif', 'WriteMode', 'append', 'DelayTime', dt);
            end

            if coef == 100 % Add a brief waiting period between two stages
                text(textpos(1), textpos(2), textpos(3), [num2str((frm-1)/fps, '%.2f') ' s'], 'color', 'w');
                scatterbar3(dorsalx(1:end-1), dorsaly(1:end-1), coef*dorsal_smd(:,frm)./dorsal_smd_r(:,frm), 0.5, 0);
                scatterbar3(ventralx(1:end-1), ventraly(1:end-1), coef*ventral_smd(:,frm)./ventral_smd_r(:,frm), 0.5, 1);
                drawnow;
                frame = getframe(1);
                im = frame2im(frame);
                [imind, cm] = rgb2ind(im, 256);
                imwrite(imind, cm, [filename '.gif'], 'gif', 'WriteMode', 'append', 'DelayTime', 3*dt);
            end

        end

    else % Draw actual movement

        hold off;
        imagesc(double(flt_g).*double(imagelogic)./double(flt_r)); caxis([0 2.5]); 
        text(textpos(1), textpos(2), textpos(3), [num2str((frm-1)/fps, '%.2f') ' s'], 'color', 'w');
        colorbar; colormap(cmp);

        hold on;
        coef = 100;
        plot3([dorsalx;ventralx], [dorsaly;ventraly], coef*ones(edgelength), 'color', 'w', 'linewidth', 1);
        scatterbar3(dorsalx(1:end-1), dorsaly(1:end-1), coef*dorsal_smd(:,frm)./dorsal_smd_r(:,frm), 0.5, 0);
        scatterbar3(ventralx(1:end-1), ventraly(1:end-1), coef*ventral_smd(:,frm)./ventral_smd_r(:,frm), 0.5, 1);

        view(al, ez); axis image; box off;
        set(gca, 'visible', 'off'); 
        xlim([1 imagerangex]); ylim([1 imagerangey]); zlim([0 imagerangez]);
        set(gcf,'units','normalized','outerposition',[0 0 1 1], 'color', 'k');
        zoom(1.8);

        drawnow;
        frame = getframe(1);
        im = frame2im(frame);
        [imind, cm] = rgb2ind(im, 256);
        imwrite(imind, cm, [filename '.gif'], 'gif', 'WriteMode', 'append', 'DelayTime', dt);

    end

end