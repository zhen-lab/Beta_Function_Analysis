function [F_denoised_d, F_denoised_v, ratio_d_v] = ratio_dorsal_ventral(imagelist, dorsal_data, ventral_data, centerline_data)

centerline_data(centerline_data<1) = 1; % The possibility of the smoothed line out of the box, regardless of rows/cols

frmnum = size(imagelist, 1);
F_denoised_d = zeros(frmnum, 1);
F_denoised_v = zeros(frmnum, 1);

for i = 1:frmnum
    
    roi_d = zeros(size(imagelist{1, 1}, 1), size(imagelist{1, 1}, 2));
    roi_v = zeros(size(imagelist{1, 1}, 1), size(imagelist{1, 1}, 2));
    
    line_frm = centerline_data(:, (2*i-1):(2*i));
    line_x = line_frm(:, 1); line_y = line_frm(:, 2);
    line_x(line_x>=size(imagelist{1, 1}, 1)) = size(imagelist{1, 1}, 1); % Restrain the points above the number of rows
    line_y(line_y>=size(imagelist{1, 1}, 2)) = size(imagelist{1, 1}, 2); % Restrain the points above the number of cols
    
    for l = 1:length(line_frm)        
        roi_d(ceil(line_x(l)), ceil(line_y(l)))=1;
        roi_v(ceil(line_x(l)), ceil(line_y(l)))=1;
    end
    
    %%%%%%%% Dorsal side %%%%%%%%
    
    dorsal_frm_y = dorsal_data{2*i-1, 1};
    dorsal_frm_x = dorsal_data{2*i, 1};
    
    for d = 1:length(dorsal_frm_x)
        roi_d(round(dorsal_frm_x(d)), round(dorsal_frm_y(d)))=1;
    end
    
    %%%%%%%% Ventral side %%%%%%%%

    ventral_frm_y = ventral_data{2*i-1, 1};
    ventral_frm_x = ventral_data{2*i, 1};
    
    for v = 1:length(ventral_frm_x)
        roi_v(round(ventral_frm_x(v)), round(ventral_frm_y(v)))=1;        
    end

    %%%%%%%% Close edges for dorsal and ventral side %%%%%%%%
    
    se = strel('disk', 6);
    closeBW_edge_d = imclose(roi_d, se);
    image_d = times(uint16(imagelist{i, 1}), uint16(closeBW_edge_d));
    closeBW_edge_v = imclose(roi_v, se);
    image_v = times(uint16(imagelist{i, 1}), uint16(closeBW_edge_v));
    
    %%%%%%%% Plot dorsal and ventral side %%%%%%%%
    
    hold off;
    figure(5); 
    subplot(1,2,1);
    imagesc(image_d); axis equal; axis off;
    subplot(1,2,2);
    imagesc(image_v); axis equal; axis off;
    
    image_d_pos = image_d(image_d > 0);
    [F_denoised_d(i), ~, ~] = pixel_intensity(image_d_pos, 0.25, 0.1);
    image_v_pos = image_v(image_v > 0);
    [F_denoised_v(i), ~, ~] = pixel_intensity(image_v_pos, 0.25, 0.1);
    
    ratio_d_v= F_denoised_d ./ F_denoised_v;

end

%%%%%%%% Plot everything %%%%%%%%

figure(6);

subplot(3,1,1); plot(F_denoised_d, 'Color', 'r', 'LineWidth', 2); xlabel('Frame number'); ylabel('\color{red}Dorsal \color{black}signal');
subplot(3,1,2); plot(F_denoised_v, 'Color', 'g', 'LineWidth', 2); xlabel('Frame number'); ylabel('\color{red}Ventral \color{black}signal');
subplot(3,1,3); semilogy(ratio_d_v, 'Color', 'b', 'LineWidth', 2); hold on; line(1:length(ratio_d_v), 1); xlabel('Frame number'); ylabel('Dorsal/Ventral');

beep;

end
