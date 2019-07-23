function [dorsal_smd, ventral_smd, dorsal_smd_r, ventral_smd_r] = activity_all(imagelist_g, imagelist_r, range, dorsal_data, ventral_data, centerline_data_spline, curvdatafiltered)

[dorsal_smd, ~] = activity_wave(imagelist_g(range(1):range(2)), dorsal_data, centerline_data_spline, curvdatafiltered, 0.35, 1, 0, 0);
[ventral_smd, ~] = activity_wave(imagelist_g(range(1):range(2)), ventral_data, centerline_data_spline, curvdatafiltered, 0.35, 1, 0, 0);
[dorsal_smd_r, ~] = activity_wave(imagelist_r(range(1):range(2)), dorsal_data, centerline_data_spline, curvdatafiltered, 0.35, 1, 0, 0);
[ventral_smd_r, ~] = activity_wave(imagelist_r(range(1):range(2)), ventral_data, centerline_data_spline, curvdatafiltered, 0.35, 1, 0, 0);

close all; 
figure; headneck = 36;
subplot(1,2,1); title('Absolute'); hold on;
plot(mean(dorsal_smd(headneck:end, :)./dorsal_smd_r(headneck:end, :)), 'r');
plot(mean(ventral_smd(headneck:end, :)./ventral_smd_r(headneck:end, :)), 'b');
hold off;
subplot(1,2,2); title('Relative'); hold on;
plot(mean(dorsal_smd(headneck:end, :)./dorsal_smd_r(headneck:end, :)/mean(dorsal_smd(headneck:end, 1)./dorsal_smd_r(headneck:end,1))), 'r');
plot(mean(ventral_smd(headneck:end, :)./ventral_smd_r(headneck:end, :)/mean(ventral_smd(headneck:end, 1)./ventral_smd_r(headneck:end,1))), 'b');
hold off;

end