%%
load('Image_processing_temp14.tif-temp15.tif.mat')
%%
filename = 'temp14.mat';
load(filename);
j = 1;
centerline_data_spline = mst_centerline_data_spline{j,1};
curvdatafiltered = mst_curvdatafiltered{j,1};
dorsal_data = mst_dorsal_data{j,1}; 
ventral_data = mst_ventral_data{j,1}; 
save([filename(1:end-4) '_updated.mat'], ...
    'dorsal_data', 'ventral_data', 'centerline_data_spline', 'curvdatafiltered', 'dorsal_smd', 'dorsal_smd_r', 'ventral_smd', 'ventral_smd_r');