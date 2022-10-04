h(1) = figure(1);
imagesc(curvdatafiltered);
colormap(viridis); 
h(2) = figure(2);
imagesc(dorsal_smd./dorsal_smd_r);
colormap(magma); 
h(3) = figure(3);
imagesc(ventral_smd./ventral_smd_r);
colormap(magma); 
