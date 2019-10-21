maxi = 15000;

curv_pos = curvdatafiltered>0.5;
curv_neg = curvdatafiltered<-0.5;

dorsal_pos = dorsal_smd .* curv_pos(1:99, :);
dorsal_neg = dorsal_smd .* curv_neg(1:99, :);
figure; colormap(jet);
subplot(2,2,1); 
image_format(dorsal_neg, 0, 1, maxi, dorsal_smd, ventral_smd); title('Dorsal Contraction', 'fontsize', 20);
subplot(2,2,3);
image_format(dorsal_pos, 0, 1, maxi, dorsal_smd, ventral_smd); title('Dorsal Extension', 'fontsize', 20);

ventral_pos = ventral_smd .* curv_pos(1:99, :);
ventral_neg = ventral_smd .* curv_neg(1:99, :);
subplot(2,2,2);
image_format(ventral_neg, 0, 1, maxi, dorsal_smd, ventral_smd); title('Ventral Extension', 'fontsize', 20);
subplot(2,2,4);
image_format(ventral_pos, 0, 1, maxi, dorsal_smd, ventral_smd); title('Ventral Contraction', 'fontsize', 20);

cbh = findobj(0, 'tag', 'Colorbar' ); 
delete(cbh);
colormap(jet);

dp_dn = [sum(sum(dorsal_pos(35:end, :)))/sum(sum(curv_pos(35:99, :))); sum(sum(dorsal_neg(35:end, :)))/sum(sum(curv_neg(35:99, :)))];
vp_vn = [sum(sum(ventral_pos(35:end, :)))/sum(sum(curv_pos(35:99, :))); sum(sum(ventral_neg(35:end, :)))/sum(sum(curv_neg(35:99, :)))];
dp_dn_vp_vn = [dp_dn ; vp_vn]