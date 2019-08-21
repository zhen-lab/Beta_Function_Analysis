function [dorsal_mean, ventral_mean] = collect_DV_data_mean(dorsal, ventral)

dorsal(dorsal==inf) = 0;
dorsal(isnan(dorsal)) = 0;
dorsal_sum = sum(dorsal,2);
dorsal_usednum = sum(dorsal>0,2);
dorsal_mean = dorsal_sum./dorsal_usednum;

ventral(ventral==inf) = 0;
ventral(isnan(ventral)) = 0;
ventral_sum = sum(ventral,2);
ventral_usednum = sum(ventral>0,2);
ventral_mean = ventral_sum./ventral_usednum;

% save('Ablated_Dorsal_Ventral_All_mean.mat', 'dorsal_mean', 'ventral_mean');

end