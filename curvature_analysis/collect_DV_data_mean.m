function [dorsal_mean, ventral_mean] = collect_DV_data_mean(dorsal, ventral)

dorsal_sum = sum(dorsal,2);
dorsal_usednum = sum(dorsal>0,2);
dorsal_mean = dorsal_sum./dorsal_usednum;

ventral_sum = sum(ventral,2);
ventral_usednum = sum(ventral>0,2);
ventral_mean = ventral_sum./ventral_usednum;

end