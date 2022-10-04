filename = uigetfile('.xlsx', 'Select excel files for curvature');
    
c = midline2curv(filename);

pos = nansum((c>0).*c, 2);
neg = abs(nansum((c<0).*c, 2));

save([filename(1:end-5) '_curv_pos_neg.mat'], 'pos', 'neg');