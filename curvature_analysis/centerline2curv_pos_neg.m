filename = uigetfile('.mat', ...
    'Select excel files for curvature', '*_curvature_curated.mat');
    
load(filename, 'curvdatafiltered');
c = curvdatafiltered;
headbody = 35;
c_body = c(headbody:end,:);

pos = nansum((c>0).*c, 2);
pos_total = nansum((c_body>0).*c_body, 'all');
neg = abs(nansum((c<0).*c, 2));
neg_total = abs(nansum((c_body<0).*c_body, 'all'));

save([filename(1:end-5) '_pos_neg.mat'], 'pos', 'pos_total', 'neg', 'neg_total');