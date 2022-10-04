%%

centerline2curv;

%%

centerline2curv_pos_neg;

%%

f = dir(fullfile(pwd, '*_pos_neg.mat'));
numsamples = numel(f);
numcurvpts = 98;
headbody = 35;
numbody = numcurvpts-headbody+1;

totalposneg = zeros(numbody*numsamples, 2);
totalposnegbody = zeros(numsamples, 2);

for idx = 1:numsamples
    
    pathname = fullfile(pwd, f(idx).name);
    load(pathname, 'pos', 'neg', 'pos_total', 'neg_total');
    totalposneg(((idx-1)*numbody+1):idx*numbody,1) = pos(headbody:end);
    totalposneg(((idx-1)*numbody+1):idx*numbody,2) = neg(headbody:end);
    totalposnegbody(idx,1) = pos_total;
    totalposnegbody(idx,2) = neg_total;
    
end

totalposneg = totalposneg/100; % curvdatafiltered was multiplied by 100 from original angle data
totalposnegbody = totalposnegbody/100;

save('Curvature_pos_neg.mat', 'totalposneg', 'totalposnegbody');

%%

minic = -15; maxic = 15; 
body = 35;

imagesc(curvdatafiltered); 
caxis([minic maxic]); colorbar; 
colormap(viridis);
set(gca, 'ticklength', [0 0], 'xticklabel', [], 'yticklabel', []);
hold on;
plot([1 size(curvdatafiltered,2)], [body body], ':w', 'linewidth', 2);
pbaspect([3 1 1]);