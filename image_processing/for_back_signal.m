function [] = for_back_signal(filename, gfp, rfp)

vel = xlsread(filename);
vel = vel(:, 1);

rPos = vel>0;
rNeg = vel<0;

framenum = size(gfp, 1);
gfpPosMean = zeros(size(gfp, 2), 1);
gfpNegMean = zeros(size(gfp, 2), 1);

for i = 1:size(gfp, 2)

    ratio = gfp(:, i) ./ rfp(:, i);
    gfpPos = ratio .* rPos(1:framenum);
    gfpNeg = ratio .* rNeg(1:framenum);
    
    gfpPosMean(i) = nanmean(gfpPos(gfpPos>0));
    gfpNegMean(i) = nanmean(gfpNeg(gfpNeg>0));

end

disp(gfpPosMean);
disp(gfpNegMean);

end