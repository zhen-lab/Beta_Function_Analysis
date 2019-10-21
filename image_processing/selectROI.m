function [ cropy1,cropy2,cropx1,cropx2] = selectROI( img )

 [n,m]=size(img);

 text(10,10, 'select ROI: upper left then lower right', 'Color', 'white'); %define ROI, which typically contains single worm
 [cropx1 cropy1] = ginput(1);
 cropx1= floor(cropx1);
 cropy1  = floor(cropy1);
 hold on, plot([1 n], [cropy1 cropy1], '-r');
 hold on, plot([cropx1 cropx1], [1 m], '-r');
 [cropx2 cropy2 ] = ginput(1);
 cropx2 = floor(cropx2);
 cropy2 = floor(cropy2);
 hold on, plot([1 n], [cropy2 cropy2], '-r');
 hold on, plot([cropx2 cropx2], [1 m], '-r');
 text(10,10, 'select ROI: upper left then lower right', 'Color', 'black');



end

