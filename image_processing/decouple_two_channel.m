
imagelist_1 = imagelist(:,1);
imagelist_2 = imagelist(:,2);

imagelist_1 = reshape(imagelist_1, 2, []);
imagelist_2 = reshape(imagelist_2, 2, []);

imagelist_g = [imagelist_1(1,:)' imagelist_2(1,:)'];
imagelist_r = [imagelist_1(2,:)' imagelist_2(2,:)'];