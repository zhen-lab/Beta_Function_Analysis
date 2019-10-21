frmnum = 10;
fixed = imagelist_r{frmnum};
moving = imagelist_g{frmnum};
movingRegistered = imagelist_g;

% imshowpair(fixed, moving);
[optimizer, metric] = imregconfig('multimodal');
optimizer.InitialRadius = optimizer.InitialRadius/coef;

tform = imregtform(moving, fixed, 'rigid', optimizer, metric);

for i = 1:size(imagelist_r, 1)
    
    movingRegistered{i} = imwarp(imagelist_g{i},tform,'OutputView',imref2d(size(imagelist_g{i})));
    imagelist{i,1} = [imagelist_r{i} movingRegistered{i}];
    
end