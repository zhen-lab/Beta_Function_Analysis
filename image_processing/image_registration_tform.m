answer = inputdlg({'Registration coefficient',...
    'Frame number to be used'}, 'Input frame number',...
    [1 35], {'1','1'});
coef = str2double(answer{1});
frmnum = str2double(answer{2});
fixed = imagelist_r{frmnum};
moving = imagelist_g{frmnum};
movingRegistered = imagelist_g;
imagelistRegistered = imagelist;

% imshowpair(fixed, moving);
[optimizer, metric] = imregconfig('multimodal');
optimizer.InitialRadius = optimizer.InitialRadius/coef;

tform = imregtform(moving, fixed, 'rigid', optimizer, metric);

for i = 1:size(imagelist_r, 1)
    
    movingRegistered{i} = imwarp(imagelist_g{i},tform,'OutputView',imref2d(size(imagelist_g{i})));
    imagelistRegistered{i,1} = [imagelist_r{i} movingRegistered{i}];
    
end