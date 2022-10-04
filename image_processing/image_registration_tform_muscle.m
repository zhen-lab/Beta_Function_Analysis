function [img_updated, movingRegistered, tform] = image_registration_tform_muscle(img_tofix, img_tomove, lr)

answer = inputdlg({'Registration coefficient','Frame number to be used'},...
    'Input frame number', [1 35], {'1','1'});
coef = str2double(answer{1});
regfrmnum = str2double(answer{2});

fixed = img_tofix{regfrmnum};
moved = img_tomove{regfrmnum};
% movingRegistered = img_tomove;
% img_updated = cell(size(img_tofix),1);

[optimizer, metric] = imregconfig('multimodal');
optimizer.InitialRadius = optimizer.InitialRadius/coef;
tform = imregtform(moved, fixed, 'rigid', optimizer, metric);

% for i = 1:size(img_tofix, 1)
%     
%     movingRegistered{i} = imwarp(img_tomove{i},tform,'OutputView',imref2d(size(img_tomove{i})));
%     img_updated{i,1} = [img_tofix{i} movingRegistered{i}];
%     
% end

movingRegistered = cellfun(@(x) imwarp(x,tform,'OutputView',imref2d(size(x))), img_tomove, 'uniformoutput', 0);
switch lr
    case 'Left'
        img_updated = cellfun(@(x,y) [x y], movingRegistered, img_tofix, 'uniformoutput', 0);
    case 'Right'
        img_updated = cellfun(@(x,y) [x y], img_tofix, movingRegistered, 'uniformoutput', 0);
end

end
