%% Clear everything

clear; close all;

%% Load images

% Load tiff files
setup_proof_reading;
fprintf('tiff files loading finished. \n');

%% Register channels without existing geometric transformation object

close all;
figure; set(gcf, 'units', 'normalized',...
    'outerposition', [0.1 0.1 0.8 0.8]);
subplot(1,2,1); 
imshowpair(imagelist_g{1,1}, imagelist_r{1,1});
image_registration_tform; 
subplot(1,2,2); 
imshowpair(movingRegistered{1,1}, imagelist_r{1,1});
title('Press return to continue...');
pause;
close all;

save([filename(1:end-4) '_tform.mat'],...
    'tform', 'optimizer', 'metric');

fprintf('registration completed, transformation object saved. \n');

%% Register channels with existing geometric transformation object

transform = uigetfile('*tform.mat', ...
    'Select geometric transformation object');
load(transform, 'tform');

movingRegistered = imagelist_g;
imagelistRegistered = imagelist;
for i = 1:size(imagelist_r, 1)
    
    movingRegistered{i} = imwarp(imagelist_g{i},tform,'OutputView',imref2d(size(imagelist_g{i})));
    imagelistRegistered{i,1} = [imagelist_r{i} movingRegistered{i}];
    
end

fprintf('registration completed. \n');

%% Track neuron and generate signal
% Navigate to data folder if loading tracked neurons is necessary!!!!!

if size(imagelist_g) == size(imagelist_r)
    frmrange = [1 size(imagelist_g,1)];
end
imagelist = imagelistRegistered;
proof_reading(imagelist, [], filename, ...
    frmrange(1), frmrange(2), 1, framestruct);

% %% Save data
% 
% prompt = {'Neuron analyzed just now: '};
% dlgtitle = ['Neuron name for ' filename];
% dims = [1 35];
% definput = {''};
% answer = inputdlg(prompt,dlgtitle,dims,definput);
% neuronname = answer{1};
% if isempty(neuronname)
%     data_path_end = [filename(1:end-4) '.mat'];
% else
%     data_path_end = [neuronname '_' filename(1:end-4) '.mat'];
% end
% 
% parts = strsplit(pathname, '\');
% data_path = fullfile(parts{1,1:end-3}, 'Alpha_Data_Raw', parts{1,end-1});
% warning('off'); mkdir(data_path);
% 
% data_path_name = fullfile(data_path, data_path_end);
% save(data_path_name, ...
%     'signal', 'signal_mirror', 'ratio', 'neuron_position_data', 'dual_position_data');
% 
% fprintf([data_path_end ' data saved. \n']);

%% Curvature analysis with three points

generate_curvature_from_pts;

%% Curate frames if necessary

p3 = [1006:1019 1177:1196];
p4 = [1:114 211:276 376:484 711:865 1281:1800];
p5 = [1:168 193:277 371:484 709:865 1212:1800];
deleted = union_several(p3, p4, p5);
filename = 'temp19.tif';

% Save data
parts = strsplit(pwd, '\');
data_path = fullfile(parts{1,1:end-3}, 'Alpha_Data_Raw', parts{1,end-1});
warning('off'); mkdir(data_path); 
data_path_name = fullfile(data_path, [filename(1:end-4) '_deleted_frames.mat']);
save(data_path_name, 'deleted');
fprintf('deleted frames info saved. \n');

%% Curvature analysis with contours

% Save images with contours

neuronpos_1 = uigetfile('.mat', 'Select neuron position file');
load(neuronpos_1, 'neuron_position_data'); pos_1 = neuron_position_data;
if max(pos_1(:,1))>size(img_stack,2)/2
    pos_1(:,1) = pos_1(:,1)-size(img_stack,2)/2;
end
neuronpos_2 = uigetfile('.mat', 'Select neuron position file');
load(neuronpos_2, 'neuron_position_data'); pos_2 = neuron_position_data;
if max(pos_2(:,1))>size(img_stack,2)/2
    pos_2(:,1) = pos_2(:,1)-size(img_stack,2)/2;
end
yl = size(img_stack,1);
close all;

for i = 1:size(imagelist,1)
    
    hold off;
    imshow(imagelist_g{i,1}); caxis([0 2000]); hold on;
    plot(dorsal_data{2*i-1,1}, dorsal_data{2*i,1}, 'r');
    plot(flipud(ventral_data{2*i-1,1}), flipud(ventral_data{2*i,1}), 'b');
    plot(centerline_data_spline(:,2*i-1), centerline_data_spline(:,2*i), 'w');
    plot(pos_1(i,1), pos_1(i,2), 'om', 'markersize', 10);
    plot(pos_2(i,1), pos_2(i,2), 'om', 'markersize', 10);
    text(10, 10, num2str(i), 'color', 'w');
    
    drawnow;
    frame = getframe(1);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);
    if i == 1
      imwrite(imind, cm, [filename '_contour_overlaid.gif'], 'gif', 'LoopCount', inf, 'DelayTime', ff/fps);
    else
      imwrite(imind, cm, [filename '_contour_overlaid.gif'], 'gif', 'WriteMode', 'append', 'DelayTime', ff/fps);
    end
    
end

%% Write registered tiff file

for i = 1:length(imagelist)
    if i==1
        imwrite(imagelist{i,1}, [filename(1:end-4) '_reg.tif']);
    else
        imwrite(imagelist{i,1}, [filename(1:end-4) '_reg.tif'], 'writemode', 'append');
    end
end

fprintf('registered recording saved. \n');
