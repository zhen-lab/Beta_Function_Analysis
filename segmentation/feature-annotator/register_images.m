function [rotated, tform] = register_images(moving, fixed, varargin)
% REGISTER_IMAGES(moving, fixed)
%
%   Registers two images.
%
% Author: Vivek Venkatachalam (vivekv2@gmail.com)

default_options = struct(...
    'fignum', 1, ...
    'guess', affine2d(eye(3)), ...
    'correlation_target', NaN ...
);

input_options = varargin2struct(varargin{:}); 
options = mergestruct(default_options, input_options);


[optimizer, metric] = imregconfig('monomodal');
tform = imregtform(moving, fixed, 'rigid', optimizer, metric, ...
    'InitialTransformation', options.guess);

tform = clean_rigid_transform(tform);

rotate = @(x, t) imwarp(x, t, 'OutputView', imref2d(size(fixed)));

rotated = rotate(moving, tform);
figure(options.fignum); subplot(121); imshowpair(moving, fixed);
figure(options.fignum); subplot(122); imshowpair(rotated, fixed);

if ~isnan(options.correlation_target)

    best_fit = get_image_overlap(rotated, fixed);
    best_tform = tform;

    for i = 1:100

        if best_fit < options.correlation_target

            new_angle = randi([0, 359]);
            guess = rotate_tform(options.guess, new_angle, size(moving)/2);

            tform = imregtform(moving, fixed, 'rigid', optimizer, ...
                metric, 'InitialTransformation', guess);

            rotated = rotate(moving, tform);

            fit = get_image_overlap(rotated, fixed);
            if fit > best_fit
                best_fit = fit;
                best_tform = tform;
            end

            figure(options.fignum); subplot(121); 
            imshowpair(moving, fixed);

            figure(options.fignum); subplot(122); 
            imshowpair(rotated, fixed);

        else

            tform = best_tform;
            rotated = rotate(moving, best_tform);

            return
        end
    end

    warning(['Image registration failed. Target correlation: %d. '...
        'Actual: %d'], options.correlation_target, best_fit);
    tform = best_tform;
    rotated = rotate(moving, best_tform);

end