function [offset, score, tform] = register_image(src, target, varargin)
% [offset, score, new_img] = REGISTER_IMAGE(src, target, kwargs)
%
%   Registers a PointFeature (parent) in an image with the specified 
%   options. This method uses image registration with maximum intensity
%   projections of a local volume around the feature's center.
%
%   ARGS:
%
%     'src': A 2D or 3D image (YXC), typically smaller than target.
%
%     'target': A 2D or 3D image (YXC), typically larger than src with src
%       contained somewhere within.
%
%   KWARGS:
%
%     'guess': An initial guess of the translation part of the
%       transformation. Origins are assumed to be the centers of the
%       images.
%
%     'guess_tform': A guess of the initial transformation.
%
%     'iterations': Number of times to randomly sample offsets while
%       searching for a match.
%
%     'threshold': A cutoff for a good enough match, if multiple iterations
%       are allowed.
%
% Author: Vivek Venkatachalam (vivekv2@gmail.com)

default_options = struct(...
    'guess', [0 0], ...
    'guess_tform', [], ...
    'plot', false, ...
    'iterations', 1, ...
    'threshold', 0.8, ...
    'clean', true, ...
    'pyramids', 1 ...
);

input_options = varargin2struct(varargin{:}); 
options = mergestruct(default_options, input_options);

S = size(src);
Rsrc = imref2d(S, ...
    [0, S(2)] - S(2)/2, ...
    [0, S(1)] - S(1)/2 ...
);

S = size(target);
Rtarget = imref2d(S, ...
    [0, S(2)] - S(2)/2, ...
    [0, S(1)] - S(1)/2 ...
);

affine2d_from_shift = @(x) affine2d([1,0,0; 0,1,0; x(2), x(1), 1]);

if isempty(options.guess_tform)
    G = options.guess;
    guess_tform = affine2d_from_shift(G);
else
    guess_tform = options.guess_tform;
end

[optimizer, metric] = imregconfig('monomodal');
%optimizer.MaximumStepLength = 1e-4;
%optimizer.RelaxationFactor = 0.1;
get_tform = @(src, target, T) ...
    imregtform(src, Rsrc, target, Rtarget, 'translation', optimizer, metric, ...
        'InitialTransform', T, ...
        'PyramidLevels', options.pyramids);

N = 0;
score = 0;
done = false;

offset_to_center = size(target)/2 + 0.5;

guessed_centers = zeros(1, 2) + offset_to_center;

best_offset = [0, 0];
best_score = 0;
best_tform = guess_tform;
best_target_masked = [];

while N < options.iterations
    
    try
        tform = get_tform(src, target, guess_tform);

        offset = [tform.T(3,2) tform.T(3,1)];

        new_img = imwarp(src + 1, Rsrc, tform, 'OutputView', Rtarget);
        target_masked = target .* cast(new_img > 0, class(target));
        new_img = new_img - 1;
        score = get_image_overlap(new_img, target_masked);
        
        if score > best_score
            best_offset = offset;
            best_score = score;
            best_tform = guess_tform;
            best_new_img = new_img;
            best_target_masked = target_masked;
        end
        
    catch e
        
        warning(e.message)
        
    end

    if score > options.threshold
        if done
            break
        else
            guess_tform = tform;
            done = true;
        end
    else
        guess_idx = sample_and_avoid(size(target), guessed_centers, ...
            size(target)*0.1, size(target)*0.9);
        guess_idx = cell2mat(guess_idx);
        guessed_centers = [guessed_centers; guess_idx];
        guess_tform = affine2d_from_shift(guess_idx - offset_to_center);
    end

    N = N + 1;
end

score = best_score;
offset = best_offset;
tform = best_tform;

if options.plot && ~isempty(best_target_masked)

    figure(1); clf; imshowpair(best_new_img, best_target_masked);
    pause();

end

