function feature = register_feature_frame(feature, image_location, t, varargin)
% feature = register_feature_frame(feature, image_location, t, 't_ref', t_ref)
%
%   Register a single frame of a given feature using the coordinates from
%   t_ref as guess. The image is then registered using image correlation
%   of max intensity projections, followed by searching for a local
%   intensity maximum.
%
% feature = register_feature_frame(feature, image_location, t, 'image_matfile', mfilename)
%
%   Obtains images from the provided matfile.

if feature.is_bad_frame(t)
    return;
end

default_options = struct(...
    'centroid', false, ...
    'filter', [], ...
    'feature_annotator_filter', false, ...
    'drift', [3, 3, 1], ...
    'imreg_padding', [3, 3, 1], ...
    'steps', 1, ...
    'fine_tune', false, ...
    'reference_features', [], ...
    'lookup_threshold', 0.1, ...
    't_ref', NaN, ...
    'parents', [], ...
    'image_matfile', [], ...
    'headless', false, ...
    'update_registration', false...
);

input_options = varargin2struct(varargin{:}); 
options = mergestruct(default_options, input_options);

refs = options.reference_features;

size_T = get_size_T(image_location);

im = load_image(image_location, 't', t);
N = ndims(im);

offset_to_center = 0.5*feature.size;
offset_to_large_center = 0.5*(feature.size) + options.drift;

if ~isempty(options.image_matfile)
    mfile = matfile(options.image_matfile.name);
    load_im0 = @(t) mfile.(options.image_matfile.channel)(:,:,:,t);
else
    load_im0 = @(t) load_image(image_location, 't', t);
end

if options.feature_annotator_filter

    load_im = @(t) feature_annotator_filter(load_im0(t));

elseif ~isempty(options.filter)

    load_im = @(t) imfilter(load_im0(t), options.filter);

else

    load_im = load_im0;

end

get_local_vol = @(t, center) ...
    get_image_section(center-offset_to_center, ...
        feature.size, ...
        load_im(t));

get_MIPs = @(t, center) ...
    get_expanded_MIPs(load_im(t), center, feature.size, [0 0 0]);

get_large_MIPs = @(t, center) ...
    get_expanded_MIPs(load_im(t), center, feature.size, ...
        options.imreg_padding);


near_features = [];
good_refs = [];
good_ref_coords =[];
t_ref = [];

if options.fine_tune && isnan(feature.modified_coordinates(t, 1))
    feature.registration_info{t}.guess_source = 'fine tuning';

    t_ref = t;
    guess = get_feature_center(feature, t);

elseif ~isnan(feature.modified_coordinates(t, 1))
    feature.registration_info{t}.guess_source = 'modified frame';

    t_ref = t;
    guess = get_modified_center(feature, t);

elseif ~isnan(options.t_ref)
    feature.registration_info{t}.guess_source = 'reference frame';

    t_ref = options.t_ref;
    guess = get_feature_center(feature, t_ref);

elseif ~isempty(refs)
    feature.registration_info{t}.guess_source = 'reference features';

    for j = 1:length(refs)

        if is_good_frame(refs{j}, t)

            good_refs(end+1) = j;
            good_ref_coords(:, end+1) = get_feature_center(refs{j}, t)';

        end

    end

    near_features = refs(good_refs);

    [lookup_guess, fit, t_ref] = lookup_feature_location(...
        feature, ...
        near_features, ...
        t);

    if ~options.headless
        disp(['Got fit : ' num2str(fit) ' from frame ' ...
              num2str(t_ref) ' for frame ' num2str(t)]); 
    end

    if t == t_ref

        guess = get_modified_center(feature, t);

    elseif fit > options.lookup_threshold

        guess = lookup_guess;

        if ~options.headless
            disp('Looked up answer!');
        end

    else

        % Find the reference coordinates from t_ref
        source = [];
        for j = 1:length(near_features)

            source(:,j) = get_feature_center(near_features{j}, t_ref)';

        end

        % Take the average displacement of the references
        old_location = get_feature_center(feature, t_ref);
        displacements = good_ref_coords - source;
        guess = old_location + mean(displacements, 2)';

    end

else
    feature.registration_info{t}.guess_source = 'previous time (default)';

    if t > 2
        t_ref = t - 1;
    else
        t_ref = 1;
    end

    guess = get_feature_center(feature, t_ref);

end

feature.registration_info{t}.initial_guess = guess;
feature.registration_info{t}.t_ref = t_ref;

if options.centroid

    vol = get_local_vol(t, guess);
    new_center = centroid(vol)-1;

    if ~isnan(new_center)

        guess = guess + (new_center - offset_to_center);

    end

else

    if ~all(isnan(get_modified_center(feature, t_ref)))
        % This is necessary to avoid accidentally a miscomputed frame when
        % a feature is modified and marked as registered.
        prev_center = get_modified_center(feature, t_ref);
    else
        prev_center = get_feature_center(feature, t_ref);
    end

    [px, py, pz] = get_MIPs(t_ref, prev_center);
    [cx, cy, cz] = get_large_MIPs(t, guess);

    offset_from_x = register_image(px, cx, 'iterations', 20, ...
        'threshold', 0.8, 'pyramids', 1, options);
    offset_from_y = register_image(py, cy, 'iterations', 20, ...
        'threshold', 0.8, 'pyramids', 1, options);
    offset_from_z = register_image(pz, cz, 'iterations', 20, ...
        'threshold', 0.8, 'pyramids', 2, options);
    
    new_x = guess(2) + offset_from_z(2);
    new_y = guess(1) + offset_from_z(1);
    new_z = guess(3) + mean([offset_from_x(2), offset_from_y(1)]);
    
    guess = [new_y, new_x, new_z];
    
    [~, ~, nz] = get_MIPs(t, guess);
    score = get_image_overlap(pz, nz);

    if ~isfield(feature.registration_info{t}, ...
                'image_registration_score') || ...
       score > feature.registration_info{t}.image_registration_score
       
        feature.registration_info{t}.image_registration_score = score;
        feature.registration_info{t}.image_registration_guess = guess;
        feature.registration_info{t}.image_registration_parent = t_ref;
        
    end

    for step = 1:options.steps

        % make a circular gaussian mask to penalize large drifts
        if ~all(isnan(options.drift))

            mask = ones(size(im));

            XN = cell(1,N);
            for j=1:N
                gv{j} = 1:size(mask,j);
            end
            [XN{:}] = ndgrid(gv{:});
            for j=1:N
                mask = mask.*exp(-((XN{j} - guess(j)).^2/ ...
                                  (2*options.drift(j)^2)));
            end 
            im = double(im).*mask;
        end

        % use a second gaussian mask to avoid nearby reference
        % features.
        if ~isempty(near_features)
            mask = ones(size(im));

            XN = cell(1,N);
            for j=1:N
                gv{j} = 1:size(mask,j);
            end
            [XN{:}] = ndgrid(gv{:});

            for k=1:length(near_features)

                % Find the distance to the reference feature.  That
                % will determine if need to include it, and also will
                % set scale for excluded region's width.
                ref = get_feature_center(near_features{k}, t);

                d = norm(ref - guess);

                feature_mask = ones(size(im));

                scale = 1.5;
                height = 5;

                s = options.drift;
                for j = 1:N

                    feature_mask = height/scale * feature_mask .* ...
                                        exp(-((XN{j} - ref(j)).^2/ ...
                                                ((scale*s(j))^2)));
                end

                mask = mask - feature_mask;
            end

            im = double(im).*mask;
        end

        % choose the brightest point

        guess = subpixel_max_ND(im);
    end

end

feature.registration_info{t}.local_maximum_guess = guess;

new_feature_location = guess - offset_to_center;
feature.coordinates(t,:) = new_feature_location;

if options.update_registration
    feature.is_registered(t) = true;
end

