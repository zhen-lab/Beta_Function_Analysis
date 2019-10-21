classdef ImwarpedArray
    % IMWARPEDARRAY Lazy array generated using the imwarp function with
    % a fixed output view.
    %
    % Author: Vivek Venkatachalam (vivekv2@gmail.com)

    properties (SetAccess = immutable)

        Preimage

        PreimageSize

        Size

        ElementClass

        % Always the last dimension
        SlicedDimension

        NSlices

        OutputView

    end

    properties

        Transforms;

        DefaultTransform = affine2d(eye(3));

    end

    methods

        function obj = ImwarpedArray(array, tforms, view)

            obj.Preimage = array;
            obj.PreimageSize = size(array);

            identity = affine2d(eye(3));
            S = size(array);
            NSlices = S(end);

            if nargin < 2
                tforms = [];
            end

            if nargin < 3
                sample_slice = get_slice(array, 1);
                sample_size = size(sample_slice);
                view = imref2d(sample_size(1:2));
            end

            obj.OutputView = view;

            if ~iscell(tforms)

                if isa(tforms, 'affine2d')
                    obj.DefaultTransform = tforms;
                end

                tforms = cell(1, NSlices);

            end

            if length(tforms) < NSlices 
                [tforms{end+1:NSlices}] = deal([]);
            end

            obj.Transforms = tforms;
            obj.NSlices = NSlices;

            sample_slice = obj.get_slice(1);
            obj.Size = [size(sample_slice) obj.NSlices];
            obj.SlicedDimension = length(obj.Size);
            obj.ElementClass = element_class(sample_slice);

        end

        function [varargout] = subsref(this, S)

            % Determine which slices we will need to transform
            requested = S.subs{this.SlicedDimension};

            data = zeros([this.Size(1:end-1) length(requested)], ...
                this.ElementClass);

            idx = num2cell(repmat(':', 1, length(this.Size)));
            for i = 1:length(requested)

                idx{end} = i;
                data(idx{:}) = this.get_slice(requested(i));

            end

            new_S = S;
            new_S.subs{this.SlicedDimension} = ':';
            varargout{1} = subsref(data, new_S);

        end

        function data = get_slice(this, t)

            if isa(this.Transforms{t}, 'affine2d')
                tf = this.Transforms{t};
            else
                tf = this.DefaultTransform;
            end

            x = get_slice(this.Preimage, t);

            data = imwarp(x, tf, 'OutputView', this.OutputView);

        end

        function tforms = get_tforms(this)
            tforms = this.Transforms;
        end

        function s = size(obj, d)
            s = obj.Size;
            if nargin > 1
                s = s(d);
            end
        end

        function n = numel(this)
            n = prod(this.Size);
        end

        function n = ndims(this)
            n = length(this.Size);
        end

        function t = element_class(this)
            t = this.ElementClass;
        end

        function array = get_preimage(this)
            array = this.Preimage;
        end

    end

    methods (Static)

        function obj = update_tform(obj, idx, new_tform)

            tforms = obj.Transforms;
            tforms{idx} = new_tform;
            obj = ImwarpedArray(obj.Preimage, tforms, obj.OutputView);

        end

        function obj = move_tforms(source, target)

            tforms = source.Transforms;

            source_size = size(source);
            target_size = size(target);

            output_view = imref2d(target_size(1:2));

            vscale = target_size(1)/source_size(1);
            hscale = target_size(2)/source_size(2);

            if vscale ~= hscale
                warning(['Attempting to clone transformation between ' ...
                    'arrays with different aspect ratios.']);
            end

            for i = 1:length(tforms)

                tforms{i}.T(3,1) = hscale * tforms{i}.T(3,1);
                tforms{i}.T(3,2) = vscale * tforms{i}.T(3,2);

            end

            obj = ImwarpedArray(target, tforms);

        end

    end

end

