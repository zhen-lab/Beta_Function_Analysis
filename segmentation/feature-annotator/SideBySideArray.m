classdef SideBySideArray
    % SIDEBYSIDEARRAY Lazily puts two arrays side by side.
    %
    % Author: Vivek Venkatachalam (vivekv2@gmail.com)

    properties

        % All arrays.
        Arrays

        % Dimension along which to concatenate arrays.
        CatDim = 2;

        % Additional spacing to place between the arrays.
        Spacing = 0;

        Size

        ElementClass

        % Always the last dimension
        SlicedDimension

        NSlices

    end

    methods

        function obj = SideBySideArray(arrays, dim, spacing)

            obj.Arrays = arrays;
            if nargin > 1
                obj.CatDim = dim;
            end
            if nargin > 2
                obj.Spacing = spacing;
            end

            sample_slice = get_slice(obj, 1);
            preimage_size = size(arrays{1});

            obj.NSlices = preimage_size(end);

            obj.Size = [size(sample_slice) obj.NSlices];
            obj.SlicedDimension = length(obj.Size);
            obj.ElementClass = element_class(sample_slice(1));

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

            s = size(get_slice(this.Arrays{1}, 1));
            s(this.CatDim) = this.Spacing;
            space = zeros(s, element_class(this.Arrays{1}));

            arrays_to_cat = {get_slice(this.Arrays{1}, t)};

            for i = 2:length(this.Arrays)
                arrays_to_cat{end+1} = space;
                arrays_to_cat{end+1} = get_slice(this.Arrays{i}, t);
            end

            data = cat(this.CatDim, arrays_to_cat{:});

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

    end


end

