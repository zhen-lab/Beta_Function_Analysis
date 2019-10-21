classdef ColorArray
    %COLORARRAY Creates an array with different color channels.
    %
    % Author: Vivek Venkatachalam (vivekv2@gmail.com)

    properties (SetAccess = immutable)

        % Set of color arrays
        Colors

        PreimageSize

        % Determines the color dimension. This is typically 3 for a 2D
        % array and 4 for a 3D array.
        CDim = 3;
        
        % Convert images to RGB or not
        ConvertToRGB = false;

        % Number of color channels
        SizeC

        Size
        
        % This will differ from Size if RGB conversion is on.
        OutputSize

        ElementClass

        NSlices
        
        SlicedDimension

    end

    properties

        % RGB values for the colors, with each column corresponding to a
        % separate color channel.
        LUT = eye(3);

    end

    methods

        function obj = ColorArray(colors, cdim, rgb, use_rgb)
            % obj = COLORARRAY(colors, cdim, rgb, use_rgb)
            %
            %   Creates an array from the color channels in colors.

            if ~iscell(colors)
                colors = {colors};
            end

            N = length(colors);
            
            obj.Colors = colors;
            obj.CDim = cdim;
            obj.SizeC = N;
            
            if nargin > 2
                obj.LUT = rgb;
            end
            obj.LUT = obj.LUT(:, 1:N);

            if nargin > 3
                obj.ConvertToRGB = use_rgb;
            end
            
            p = size(colors{1});
            
            image_size = [p(1:cdim-1) N p(cdim:end)];
                        
            obj.PreimageSize = p;
            obj.Size = image_size;
            obj.SlicedDimension = length(obj.Size);
            
            if obj.SlicedDimension == cdim
                obj.SlicedDimension = cdim + 1;
            end
            
            obj.ElementClass = element_class(colors{1});

            obj.OutputSize = image_size;
            if obj.ConvertToRGB
                obj.OutputSize(obj.CDim) = 3;
            end
            
        end

        function [varargout] = subsref(this, S)

            % Determine which slices we will need to transform
            requested = S.subs{this.SlicedDimension};

            data = zeros([this.OutputSize(1:end-1) length(requested)], ...
                this.ElementClass);

            idx = num2cell(repmat(':', 1, length(this.OutputSize)));
            for i = 1:length(requested)

                idx{end} = i;
                data(idx{:}) = this.get_slice(requested(i));

            end

            new_S = S;
            new_S.subs{this.SlicedDimension} = ':';
            varargout{1} = subsref(data, new_S);

        end
        
        function data = get_slice(this, t)
            
            if this.ConvertToRGB
                data = this.get_rgb_slice(t);
            else
                data = this.get_raw_slice(t);
            end
            
        end

        function data = get_raw_slice(this, t)

            assert(numel(t)==1, ...
                'get_slice can only be called on single slices');

            idx = num2cell(repmat(':', 1, this.SlicedDimension));
            idx{end} = 1;

            data_size = [this.Size(1:end-1) 1];

            data = zeros(data_size, this.ElementClass);

            slice_idx = idx;
            for c = 1:this.SizeC
                slice_idx{this.CDim} = c;
                data(slice_idx{:}) = get_slice(this.Colors{c}, t);
            end

        end
        
        function data = get_rgb_slice(this, t)
            
            cdata = get_raw_slice(this,t);
            data = multiply_nd(this.LUT, cdata, this.CDim);
            
        end
        
        function colors = get_colors(this)
            
            colors = this.Colors;
            
        end

        function s = size(obj, d)
            s = obj.OutputSize;
            
            if nargin > 1
                if d <= length(s)
                    s = s(d);
                else
                    s = 1;
                end
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

    methods (Static)


    end

end
