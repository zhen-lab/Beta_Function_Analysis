classdef RangeScaledArray < LazyArray
    % RANGESCALEDARRAY Lazily transforms the output range and datatype of
    % an array.
    %
    % Author: Vivek Venkatachalam (vivekv2@gmail.com)

    properties (SetAccess = immutable)

        % This is the range of values to be placed in the full output
        % range.
        Range

        BitDepth;


    end

    methods

        function obj = RangeScaledArray(array, input_range, bitdepth)
            
            if nargin < 3
                
                range_cutoff = 0.75;
                pixel_count = 0.9;
            
                vals = get_sample_slices(array);

                high_cutoff = max(vals);
                low_cutoff = min(range_cutoff*high_cutoff, ...
                    quantile(vals, pixel_count));

                input_range = [low_cutoff, high_cutoff];
            end

            if nargin < 4
                bitdepth = 8;
            end

            bitrange = ceil(log2(double(range(input_range))));
            bitrange = max(bitrange, bitdepth);
            bit_offset = bitrange - bitdepth;

            obj@LazyArray(array, @RangeScaledArray.downsample, ...
                {input_range(1), bit_offset, bitdepth});

            obj.Range = input_range;
            obj.BitDepth = bitdepth;

        end
        

    end

    methods (Static)

        function y = downsample(x, cutoff, bit_offset, bitdepth)

            if bitdepth == 8
                y = uint8(bitshift(x-cutoff, -bit_offset));
            elseif bitdepth == 16
                y = uint16(bitshift(x-cutoff, -bit_offset));
            end

        end

    end


end

