classdef CachedArray < handle

    properties (SetAccess=protected)

        % Original array. This has to support standard '()' indexing.
        BaseArray

        Size

        Type

        % Dimension that caching will occur along. This is currently always
        % the last dimension.
        SlicedDimension

        % Number of elements along the cached dimension
        NSlices

        % All data that has been accessed will be cached here
        LocalData

        % This bitvector keeps track of which slices have been cached. It
        % should match the size of the last dimension of the array
        IsCached

        % This bitvector keeps track of which slices have been modified and
        % not written to disk.
        IsModified

    end

    methods

        function obj = CachedArray(x)
            % y = CACHEDARRAY(x)
            %
            %   Read-only array that caches values in memory as they are
            %   read from disk.
            %
            %   Caching occurs along the last dimension, so if you request
            %   y(1,5,8:10), y(:,:,8:10) will be loaded into memory and
            %   y(1,5,8:10) will be returned.

            obj.BaseArray = x;

            % To determine the class of elements that are stored, we will
            % use 'class' for built-in arrays (which should almost never be
            % cached) and 'type' for others.
            try
                obj.Type = type(x);
            catch
                obj.Type = class(x);
            end
            
            obj.Size = size(obj.BaseArray);
            obj.SlicedDimension = length(obj.Size);
            obj.NSlices = obj.Size(obj.SlicedDimension);
            obj.LocalData = zeros(obj.Size, obj.Type);
            obj.IsCached = false(obj.NSlices, 1);
            obj.IsModified = false(obj.NSlices, 1);

        end

        function [varargout] = subsref(obj, S)

                assert(length(S) == 1, ...
                    'Cached arrays support only one level of indexing.');
                assert(strcmp(S(1).type, '()'), ...
                    'Cached arrays only support () indexing.');
                assert(length(S.subs)==length(obj.Size), ...
                    'Cached arrays require full indexing.');

                requested = false(obj.NSlices, 1);
                requested(S.subs{obj.SlicedDimension}) = true;

                needed = requested & ~obj.IsCached;

                % Convert logical indexing into numeric indexing.
                x = 1:length(needed);
                needed = x(needed);

                if any(needed)
                    idx = num2cell(repmat(':', 1, length(obj.Size)));
                    idx{obj.SlicedDimension} = needed;
                    obj.LocalData(idx{:}) = obj.BaseArray(idx{:});
                    obj.IsCached(needed) = true;
                end

                varargout{1} = subsref(obj.LocalData, S);

        end
        
        function s = size(obj, d)
            s = obj.Size;
            if nargin == 2
                s = s(d);
            end
        end
        
        function n = numel(obj)
            n = prod(obj.Size);
        end
        
        function n = ndims(this)
            n = length(this.Size);
        end
        
        function t = type(obj)
            t = obj.Type;
        end

    end

end