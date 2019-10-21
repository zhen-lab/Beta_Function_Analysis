classdef MatfileArray
    
    properties (SetAccess=protected)
        
        % matfile object that corresponds to where the data is located on
        % disk
        Matfile

        % Field in matfile this cached array will access
        Field
        
        Size

        Type
        
    end
    
    methods
        
        function obj = MatfileArray(filename, field)
            % x = MATFILEARRAY(filename, field)
            %
            %   Creates read-only reference to the array contained in a
            %   specified field of a .mat file.
            
            obj.Matfile = matfile(filename);
            obj.Field = field;
            obj.Size = size(obj.Matfile, field);
            
            % Determine the type and size of the array.
            mfile_info = whos(obj.Matfile);

            for i = 1:length(mfile_info)
                if strcmp(mfile_info(i).name, field)
                    array_info = mfile_info(i);     
                end
            end

            obj.Size = array_info.size;
            obj.Type = array_info.class;
            
        end
        
        function [varargout] = subsref(obj, S)

            varargout{1} = obj.Matfile.(obj.Field)(S.subs{:});

        end
        

        function s = size(this, d)
                
            s = this.Size;
            
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