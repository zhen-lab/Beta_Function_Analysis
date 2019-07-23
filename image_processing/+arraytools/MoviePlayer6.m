classdef MoviePlayer6 < hgsetget
    
    properties
        DataX
        DataY
        DataZ
        T
        LUT = [0 500]
        Size_T
    end
    
    properties (Access=protected)
        Figure
        Slider
        AxesX
        AxesY
        AxesZ
    end
    
    methods
        
        function obj = MoviePlayer6(x)
            % obj = MOVIEPLAYER(x)
            %
            %   Creates a movie player for array x.
            
            import arraytools.*;
            
            DataX = CachedArray(LambdaArray(MatfileArray(x, 'images'), ...
                @max_intensity_x_side_camera));
            DataY = CachedArray(LambdaArray(MatfileArray(x, 'images'), ...
                @max_intensity_y_side_camera));
            DataZ = CachedArray(LambdaArray(MatfileArray(x, 'images'), ...
                @max_intensity_z_side_camera));
            
            obj.Size_T = size(DataZ, 3);
            
            obj.Figure = figure('Visible','off',...
                'Position',[100,100,1100,1100]);
            
            obj.Slider = uicontrol('Style','slider',...
                'Units', 'normalized', ...
                'Position',[0.05, 0.05, 0.9, 0.04], ...
                'Min', 1, ...
                'Max', obj.Size_T, ...
                'Value', 1, ...
                'SliderStep', [1/obj.Size_T 10/obj.Size_T], ...
                'Callback', @(h,v) obj.slider_callback(h) ...
            );
            
        
            obj.AxesZ = axes(...
                'Position', [.05, .25, .8, .7], ...
                'parent', obj.Figure ...
            );
            obj.AxesX = axes(...
                'Position', [.85, .25, .1, .7], ...
                'parent', obj.Figure ...
            );
            obj.AxesY = axes(...
                'Position', [.05, .15, .8, .1], ...
                'parent', obj.Figure ...
            );
        
            set(obj.Figure, 'visible', 'on');
            
            obj.DataX = DataX;
            obj.DataY = DataY;
            obj.DataZ = DataZ;
            
            obj.T = 1;
            
        end
        
        function this = set.T(this, val)
            
            imshow(this.DataX(:,:,val), this.LUT, 'Parent', this.AxesX);
            imshow(this.DataY(:,:,val), this.LUT, 'Parent', this.AxesY);
            imshow(this.DataZ(:,:,val), this.LUT, 'Parent', this.AxesZ);
            this.T = val;
            
            title(num2str(this.T));
            
            set(this.Slider, 'Value', val);
            
        end
        
    end
    
    methods (Access=protected)
        
        function slider_callback(this, slider)
            
            this.T = round(get(slider, 'Value'));
            
        end
        
    end
    
end