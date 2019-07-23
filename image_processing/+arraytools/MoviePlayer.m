classdef MoviePlayer < hgsetget
    
    properties
        Data
        T
        LUT = [0 500]
        Size_T
    end
    
    properties (Access=protected)
        Figure
        Slider
        Axes
    end
    
    methods
        
        function obj = MoviePlayer(x)
            % obj = MOVIEPLAYER(x)
            %
            %   Creates a movie player for array x.
            
            obj.Size_T = size(x, 3);
            
            obj.Figure = figure('Visible','off',...
                'Position',[100,100,500,550]);
            
            obj.Slider = uicontrol('Style','slider',...
                'Units', 'normalized', ...
                'Position',[0.05, 0.05, 0.9, 0.04], ...
                'Min', 1, ...
                'Max', obj.Size_T, ...
                'Value', 1, ...
                'SliderStep', [1/obj.Size_T 10/obj.Size_T], ...
                'Callback', @(h,v) obj.slider_callback(h) ...
            );
            
        
            obj.Axes = axes(...
                'Position', [.05, .15, .9, .8], ...
                'parent', obj.Figure ...
            );
        
            set(obj.Figure, 'visible', 'on');
            
            obj.Data = x;
            
            obj.T = 1;
            
        end
        
        function this = set.T(this, val)
            
            imshow(this.Data(:,:,val), this.LUT, 'Parent', this.Axes);
            this.T = val;
            
            set(this.Slider, 'Value', val);
            
        end
        
    end
    
    methods (Access=protected)
        
        function slider_callback(this, slider)
            
            this.T = round(get(slider, 'Value'));
            
        end
        
    end
    
end