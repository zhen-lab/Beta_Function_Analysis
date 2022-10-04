
% button = length(questdlg('Load new data?','','Yes (TIF)','Yes (MAT) ','No', 'Yes (TIF)') ) ;
% 
% if button == 10
%     
%     [filename, pathname] = uigetfile({'*.mat'});
%     load([pathname filename]);
% 
% elseif button == 9
% 
%     do_dialog = 1;
% 
%     if do_dialog

        if exist('pathname', 'var')
            
            try
                
                if isdir(pathname)
                
                    cd(pathname);
                
                end
                
            end
            
        end

        [filename, pathname]  = uigetfile({'*.tif'});
        
        if ~(filename==0)

            fname = [pathname filename];

            %load the nd2 data

            if ~exist('data','var')

                data=bfopen(fname);

            end

            [num_series, num]=size(data);

            imagelist=data{num_series,1};

            [m,n]=size(imagelist{1,1});

            img_stack=zeros(m,n,length(imagelist));

            %if ~exist('calibration_matrix_green','var')

            %   disp('Please load green calibration matrix');

            %   [filename,pathname]  = uigetfile({'*.mat'});

            %   f_calibration_green=[pathname filename];

            %   load(f_calibration_green);
            % end

            for idx=1:size(img_stack,3)

                img_stack(:,:,idx)=imagelist{idx,1};

            end

    %         if ~exist('frames','var')
    %             
    %             frames = input('Please enter start and end frames for analyzing the data:','s');
    %             frames = str2num(frames);
    %         
    %         end
    % 
    %         istart = frames(1);
    %         iend = frames(2); 
            istart = 1; iend = idx;
            numframes=iend-istart+1; 

            %number of time series

    %     end
    
            framestruct = questdlg('Split view or alternating sampling?', 'Image', 'Split', 'Alternating', 'Split');
            switch framestruct
                case 'Split'
                    % Left and right screen
                    [imagelist_r, imagelist_g] = split_two_screens(imagelist);
                    range = [istart iend];
                case 'Alternating'
                    % Alternating frames
                    channelorder = questdlg('GFP or RFP channel first?',...
                        'Channel order',...
                        'GFP','RFP','GFP');
                    switch channelorder
                        case 'GFP'
                            imagelist_g = imagelist(1:2:end,1);
                            imagelist_r = imagelist(2:2:end,1);
                        case 'RFP'
                            imagelist_g = imagelist(2:2:end,1);
                            imagelist_r = imagelist(1:2:end,1);
                    end
                    range = [istart iend/2];
            end
%             [imagelist_r, imagelist_g] = split_two_screens(imagelist);
%             range = [istart iend];
        else
            fprintf('user canceled selection. \n');
        end
% end