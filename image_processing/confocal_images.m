function [] = confocal_images(filename, frames_t)

frames_z = 20;
M = matfile(filename);

% Initialize the first frame
imwrite(M.images(:,:,1,1), [filename(1:end-4) '_' num2str(frames_t(1)) '_' num2str(frames_t(2)) '.tif']);

% Append the remaining frames
for j = frames_t(1):frames_t(2)
    for i = 1:frames_z      
        if ~(i==1 && j==frames_t(1))        
            imwrite(M.images(:,:,i,j), [filename(1:end-4) '_' num2str(frames_t(1)) '_' num2str(frames_t(2)) '.tif'], 'writemode', 'append');        
        end        
    end    
end

end