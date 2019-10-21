function video_from_array(array, filename, varargin)
% VIDEO_FROM_ARRAY(array, filename)
%
%   This will take an array and write it out to a video.
%
%  VIDEO_FROM_ARRAY(array, filename, frames)
%
%   This will take an array and write the specified frames to a video.

S = size(array);

size_T = S(end);

default_options = struct(...
    'frames', 1:size_T, ...
    'framerate', 30 ...
);

input_options = varargin2struct(varargin{:}); 
options = mergestruct(default_options, input_options);


writer = VideoWriter(filename, 'MPEG-4');
writer.FrameRate = options.framerate;

open(writer);
for t = options.frames
    
    frame = get_slice(array, t);
    writeVideo(writer, frame);
    
end
close(writer);