function tform = rotate_tform(tform, theta, center, varargin)
% tform = ROTATE_TFORM(tform, theta, center)
%
%   Rotates a transformation by a specified angle about a specified center.
%
% Author: Vivek Venkatachalam (vivekv2@gmail.com)

default_options = struct(...
);

input_options = varargin2struct(varargin{:}); 
options = mergestruct(default_options, input_options);

theta_0 = atan2(tform.T(2), tform.T(1));
a = theta;


R = [cos(a), -sin(a); sin(a), cos(a)];


b = center*(eye(2)-R);

tform2 = affine2d([R, [0;0]; b, 1]);

%tform = compose_transforms(tform, tform2);

tform = clean_rigid_transform(tform2);
