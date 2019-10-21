function out = interp_nans(in, method)
% out = INTERP_NANS(in, method)
%
%   Interpolate NaN values.
%
%
% Author: Vivek Venkatachalam (vivekv2@gmail.com)

idx = 1:length(in);
good_idx = idx(~isnan(in));
good_in = in(good_idx);

out = in;

if ~isempty(good_idx)
    interp_range = good_idx(1):good_idx(end);
    out(interp_range) = interp1(good_idx, good_in, interp_range, method);
end
