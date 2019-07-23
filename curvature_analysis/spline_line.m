function [smooth_line, cv2, splen] = spline_line(raw_line, spline_p, numcurvpts)
    
    xy = circshift(raw_line, [0 1])'; % Transpose, and x coordinates on the first row, y on the second row (raw_line in [y x])
    df = diff(xy,1,2); % 2nd-1st, 3rd-2nd, ...
    t = cumsum([0, sqrt([1 1]*(df.*df))]); % Cumulative sum of distances
    cv = csaps(t,xy,spline_p);

    cv2 =  fnval(cv, t)';
    df2 = diff(cv2,1,1); df2p = df2';

    splen = cumsum([0, sqrt([1 1]*(df2p.*df2p))]);
    
    smooth_line = interp1(splen+.00001*(0:length(splen)-1), cv2, 0:(splen(end)-1)/(numcurvpts-1):(splen(end)-1));

end

