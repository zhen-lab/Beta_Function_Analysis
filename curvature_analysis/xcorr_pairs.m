function [r, lags, lagtimeall] = xcorr_pairs(ratio, curv, frphase)

ratiopruned = ratio(frphase,:);
curvpruned = curv(frphase,:);
ratiosmd = smoothdata(ratiopruned, 'rloess');
ratiosmdnorm = normalize_signal(ratiosmd);
curvsmdnorm = normalize_signal(curvpruned); % it is already smoothened

loops = size(ratio, 2);

for idx = 1:loops
    
    if idx<loops
        [rnn_1(:, idx), lags] = xcorr(ratiosmdnorm(:, idx+1), curvsmdnorm(:, idx), 'normalized');
        lagtime_1(idx) = lags(rnn_1(:, idx)==max(rnn_1(:, idx)));
    end
    
        [rnn(:, idx), lags] = xcorr(ratiosmdnorm(:, idx), curvsmdnorm(:, idx), 'normalized');
        lagtime_2(idx) = lags(rnn(:, idx)==max(rnn(:, idx)));
                
    if idx>1
        [rnn__1(:, idx-1), lags] = xcorr(ratiosmdnorm(:, idx-1), curvsmdnorm(:, idx), 'normalized');
        lagtime_3(idx-1) = lags(rnn__1(:, idx-1)==max(rnn__1(:, idx-1)));
    end
        
end

r = cell(3, 1);
r{1} = rnn_1; 
r{2} = rnn; 
r{3} = rnn__1;

lagtimeall = cell(3, 1);
lagtimeall{1} = lagtime_1';
lagtimeall{2} = lagtime_2';
lagtimeall{3} = lagtime_3';

disp(num2str([lagtimeall{1}; nan; lagtimeall{2}; nan; lagtimeall{3}; nan]));

end