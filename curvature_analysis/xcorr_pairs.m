function [r, lags, lagtimeall] = xcorr_pairs(ratio_rev, curv_rev)
 
ratio_rev_norm = ratio_rev;
curv_rev_norm = curv_rev;
loops = size(ratio_rev, 2);

% for i = 1:loops
%     
% ratio_rev_norm(:, i) = normalize_signal(ratio_rev(:, i));
% curv_rev_norm(:, i) = normalize_signal(curv_rev(:, i));
% 
% end

for i = 1:loops
    
    if i<loops
        [rnn_1(:, i), lags] = xcorr(ratio_rev_norm(:, i+1), curv_rev_norm(:, i), 'coeff');
        lagtime_1(i) = lags(rnn_1(:, i)==max(rnn_1(:, i)));
    end;
    
        [rnn(:, i), lags] = xcorr(ratio_rev_norm(:, i), curv_rev_norm(:, i), 'coeff');
        lagtime_2(i) = lags(rnn(:, i)==max(rnn(:, i)));
                
    if i>1
        [rnn__1(:, i-1), lags] = xcorr(ratio_rev_norm(:, i-1), curv_rev_norm(:, i), 'coeff');
        lagtime_3(i-1) = lags(rnn__1(:, i-1)==max(rnn__1(:, i-1)));
    end;
        
end

r = cell(3, 1);
r{1} = rnn_1; r{2} = rnn; r{3} = rnn__1;

lagtimeall = cell(3, 1);
lagtimeall{1} = lagtime_1';
lagtimeall{2} = lagtime_2';
lagtimeall{3} = lagtime_3';
disp(num2str([lagtime_1(1) lagtime_2(end-1) lagtime_3(end)]'));
end