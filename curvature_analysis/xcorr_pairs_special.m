function [correall, lagsall, lagtimeall] = xcorr_pairs_special(ratio, curv, frphase, prom)

curvnum = size(curv,2);

if curvnum>=3 
    
    ratiopruned = ratio(frphase,:);
    curvpruned = curv(frphase,:);
    ratiosmd = smoothdata(ratiopruned, 'rloess');
    ratiosmdnorm = normalize_signal(ratiosmd);
    curvsmdnorm = normalize_signal(curvpruned); % curv data is already smoothened
    compnum = 3; % compare three pairs: n-1, n, n+1
    cmp = gray(compnum*10);
    % curv number should be 2 more than ratio number used for complete comparisons
    correall = zeros(2*length(frphase)-1,compnum,curvnum-2);
    lagsall = zeros(compnum,2*length(frphase)-1,curvnum-2);
    lagtimeall = zeros(compnum,curvnum-2);

    if isequal(size(ratio),size(curv)) % this is the situation for DAs
        figure; hold on;
        for nm = 1:curvnum-2
            for idx = 1:compnum
                [corre, lags] = xcorr(ratiosmdnorm(:,nm+1), curvsmdnorm(:,nm+idx-1), 'normalized');
                correall(:,idx,nm) = corre;
                lagsall(idx,:,nm) = lags;
                % Obtain the local max nearest to point zero
                TF = islocalmax(corre, 'minprominence', prom);
                if any(TF)
                    [~, I] = min(abs(lags(TF)));
                    maxidx = find(TF,I); 
                    disp(num2str(maxidx));
                    lagtime = lags(maxidx(end)); 
                    fprintf(['local max closest to 0 is: ' num2str(maxidx(end)) '\n']);
                    lagtimeall(idx,nm) = lagtime;
                    plot(lags, corre, 'color', cmp(end-idx*8,:));
                    plot(lagtime, corre(lagtime+length(frphase)), 'or');
                    fprintf(['comparison #: ' num2str(idx) ' completed.\n']);
                else
                    fprintf('no local max is found at minprominence.\n');
                end
            end
        end
        hold off;
        disp(num2str(lagtimeall));
    elseif  (size(ratio,2)==1&&curvnum==3) % this is the situation for DBs
        figure; hold on;
        for idx = 1:compnum
            [corre, lags] = xcorr(ratiosmdnorm, curvsmdnorm(:,idx), 'normalized');
            correall(:,idx,:) = corre;
            lagsall(idx,:,:) = lags;
            % Obtain the local max nearest to point zero
            TF = islocalmax(corre, 'minprominence', prom);
            if any(TF)
                [~, I] = min(abs(lags(TF)));
                maxidx = find(TF,I); 
                disp(num2str(maxidx));
                lagtime = lags(maxidx(end)); 
                fprintf(['local max closest to 0 is: ' num2str(maxidx(end)) '\n']);
                lagtimeall(idx,:) = lagtime;
                plot(lags, corre, 'color', cmp(end-idx*8,:));
                plot(lagtime, corre(lagtime+length(frphase)), 'or');
                fprintf(['comparison #: ' num2str(idx) ' completed.\n']);
            else
                fprintf('no local max is found at minprominence.\n');
            end
        end
        hold off;
        disp(num2str(lagtimeall));
    else
        fprintf('ratio and curv data mismatch.\n');
    end
    
else
    fprintf('curv data has too few columns for complete comparisons.\n');
end

end