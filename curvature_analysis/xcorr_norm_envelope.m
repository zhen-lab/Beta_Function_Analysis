%% No deletion for most recordings

deleted = [];

%% Delete frames if necessary

% temp6 of DD_Calcium
a = [685:712 1000:1006 1157:1168 1246:1259 1293:1305];
b = [3:15 619:700 1012:1028 1224:1238 1382:1386];
deleted = union(a,b);

%% Normalize against wave envelope

range = setdiff(1:(length(curvdatafiltered)-1), deleted);
r_gfp = signal{1,1}(range);
r_rfp = signal_mirror{1,1}(range);
cv = curvdatafiltered(range);
rt = smoothdata(r_gfp./r_rfp,...
    'rloess');
[rtu, rtl] = envelope(rt);
[cvu, cvl] = envelope(cv);

% rtcorrect = (rt-rtl);
% cvcorrect = (curvdatafiltered-cvl);
% rtcorrect = rt./(rtu-rtl);
% cvcorrect = cv./(cvu-cvl);
rtcorrect = (rt-rtl)./(rtu-rtl);
cvcorrect = (cv-cvl)./(cvu-cvl);

[r, lags] = xcorr(rtcorrect, cvcorrect);
disp(num2str(lags(r==max(r))));
figure;
subplot(411); 
yyaxis left; plot(rt); ylabel('activity');
yyaxis right; plot(cv); ylabel('curvature');
subplot(412);
yyaxis left; plot(rtcorrect); ylabel('activity');
yyaxis right; plot(cvcorrect); ylabel('curvature');
subplot(413);
plot(vel_ap_sign_smd_mean(range));
subplot(414); 
plot(lags, r);
