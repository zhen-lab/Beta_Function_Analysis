%%
close all;


%%
d = dorsal_smd./dorsal_smd_r;
v = ventral_smd./ventral_smd_r;
dn = normalize_signal(d);
vn = normalize_signal(v);
cnd = normalize_signal(curvdataBody);
cnv = normalize_signal(-curvdataBody);
bd = 35;
dns = dn(:,1); dne = dn(:,end);
vns = vn(:,1); vne = vn(:,end);
cnds = cnd(2:end,1); cnde = cnd(2:end,end);
cnvs = cnv(2:end,1); cnve = cnv(2:end,end);

figure; hold on; 
plot([bd bd],[0 1],'--k');
plot(dns,'r', 'linewidth', 4); 
plot(vns,'color', [0 0.5 1], 'linewidth', 4); 
plot(cnds, ':k', 'linewidth', 2); 
set(gca, 'visible', 'off');

figure; hold on;
plot([bd bd],[0 1],'--k');
plot(dne,'r', 'linewidth', 4); 
plot(vne,'color', [0 0.5 1], 'linewidth', 4); 
plot(cnde, ':k', 'linewidth', 2); 
set(gca, 'visible', 'off');


%%
diffdc = sum(([dns(bd:end) dne(bd:end)]-[cnds(bd:end) cnde(bd:end)]).^2);
diffvc = sum(([vns(bd:end) vne(bd:end)]-[cnvs(bd:end) cnve(bd:end)]).^2);
disp([diffdc diffvc]');
