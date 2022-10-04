function [] = plot_heatmap_perievent(ratio, filename, r)

opengl software;

f = figure(1);
cl = 1.12; cw = 0.4;
subplot(2,1,1);
imagesc(ratio'); colormap(magma);
c1 = colorbar; set(c1, 'location', 'eastoutside', 'units', 'normalized');
ylim([0 1.1*max(max(ratio))]);
pos = get(c1, 'position');
pos(1) = cl*pos(1);
pos(3) = cw*pos(3);
set(c1, 'position', pos);
% title('Heatmap of DB activities', 'fontsize', 15);
% set(gca, 'YDir', 'normal');
% set(gca, 'ytick', 1:size(ratio, 2)); % To display the head-tail direction from bottom to top
% set(gca, 'yticklabel', {'2' '3' '4' '5' '6' '7' '8/9'}); 
% fps = 17.22; 
% numframes = size(ratio, 1);
% ylabel('DA neuron number');  
% set(gca, 'xtick', [], 'ticklength', [0.001 0.001], 'linewidth', 2.5); 
% box on;
% xlim([1 numframes+1]);
axis tight; 
set(gca, 'visible', 'off'); set(gcf, 'color', 'w');

subplot(2,1,2); 
%fill([1 1 958 958], [0 10 10 0], [0.8 0.8 0.8], 'edgecolor', 'none');
hold on;
% for DA and DB
% coef = 2; 
% cmp = lbmap(coef*size(ratio, 2), 'blue');
% for DD
coef = 4; 
cmp = lbmap(coef*size(ratio, 2), 'blue');
c2 = colorbar; set(c2, 'location', 'eastoutside', 'units', 'normalized');
ylim([0 1.1*max(max(ratio))]);
pos = get(c2, 'position');
pos(1) = cl*pos(1);
pos(3) = cw*pos(3);
set(c2, 'position', pos, 'visible', 'off');
lm = ylim;

addcurv = questdlg('Add curvature data?', 'Curvdata', 'Yes', 'No', 'No');
if isequal(addcurv, 'Yes')
    
    for i=1:size(r,1)
        fill([r(i,1) r(i,2) r(i,2) r(i,1)],...
            [lm(1) lm(1) lm(2) lm(2)], 0.8*[1 1 1], 'edgecolor', 'none');
    end
    for i=1:size(ratio,2)
        curv = uigetfile('*_curvature_curated.mat', 'Select curvature data');
        load(curv, 'curvdatafiltered');
        yyaxis left; 
        plot(ratio(:, i), '-', 'color', cmp(2*i, :), 'linewidth', 2);   
        set(gca, 'yticklabel', [], 'xticklabel', [], 'ticklength', [0 0]);
        yyaxis right; 
        plot(curvdatafiltered/100, ':', 'color', cmp(2*i, :), 'linewidth', 2);
        yticks(-2:1:2);
        xlim([0 size(ratio,1)]); ylim([-2 2]);
    end
    set(gcf, 'color', 'w');
    hold off;
else
    
    for i=1:size(r,1)
        fill([r(i,1) r(i,2) r(i,2) r(i,1)],...
            [lm(1) lm(1) lm(2) lm(2)], 0.8*[1 1 1], 'edgecolor', 'none');
    end
    for i=1:size(ratio, 2)
        % for DA and DB
    %     plot(ratio(:, i), 'color', cmp(i+coef, :), 'linewidth', 2);
        % for DD
        plot(ratio(:, i), 'color', cmp(i, :), 'linewidth', 2);
    end
    axis tight;
    set(gca, 'visible', 'off'); set(gcf, 'color', 'w');
    hold off;
end
           
% title('Peri-event plot of DB activities', 'fontsize', 15);
% xlim([1 size(ratio, 1)+1]); ylim([0 1.1*max(max(ratio))]);
% set(gca, 'xtick', [], 'ticklength', [0.001 0.001], 'linewidth', 2.5); 
% set(gca,'XTICK',1:10*fps:numframes);
% x_tick=get(gca,'XTICK');
% set(gca,'XTICKLABEL',(x_tick-1)/fps, 'ticklength', [0.001 0.001], 'linewidth', 2.5, 'layer', 'top');xlabel('Time/s'); ylabel('Signal');
% set(gcf, 'color', [1 1 1 ]);


savefig(f, [filename '.fig']);

end
