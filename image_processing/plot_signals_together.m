thickness = 1;
fps = 17.22;
numframes = size(gfp, 1);
% period = [61 104; 408 423; 473 642];
% period = [148 350];
% period = [1 1094; 1602 1650; 1745 1800];
% period = [1 483; 561 1800];
% period = [1 58; 108 248; 264 572; 588 755; 782 815; 855 1055; 1091 1583; 1626 1771];
% period = [1 141; 180 268; 329 436; 543 1800];
% period = [45 373];
% period = [1 344];
% period = [1 678; 770 865];
% period = [1 592; 646 1178];
figure;
for i = 1: size(gfp, 2)
    
    subplot(size(gfp, 2), 1, (size(gfp, 2) - i) +1);
    plot(gfp(:, i), 'color', [0 0.8 0], 'linewidth', thickness); hold on;
    plot(rfp(:, i), 'color', [1 0.1 0], 'linewidth', thickness);
%     plot(curv_norm(:, i), 'k:');
    box off;
    ylabel(['DA ' num2str(i+1)]);
    ylim([0 500]);
    
    if i == 1
        set(gca, 'ticklength', [0.005 0.005], 'xcolor', [0 0 0]);
        set(gca,'XTICK',1:10*fps:1800);
        x_tick=get(gca,'XTICK');
        set(gca,'XTICKLABEL',(x_tick-1)/fps);
        xlabel('Time/s');
    else
        set(gca, 'ticklength', [0.005 0.005], 'xcolor', [1 1 1]);
    end
    
    if i == size(gfp, 2)
        ylim([0 1000]);
        for j = 1:size(period, 1);
            plot([period(j, 1) period(j, 2)], [max(ylim) max(ylim)], 'color', [0 0.6 1], 'linewidth', 3*thickness);
        end
        title('GCaMP6 and wCherry Signal', 'fontweight', 'bold');
        ylabel('DB 8/9');
    end
%     axis tight;
    xlim([1 1800]);
    
%     subplot(size(gfp, 2), 2, 2*(size(gfp, 2) - i +1));
%     plot(gfp_delta_over_median(:, i), 'color', [0 0.8 0], 'linewidth', thickness); hold on;
%     plot(ratio_delta_over_median(:, i), 'color', [0 0.4 1], 'linewidth', thickness); 
%     plot(ratio_delta_over_median_old(:, i+1), 'color', [0 0 0], 'linewidth', thickness);
%     plot(curv_norm(:, i+2), 'linestyle', ':', 'color', 'm', 'linewidth', thickness); 
%     box off;
%     if i == 1
%         set(gca, 'ticklength', [0.005 0.005], 'xcolor', [0 0 0]);
%         set(gca,'XTICK',1:10*fps:numframes);
%         x_tick=get(gca,'XTICK');
%         set(gca,'XTICKLABEL',(x_tick-1)/fps);
%         xlabel('Time/s');
%     else
%         set(gca, 'ticklength', [0.005 0.005], 'xcolor', [1 1 1]);
%     end
%     if i == size(gfp, 2)
%         plot([1 958], [max(ylim) max(ylim)], 'r', 'linewidth', 3*thickness);
% %         plot([408 423], [max(ylim) max(ylim)], 'r', 'linewidth', 3*thickness);
% %         plot([473 642], [max(ylim) max(ylim)], 'r', 'linewidth', 3*thickness);
%         title('Curvature Together with Normalized Signal Using Different Methods', 'fontweight', 'bold');
%     end
%     axis tight;
%     xlim([1 numframes]);

end

set(gcf, 'color', [1 1 1]);