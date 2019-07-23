function [] = ratio_dorsal_ventral_norm(ratio_d, ratio_v)
   
    if size(ratio_d == ratio_v)

        fps = 20; 
        numframes = length(ratio_d);
        d_v = ratio_d./ratio_v;
        
        figure(7);
        
        subplot(2,1,2);
        
        range = 1;
        semilogy(d_v, 'Color', [1 1 1]*0.2, 'LineWidth', 2); ylim([10^(-range) 10^(range)]);
%         h = findobj(gcf, 'Type', 'axes');
%         set(h, 'YTick', [10^(-range) 10^(-range/2) 10^0 10^(range/2) 10^(range)], 'YTickLabel', [10^(-range) 10^(-range/2) 10^0 10^(range/2) 10^(range)]);
        hold on;
        area(d_v, 1, 'FaceColor', [1 1 1]*0.8, 'EdgeColor', 'none');
        plot([0 length(ratio_d)], [1 1], 'Color', 'k', 'LineStyle', '--');
        hold on;
        t = text(length(ratio_d)-20, max(d_v), ['D/V = ' num2str(mean(d_v), 2)], 'FontSize', 8);
        set(t, 'HorizontalAlignment', 'right');
        ylabel('Dorsal/Ventral'); xlabel('Time/s');
        
        set(gca,'XTICK',1:10*fps:numframes);
        x_tick=get(gca,'XTICK');
        set(gca,'XTICKLABEL',(x_tick-1)/fps);
        box off; set(gca, 'TickLength', [0 0]);
        xlim([1 size(ratio_d, 1)]);
        
        
        subplot(2,1,1); 
        plot(ratio_d, 'Color', [1 0.2 0], 'LineWidth', 2); 
        hold on; 
        plot(ratio_v, 'Color', [0 0.6 1], 'LineWidth', 2); 
        ylabel('Dorsal and Ventral signal [a.u.]');
        box off; set(gca, 'XColor', get(gca, 'Color'));
        xlim([1 size(ratio_d, 1)]); ylim([0 2]);

        set(gcf, 'Color', get(gca, 'Color'));
    
    end

end