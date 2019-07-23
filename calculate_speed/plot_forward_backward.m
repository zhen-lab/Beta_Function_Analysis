function [] = plot_forward_backward(forward, backward)

cmp = lbmap(20, 'redblue');
clr_for = cmp(8,:); clr_back = [0 0.5 1];
neurons_separated = size(forward, 2);

for j = 1:neurons_separated
    
    s(1) = figure; hold on; 
    for i = 1:size(forward,1)
        if ~isempty(forward{i,j})
            plot(forward{i,j}, 'color', clr_for);            
        end
    end
    for i = 1:size(backward,1)
        if ~isempty(backward{i,j})
            plot(backward{i,j}, 'color', clr_back);
        end
    end
    set(gca, 'visible', 'off')
    hold off

    s(2) = figure; hold on; 
    for i = 1:size(forward,1)
        if ~isempty(forward{i,j})
            plot(forward{i,j}/forward{i,j}(1), 'color', clr_for);
        end
    end
    for i = 1:size(backward,1)
        if ~isempty(backward{i,j})
            plot(backward{i,j}/backward{i,j}(1), 'color', clr_back);
        end
    end
    xl = xlim;
    plot(xl, [1 1], 'k');
    set(gcf, 'color', 'w')
    set(gca, 'visible', 'off')
    hold off
    
    if neurons_separated>1
        switch j 
            case 1
                nn = 'Anterior';
            case 2
                nn = 'Posterior';
        end
    else
        nn = 'All';
    end
    
    savefig(s, [nn '_Neuron_Forward_Backward.fig']);
    
end

end