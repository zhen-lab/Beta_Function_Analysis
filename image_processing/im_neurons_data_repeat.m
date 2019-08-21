%% Set up inputs

input = ratio_norm;
ymin = 0;
ymax = 2.5;

s_exp = []; 
s_ctrl = [];

%% Plot green lines

close all; 
s_exp = im_neurons_data(input, [ymin ymax], [ymin ymax], [0 1]);

%% Plot black lines

close all; 
s_ctrl = im_neurons_data(input, [ymin ymax], [ymin ymax], [1 0]);

%% Save figures

parts = strsplit(pwd, '\');

if ~isempty(s_exp)
    savefig(s_exp, ['Experimental_' parts{end} '.fig']);
    for i = 1:size(s_exp, 2)
        saveas(s_exp(i), ['Experimental_' parts{end} '.tif'], 'tiffn');
    end
    fprintf('figures saved for EXPERIMENTAL group. \n');
else
    if ~isempty(s_ctrl)
        savefig(s_ctrl, ['Control_' parts{end} '.fig']);
        for i = 1:size(s_ctrl, 2)
            saveas(s_ctrl(i), ['Control_' parts{end} '.tif'], 'tiffn');
        end
        fprintf('figures saved for CONTROL group. \n');
    end
end