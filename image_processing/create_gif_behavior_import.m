function [avgint] = create_gif_behavior_import(close_para, fill_para, thresh, adpt_max, adpt_min)

setup_proof_reading;
fprintf('recording has been loaded \n');

avgint = create_gif_behavior(filename, imagelist, close_para, fill_para, thresh, adpt_max, adpt_min);

end