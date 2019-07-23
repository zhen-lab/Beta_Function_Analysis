function [seg_cls, curv_cor] = closest_seg_curv(pos_midline, pos_neuron, curv_all)
    
seg_cls = zeros(size(pos_neuron, 2)/2, size(pos_midline, 2)/2);
frm_all = size(pos_midline, 2)/2;

for frm_num = 1:frm_all

    fn_x = 2 * frm_num - 1;
    fn_y = fn_x + 1;
    
    neu_frm_num = reshape(pos_neuron(frm_num, :), 2, [])'; % All neuron positions at given frm_num

    seg_cls(:, frm_num) = dsearchn(pos_midline(:, fn_x:fn_y), neu_frm_num);

end

curv_cor = zeros(size(seg_cls, 1), size(seg_cls, 2));

for curv_num = 1:size(seg_cls, 2)
    
    curv_sub = curv_all(:, curv_num);
    curv_cor(:, curv_num) = curv_sub(seg_cls(:, curv_num)); % Submatrix of curvatures corresponding to closest segment
    
end

curv_cor = curv_cor';