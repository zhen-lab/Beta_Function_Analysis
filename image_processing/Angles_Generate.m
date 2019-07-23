f = dir(fullfile(pwd, '*.xlsx'));

for i = 1:length(f)
    dataset = xlsread(f(i).name);
    angles = dataset(:, 3:35);
    angles_vec = reshape(angles, [], 1);
    xlswrite([f(i).name(1:end-5) '_angles.xlsx'], angles_vec);
end