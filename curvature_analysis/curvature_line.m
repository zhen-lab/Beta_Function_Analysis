k = gfp';

for i = 1:length(np)
line = np(i, :);
line = reshape(line, 2, [])';
k(:, i) = LineCurvature2D(line);
end

k = k';
k_norm = normalize_signal(k);