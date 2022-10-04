function [wavelength_all] = collect_wavelength_data(folder)

f = dir(fullfile(folder, '*.mat'));
numsamples = numel(f);
wavelength_all = zeros(numsamples, 1);

for idx = 1 : numsamples

    pathname = fullfile(folder, f(idx).name);
    load(pathname, 'wavelength_over_wormlengthcorrected');
    wavelength_all(idx) = wavelength_over_wormlengthcorrected;
   
end

parts = strsplit(pwd, '\');
partlast = parts{1,end};
fname = [partlast '_wavelength_all.mat'];

save(fname, 'wavelength_all');

end