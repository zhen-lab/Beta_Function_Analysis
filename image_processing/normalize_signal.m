function [signal_norm] = normalize_signal(signal)

    signal_max = max(signal);
    signal_max = repmat(signal_max, size(signal, 1), 1);

    signal_min = min(signal);
    signal_min = repmat(signal_min, size(signal, 1), 1);

    signal_del = signal - signal_min;
    signal_rng = signal_max - signal_min;
    signal_norm = signal_del ./ signal_rng;

end