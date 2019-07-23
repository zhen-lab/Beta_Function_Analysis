function [Y_ratio, Y_curv] = plot_fourier_norm(ratio_rev_lin_3_9_norm, curv_rev_lin_2_7_norm, fps)

Fs=fps;

L_ratio=length(ratio_rev_lin_3_9_norm);
NFFT = 2^nextpow2(L_ratio); % Next power of 2 from length of y
Y = fft(ratio_rev_lin_3_9_norm,NFFT)/L_ratio;
f = Fs/2*linspace(0,1,NFFT/2+1);

% Plot single-sided amplitude spectrum.
plot(f,2*abs(Y(1:NFFT/2+1)), 'Color', [0.9,0.1,0]);
title('Fourier transformation of \color{red}activity and curvature \color{black}of DA neurons in L1 in backward movement');
xlabel('Frequency/Hz');
ylabel('|Y(f)|');
hold on;
Y_ratio = 2*abs(Y(1:NFFT/2+1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L_curv=length(curv_rev_lin_2_7_norm);
NFFT = 2^nextpow2(L_curv); % Next power of 2 from length of y
Y = fft(curv_rev_lin_2_7_norm,NFFT)/L_curv;
f = Fs/2*linspace(0,1,NFFT/2+1);

% Plot single-sided amplitude spectrum.
plot(f,2*abs(Y(1:NFFT/2+1)), 'Color', [0,0.5,0.9]);
title('Fourier transformation of \color{red}activity and curvature \color{black}of DA neurons in L1 in backward movement');
xlabel('Frequency/Hz');
ylabel('|Y(f)|');
Y_curv = 2*abs(Y(1:NFFT/2+1));
end

