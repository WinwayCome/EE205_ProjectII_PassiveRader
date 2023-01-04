function cor = rangeDoppler_fft(fileName, i)
    load(fileName);
    f = 12500000;
    t_s = 1/f_s;
    
    t1 = 0 : t_s : 0.5 - t_s;
    
    e_shift = exp(1j .* 2 .* pi .* 3*10^6 .* t1);
    [b1, a1] = butter(15, 9 * 10 ^ 6 / (f_s / 2));
    seq_sur_shift = seq_sur .* e_shift;
    seq_ref_shift = seq_ref .* e_shift;
    
    
    seq_sur_fil = filter(b1, a1, seq_sur_shift);
    seq_ref_fil = filter(b1, a1, seq_ref_shift);
    seq_ref_conj = conj(seq_ref_fil);
    tau_max = 0;
    fD_max = 0;

    cor = zeros(7, 51);

    N = f;
    tau_ind = 0;
    fD_ind = 0;
    cor_max = 0;
    freq_range = N / 2 - 25 : 1 : N / 2 + 25;
    for m = 0 : 6
        tmp = [zeros(1, m) seq_ref_conj(1 : N - m)];
        r_n = seq_sur_fil .* tmp;
        cor_tmp = abs((fftshift(fft(r_n, N))));
        cor(m + 1, :) = cor_tmp(freq_range);
    end

    for m = 1 : 7
        for n = 1 : 51
           if cor(m, n) > cor_max
                cor_max = cor(m, n);
                tau_ind = m;
                fD_ind = n;
           end
        end
    end
    tau_ind = tau_ind - 1;
    
    Range = 0 : 12 : 72;
    fD = -40 : 1.6 : 40;
    [A B] = meshgrid(fD, Range);
    for m = 1 :7
        for n = 1 :51
            if cor(m, n) < 0.005
                cor(m, n) = 0;
            end
        end
    end
    figure0 = figure;
    surf(A, B, cor)
    view(0, 90)
    colorbar
    xlabel("Doppler Frequency (Hz)"), 
    ylabel("Range (m)"),
    title(sprintf("Range-Doppler Spectrum [%f s 12m]", (i - 1) * 0.5));
    saveas(figure0, sprintf('/Users/weiter/Documents/MATLAB/Project_2/figure/data_%d_range_doppler_fft.jpg', i));
end
