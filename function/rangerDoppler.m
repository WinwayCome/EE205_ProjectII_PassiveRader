function cor = rangeDoppler(fileName, i)
    load(fileName);
    f = 12500000;
    ts = 1/f_s;
    N = f_s;
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

    cor = zeros(7, 41);
    t2 = 0 + (i - 1) * 0.5 : ts : 0.5 + (i - 1) * 0.5;

    cor_max = 0.0;
    for m = 0 : 6
        for n = -40 : 2 : 40
            cor_tmp = 0;
            fD = n;
            e = exp(-1j .* 2.*pi.* fD .* t2);
            e = e(m + 1: f);
            cor_tmp = sum(e .* seq_sur_fil(m + 1 : f) .* seq_ref_conj(1 : f - (m + 1) + 1));
            cor_tmp = abs(cor_tmp);
            if cor_tmp > cor_max
                cor_max = cor_tmp;
                i_max = m;
                fD_max = fD;
            end
            cor(m + 1, fD / 2 + 21) = cor_tmp;
        end
    end
    tau_max = i_max .* 1 / f_s;
    Range = 0 : 12 : 72;
    fD = -40 : 2 : 40;
    [A B] = meshgrid(fD, Range);
    for m = 1 :7
        for n = 1 :41
            if cor(m, n) < 0.005
                cor(m, n) = 0;
            end
        end
    end
    figure0 = figure;
    surf(A, B, cor)
    view(0, 90)
    colorbar
    saveas(figure0, sprintf('/Users/Administrator/Documents/MATLAB/Project_2/figure/data_%d_range_doppler.jpg', i));
end
