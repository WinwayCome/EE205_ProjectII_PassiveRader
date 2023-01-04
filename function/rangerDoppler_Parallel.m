function cor = rangerDoppler_Parallel(fileName, i)
    load(fileName);
    f = 12500000;
    ts = 1/f_s;
    t_s = 1/f_s;
    
    t1 = 0 : t_s : 0.5 - t_s;
    
    e_shift = gpuArray(exp(1j .* 2 .* pi .* 3*10^6 .* t1));
    [b1, a1] = butter(15, 9 * 10 ^ 6 / (f_s / 2));
    seq_sur_p = gpuArray(seq_sur);
    seq_ref_p = gpuArray(seq_ref);
    seq_sur_shift = gpuArray(seq_sur_p .* e_shift);
    seq_ref_shift = gpuArray(seq_ref_p .* e_shift);
    
    
    seq_sur_fil = gpuArray(filter(b1, a1, seq_sur_shift));
    seq_ref_fil = gpuArray(filter(b1, a1, seq_ref_shift));
    seq_ref_conj = gpuArray(conj(seq_ref_fil));

    cor = gpuArray(zeros(7, 41));
    t2 = gpuArray(0 + (i - 1) * 0.5 : ts : 0.5 + (i - 1) * 0.5);

    for m = 0 : 6
        link1 = gpuArray(seq_sur_fil(m + 1 : f));
        link2 = gpuArray(seq_ref_conj(1 : f - (m + 1) + 1));
        link3 = zeros(1, 41);
        parfor n = 1 : 41
            fD = (n-1) * 2 - 40;
            e = gpuArray(exp(-1j .* 2.*pi.* fD .* t2));
            e = e(m + 1: f);
            cor_tmp = sum(e .* link1 .* link2);
            cor_tmp = abs(cor_tmp);
            link3(n) = cor_tmp;
        end
        cor(m+1, :) = link3;
    end

    Range = 0 : 12 : 72;
    fD = -40 : 2 : 40;
    [A, B] = meshgrid(fD, Range);
    cor(cor < 0.005) = 0;
    figure0 = figure;
    surf(A, B, cor)
    view(0, 90)
    colorbar
    saveas(figure0, sprintf('C:/Users/23098/Desktop/信号与系统/ProjectII/fig/data_%d_range_doppler.jpg', i));
end
