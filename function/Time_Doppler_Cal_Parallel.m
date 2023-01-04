function g = Time_Doppler_Cal_Parallel(Col_Agg)
    tic
    Col_Agg_Parallel = gpuArray(Col_Agg);
    g = gpuArray(zeros(size(Col_Agg_Parallel, 1), size(Col_Agg_Parallel, 2)));
    f_delta = gpuArray(zeros(1, 20));
    n = size(Col_Agg_Parallel, 3);
    parfor i = 1:n
        link = gpuArray(Col_Agg_Parallel(:,:,i));
        [~, I] = max(link, [], 'all', 'linear');
        [rol, col] = ind2sub(size(link), I);
	    f_delta(i) = -40+2.*(col-1);
        g(i,:) = mapminmax(link(rol, :), 0, 5);
    end
    g(g < 2.5) = 0;
    T = 0:0.5:9.5;
    fD = -40:2:40;
    [A, B] = meshgrid(fD, T);
    surf(A, B, g), view(0, 90), colorbar
    saveas(gcf, 'Time_Doppler', 'fig');
    figure;
    f = 2120000000;
    c = 300000000;
    d = 70;
    x = 0;
    x_sum = zeros(1,20);
    v = 0;
    v_rela = f_delta./f.*c;
    v_real = zeros(1, 20);
    for i = 1:20
        d = d-x
        b = atan(250/d);
        v_real(i) = v_rela(i)./cos(b)
        x = v_real(i).*0.5
        x_sum(i) = x;
        v = v_real(i)
    end
    t = 0:0.5:9.5;
    figure;
    plot(t, v_rela);
    figure;
    plot(t, v_real);
    figure;
    plot(t, x_sum);
    toc
end