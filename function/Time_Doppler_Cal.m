function g = Time_Doppler_Cal(Col_Agg)
    tic
    g = zeros(size(Col_Agg, 1), size(Col_Agg, 2));
    f_delta = zeros(1, 20);
    n = size(Col_Agg, 3);
    for i = 1:n
        [~, I] = max(Col_Agg(:,:,i), [], 'all', 'linear');
        [rol, col] = ind2sub(size(Col_Agg(:,:,i)), I);
	    f_delta(i) = -40+2.*(col-1);
        g(i,:) = mapminmax(Col_Agg(rol, :, i), 0, 5);
        for k = 1:length(g(i,:))
            if g(i, k) < 2.5
                g(i, k) = 0;
            end
        end
    end
    T = 0:0.5:9.5;
    fD = -40:2:40;
    [A, B] = meshgrid(fD, T);
    surf(A, B, g), view(0, 90), colorbar
    saveas(gcf, 'Time_Doppler', 'fig');
    figure;
    f = 2120000000;
    c = 300000000;
    v_rela = f_delta./f.*c;
    t = 0:0.5:9.5;
    t = [t(1:2), t(4: 20)];
    v_rela = [v_rela(1:2), v_rela(4:20)];
    P = polyfit(t, v_rela, 4);
    t = 0:0.1:9.5;
    v_rela = P(1).* t.^4 + P(2).*t.^3 + P(3).*t.^2 + P(4).*t + P(5);
    v_real = zeros(1, length(v_rela));
    d = 70;
    x = 0;
    x_sum = zeros(1,20);
    for i = 1:length(v_rela)
        if i == 1
            d = d-x;
            b = atan(250/d);
            cos(b)
            v_real(i) = v_rela(i)./cos(b);
            x_sum(i) = x;
        else
            d = d-x;
            b = atan(250/d);
            cos(b)
            v_real(i) = v_rela(i)./cos(b);
            x = (v_real(i)+v_real(i-1))./2.*0.1;
            x_sum(i) = x_sum(i-1)+x;
        end
    end
    figure;
    plot(t, v_rela), title("v rela");
    figure;
    plot(t, v_real), title("v real");
    figure;
    plot(t, x_sum), title("x sum");
    toc
end