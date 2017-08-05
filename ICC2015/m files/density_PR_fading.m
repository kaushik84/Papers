function out = density_PR_fading(x_pts, K, tp, it, pl, np, nc, diff, operation)
    %% This function calculates numerically the probability and expectation

    out = 0;
    bins = length(x_pts);
    b_40 = zeros(1, bins);
    t_40 = zeros(1, bins);
    for i = 1:bins
        x_sym = sym(x_pts(i));
        t_40(i) = 4 * it * K^2 * nc * pl^2 * tp * x_sym/...
            (it * K * nc * pl + 2 * np * x_sym + K * pl * tp * x_sym)^2;
        b_40(i) = vpa(log(abs(hypergeom([(2 + K)/4, (4 + K) /4], K/2, t_40(i)))));
    end
    a_40 = (K + 2)/4 * log(it * nc./(tp * x_pts))  + (K + 2)/4 *...
        log( it * K^2 * nc * pl^2 * tp * x_pts...
        ./(it * K * nc * pl + 2 * np * x_pts + K * pl * tp * x_pts).^2);

    pdf_scaled_inv_ncx2_40 = diff * np./(it * nc * pl)...
        .* exp(a_40 + b_40); 

    if strcmp(operation,'dis')
        out  = sum(pdf_scaled_inv_ncx2_40);
    elseif strcmp(operation,'exp')
        out  = sum(x_pts .* pdf_scaled_inv_ncx2_40);
    else
        %% Error
    end
end