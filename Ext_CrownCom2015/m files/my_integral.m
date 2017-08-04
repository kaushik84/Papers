function[value] = my_integral(p, e, ts, te, lower_limit, upper_limit, Num_sample_pts)

    sample_pts = linspace(lower_limit,upper_limit, Num_sample_pts);
    func_value = e .* exp( ((-erfcinv(2 * sample_pts))).^2 - (e./(1 + sqrt(2) *...
        ts * (-erfcinv(2 * sample_pts))) - p).^2 /(sqrt(2) * p * te)^2)...
        * ts ./ (p * te * ( 1 + sqrt(2) * ts * (-erfcinv(2 * sample_pts))).^2);
    value = sum(func_value .* (sample_pts(end) - sample_pts(end -1)));
end