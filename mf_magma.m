function CF=mf_magma(T)
    clamp = @(x,b,a) max(min(x,a),b); % restricts input X between a and b
    t = clamp(T,700,900);% , 700, 900);
    t1 = (t - 773.8) / 72.49;
    t2 = t1 .* t1;
    t3 = t2 .* t1;
    mf = 0.1111 * t3 - 0.2756 * t2 + 0.3087 * t1 + 0.7167;
    CF=clamp(mf,0.0, 1.);
end

