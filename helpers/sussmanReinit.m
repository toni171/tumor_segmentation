function phi = sussmanReinit(phi,n_iters,dt)
    sign_phi = phi ./ sqrt(phi.^2 + eps);
    idx_pos = sign_phi > 0;
    idx_neg = sign_phi < 0;
    for i = 1:n_iters
        a = phi - circshift(phi,[0 1]);
        b = circshift(phi,[0 -1]) - phi;
        c = phi - circshift(phi,[1 0]);
        d = circshift(phi,[-1 0]) - phi;
        a_p = max(a,0); a_n = min(a,0);
        b_p = max(b,0); b_n = min(b,0);
        c_p = max(c,0); c_n = min(c,0);
        d_p = max(d,0); d_n = min(d,0);
        grad_plus = sqrt( max(a_p.^2, b_n.^2) + max(c_p.^2, d_n.^2) );
        grad_minus = sqrt( max(a_n.^2, b_p.^2) + max(c_n.^2, d_p.^2) );
        D = zeros(size(phi));
        D(idx_pos) = grad_plus(idx_pos) - 1;
        D(idx_neg) = grad_minus(idx_neg) - 1;
        phi = phi - dt .* sign_phi .* D;
    end
end