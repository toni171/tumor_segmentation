function curv = computeCurvature(phi_x,phi_y)
    norm_grad = sqrt(phi_x.^2 + phi_y.^2 + eps);
    n_x = phi_x ./ norm_grad;
    n_y = phi_y ./ norm_grad;
    curv = divCentral(n_x,n_y);
end