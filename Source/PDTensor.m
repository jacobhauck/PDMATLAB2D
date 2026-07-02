function Lambda = PDTensor(E_1, E_2, nu_12, eta_12_11, eta_12_22)
    % Enforce Cauchy's relation
    nu_21 = E_2 / E_1 * nu_12;
    G_12 = E_2 * nu_12 / (1 - nu_12*nu_21 - eta_12_11*eta_12_22);

    % Compute classical Cauchy elasticity tensor first
    s = [
        1/E_1,         -nu_12/E_1,       eta_12_11/E_1;
        -nu_12/E_1,       1/E_2,         eta_12_22/E_2;
        eta_12_11/E_1, eta_12_22/E_2, 1/G_12
    ];
    det_s = det(s);

    cauchy = zeros(2, 2, 2, 2);  % Only need the independent components
    cauchy(1, 1, 1, 1) = (E_2 - G_12 * eta_12_22^2) / (E_2^2*G_12*det_s);
    cauchy(1, 1, 2, 2) = (E_2*nu_12 + G_12*eta_12_11*eta_12_22) / (E_1*E_2*G_12*det_s);
    cauchy(1, 1, 1, 2) = -(eta_12_22*nu_12 + eta_12_11) / (E_1*E_2*det_s);
    cauchy(2, 2, 2, 2) = (E_1 - G_12*eta_12_11^2) / (E_1^2*G_12*det_s);
    cauchy(2, 2, 1, 2) = -(E_1*eta_12_22 + E_2*eta_12_11*nu_12) / (E_1^2*E_2*det_s);
    cauchy(1, 2, 1, 2) = (E_1 - E_2*nu_12^2) / (E_1^2*E_2*det_s);

    % Set the five independent components
    Lambda = zeros(2, 2, 2, 2);
    Lambda(1, 1, 1, 1) = 10*cauchy(1, 1, 1, 1) - 20*cauchy(1, 1, 2, 2) + 2*cauchy(2, 2, 2, 2);
    Lambda(1, 1, 1, 2) = 20*cauchy(1, 1, 1, 2) - 12*cauchy(2, 2, 1, 2);
    Lambda(1, 1, 2, 2) = (-10*cauchy(1, 1, 1, 1) + 76*cauchy(1, 1, 2, 2) - 10*cauchy(2, 2, 2, 2)) / 3;
    Lambda(2, 2, 1, 2) = -12*cauchy(1, 1, 1, 2) + 20*cauchy(2, 2, 1, 2);
    Lambda(2, 2, 2, 2) = 2*cauchy(1, 1, 1, 1) - 20*cauchy(1, 1, 2, 2) + 10*cauchy(2, 2, 2, 2);

    % Now use symmetry to set the other 11
    Lambda(2, 1, 1, 1) = Lambda(1, 1, 1, 2);
    Lambda(1, 2, 1, 1) = Lambda(1, 1, 1, 2);
    Lambda(1, 1, 2, 1) = Lambda(1, 1, 1, 2);
    Lambda(2, 2, 1, 1) = Lambda(1, 1, 2, 2);
    Lambda(2, 1, 2, 1) = Lambda(1, 1, 2, 2);
    Lambda(2, 1, 1, 2) = Lambda(1, 1, 2, 2);
    Lambda(1, 2, 2, 1) = Lambda(1, 1, 2, 2);
    Lambda(1, 2, 1, 2) = Lambda(1, 1, 2, 2);
    Lambda(2, 2, 2, 1) = Lambda(2, 2, 1, 2);
    Lambda(2, 1, 2, 2) = Lambda(2, 2, 1, 2);
    Lambda(1, 2, 2, 2) = Lambda(2, 2, 1, 2);
end