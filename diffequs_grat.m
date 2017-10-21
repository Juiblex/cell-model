function dy = diffequs_grat(t, y, k)
%% Differential equations of the model with a gratuitous protein

dy = zeros(17, 1); %s_i, a, r, e_t, e_m, q, mRNAs, complexes, gratuitous protein

gamma = k.gamma_max*y(2)/(k.K_gamma+y(2));
lambda = gamma*(y(11)+y(12)+y(13)+y(14)+y(17))/k.M;
nu_imp = y(4)*k.v_t*k.s/(k.K_t+k.s);
nu_cat = y(5)*k.v_m*y(1)/(k.K_m+y(1));
nu_r = y(11)*gamma/k.n_r;
nu_et = y(12)*gamma/k.n_x;
nu_em = y(13)*gamma/k.n_x;
nu_q = y(14)*gamma/k.n_x;
nu_p = y(17)*gamma/k.n_x;
omega_r = k.w_r*y(2)/(k.theta_r+y(2));
omega_e = k.w_e*y(2)/(k.theta_x+y(2));
omega_q = k.w_q*y(2)/(k.theta_x+y(2))*1/(1+(y(6)/k.K_q)^k.h_q);
omega_p = k.w_p*y(2)/(k.theta_x+y(2));

%nutriment
dy(1) = nu_imp - nu_cat - lambda*y(1);
%energy
dy(2) = k.n_s*nu_cat - lambda*y(2) - (k.n_r*nu_r+k.n_x*(nu_et+nu_em+nu_q+nu_p));

%ribosome
dy(3) = nu_r - lambda*y(3)...
    + (nu_r - k.k_b*y(3)*y(7) + k.k_u*y(11))...
    + (nu_et - k.k_b*y(3)*y(8) + k.k_u*y(12))...
    + (nu_em - k.k_b*y(3)*y(9) + k.k_u*y(13))...
    + (nu_q - k.k_b*y(3)*y(10) + k.k_u*y(14))...
    + (nu_p - k.k_b*y(3)*y(16) + k.k_u*y(17));
% other proteins
dy(4) = nu_et - lambda*y(4);
dy(5) = nu_em - lambda*y(5);
dy(6) = nu_q - lambda*y(6);

% mRNAs
dy(7) = omega_r - (lambda+k.d_m)*y(7) + nu_r - k.k_b*y(3)*y(7) + k.k_u*y(11);
dy(8) = omega_e - (lambda+k.d_m)*y(8) + nu_et - k.k_b*y(3)*y(8) + k.k_u*y(12);
dy(9) = omega_e - (lambda+k.d_m)*y(9) + nu_em - k.k_b*y(3)*y(9) + k.k_u*y(13);
dy(10) = omega_q - (lambda+k.d_m)*y(10) + nu_q - k.k_b*y(3)*y(10) + k.k_u*y(14);

% complexes
dy(11) = -lambda*y(11) + k.k_b*y(3)*y(7) - k.k_u*y(11) - nu_r;
dy(12) = -lambda*y(12) + k.k_b*y(3)*y(8) - k.k_u*y(12) - nu_et;
dy(13) = -lambda*y(13) + k.k_b*y(3)*y(9) - k.k_u*y(13) - nu_em;
dy(14) = -lambda*y(14) + k.k_b*y(3)*y(10) - k.k_u*y(14) - nu_q;

% gratuitous protein
dy(15) = nu_p - lambda*y(15);
dy(16) = omega_p - (lambda+k.d_m)*y(16) + nu_p - k.k_b*y(3)*y(16) + k.k_u*y(17);
dy(17) = -lambda*y(17) + k.k_b*y(3)*y(16) - k.k_u*y(17) - nu_p;
