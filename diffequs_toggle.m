function dy = diffequs_base(t, y, k)
%% Differential equations of the model with the toggle switch

dy = zeros(20, 1); %s_i, a, r, e_t, e_m, q, mRNAs, complexes, LacI, TetR

gamma = k.gamma_max*y(2)/(k.K_gamma+y(2));
lambda = gamma*(y(11)+y(12)+y(13)+y(14)+y(19)+y(20))/k.M;
nu_imp = y(4)*k.v_t*k.s/(k.K_t+k.s);
nu_cat = y(5)*k.v_m*y(1)/(k.K_m+y(1));
nu_r = y(11)*gamma/k.n_r;
nu_et = y(12)*gamma/k.n_x;
nu_em = y(13)*gamma/k.n_x;
nu_q = y(14)*gamma/k.n_x;
nu_L = y(19)*gamma/k.n_L;
nu_T = y(20)*gamma/k.n_T;
omega_r = k.w_r*y(2)/(k.theta_r+y(2));
omega_e = k.w_e*y(2)/(k.theta_x+y(2));
omega_q = k.w_q*y(2)/(k.theta_x+y(2))*1/(1+(y(6)/k.K_q)^k.h_q);
omega_L = k.w_TL*y(2)/(k.theta_x+y(2))*1/(1+(y(16)/k.K_L)^k.h_L); % LacI and TetR inhibit
omega_T = k.w_TL*y(2)/(k.theta_x+y(2))*1/(1+(y(15)/k.K_T)^k.h_T); % each other's transcription

%nutriment
dy(1) = nu_imp - nu_cat - lambda*y(1);
%energy
dy(2) = k.n_s*nu_cat - lambda*y(2)...
    - (k.n_r*nu_r + k.n_x*(nu_et+nu_em+nu_q) + k.n_L*nu_L + k.n_T*nu_T);

%ribosome
dy(3) = nu_r - lambda*y(3)...
    + (nu_r - k.k_b*y(3)*y(7) + k.k_u*y(11))...
    + (nu_et - k.k_b*y(3)*y(8) + k.k_u*y(12))...
    + (nu_em - k.k_b*y(3)*y(9) + k.k_u*y(13))...
    + (nu_q - k.k_b*y(3)*y(10) + k.k_u*y(14))...
    + (nu_L - k.k_b*y(3)*y(17) + k.k_u*y(19))...
    + (nu_T - k.k_b*y(3)*y(18) + k.k_u*y(20));
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

% LacI and TetR
dy(15) = nu_L - (lambda)*y(15);
dy(16) = nu_T - (lambda)*y(16);

dy(17) = omega_L - (lambda+k.d_m)*y(17) + nu_L - k.k_b*y(3)*y(17) + k.k_u*y(19);
dy(18) = omega_T - (lambda+k.d_m)*y(18) + nu_T - k.k_b*y(3)*y(18) + k.k_u*y(20);

dy(19) = -lambda*y(19) + k.k_b*y(3)*y(17) - k.k_u*y(19) - nu_L;
dy(20) = -lambda*y(20) + k.k_b*y(3)*y(18) - k.k_u*y(20) - nu_T;