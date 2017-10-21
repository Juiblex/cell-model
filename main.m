clear all;

%% Parameters

k = struct();
% nutrient transport
k.s = 10000;
k.v_t = 726;
k.K_t  = 1000;
% metabolism
k.v_m = 5800;
k.K_m = 1000;
k.n_s = 0.5;
% protein lengths
k.n_r = 7459; % ribosome
k.n_x = 300;
% ribosome binding/unbinding
k.k_b = 1;
k.k_u = 1;
% transcription
k.w_r = 930; % transcription rates
k.w_e = 4.14;
k.w_q = 948.93;
k.w_p = k.w_e;
k.theta_r = 426.87; % transcription thresholds
k.theta_x = 4.38;
k.K_q = 152219; % q-protein autoinhibition
k.h_q = 4;
k.d_m = 0.1; % mRNA degradation rate
% translation
k.gamma_max = 1260;
k.K_gamma = 7;
k.M = 1e8;
% chloramphenicol
k.c_m = 4;
k.k_cm = 0.00599; %% binding rate
% LacI/TetR
k.w_TL = 4.14;
k.d_p = 0.05; % degradation rate for LacI and TetR
k.n_L = 600; %LacI is heavier
k.n_T = 300;
k.K_L = 600;
k.K_T = 500;
k.h_L = 4;
k.h_T = 4;

Tend = 10000;

%% Figure 1B - varying nutrient efficiency and chloramphenicol

% c_m = [0 2 4 8 12];
% n_s = [0.25 0.5 1 2 5];
% Y0_c = [0 10 10 10 10 10 0 0 0 0 0 0 0 0 0 0 0 0];
% 
% figure();
% 
% for i=1:length(n_s)
%     rmf = zeros(length(n_s), 1);
%     lambda = zeros(length(n_s), 1);
%     for j=1:length(c_m)
%        k.n_s = n_s(i);
%        k.c_m = c_m(j);
%        [T, Y] = ode23s(@diffequs_chlor, [0 Tend], Y0_c, [], k);
%        rib_mass = Y(end, 3) + sum(Y(end,11:18)); % free + bound ribosomes
%        rmf(j) = rib_mass/(rib_mass+Y(end, 4)+Y(end, 5)+Y(end, 6));
%        gamma = k.gamma_max*Y(end, 2)/(k.K_gamma+Y(end, 2));
%        lambda(j) = gamma*(Y(end, 11)+Y(end, 12)+Y(end, 13)+Y(end, 14))/k.M;
%     end
%     plot(rmf, lambda, '--o', 'DisplayName', strcat('n_s=', num2str(n_s(i)))); hold on;
% end
% xlabel('Growth rate (min^{-1})');
% ylabel('Ribosome mass fraction');
% title('Nutrient efficiency and chloramphenicol');
% legend('show');

%% Figure 2A - relative transcription rates

% k.c_m = 0;
% Y0_b = [0 10 10 10 10 10 0 0 0 0 0 0 0 0];
% opts = odeset('NonNegative', 1:length(Y0_b));
% 
% % we change intracellular activity levels via n_s
% %n \in 0.01..50 gives us a \in 0.05..5e7
% n_s = logspace(log10(0.01), log10(50), 50);
% 
% for i = 1:length(n_s)
%     k.n_s = n_s(i);
%     [T, Y] = ode15s(@diffequs_base, [0 Tend], Y0_b, opts, k);
%     omega_r(i) = k.w_r*Y(end,2)/(k.theta_r+Y(end,2));
%     omega_e(i) = k.w_e*Y(end,2)/(k.theta_x+Y(end,2));
%     omega_q(i) = k.w_q*Y(end,2)/(k.theta_x+Y(end,2))*1/(1+(Y(end,6)/k.K_q)^k.h_q);
%     a(i) = Y(end, 2);
% end
% 
% total = omega_r+omega_e+omega_q;
% figure();
% semilogx(n_s, omega_r./total, 'r', 'DisplayName', 'ribosome', 'LineWidth', 2); hold on;
% semilogx(n_s, omega_e./total, 'g', 'DisplayName', 'enzymes', 'LineWidth', 2); hold on;
% semilogx(n_s, omega_q./total, 'b', 'DisplayName', 'housekeeping', 'LineWidth', 2); hold on;
% xlabel('Nutrient efficiency');
% ylabel('Relative transcription rates');
% axis([0.01, 50, 0, 1]);
% legend('show');

%% Figure 2BCD - deletion strains
% 
% s_ext = logspace(0, 6, 13);
% Y0_p = [0 10 10 10 10 10 0 0 0 0 0 0 0 0 0 0 0];
% opts = odeset('NonNegative', 1:length(Y0_p));
% 
% %compute x_end for the wild-type
% for i=1:length(s_ext)
%     k.s = s_ext(i);
%     [T, Y_wt] = ode15s(@diffequs_grat, [0 Tend], Y0_p, opts, k);
%     enz_wt(i) = Y_wt(end, 4) + Y_wt(end, 5);
%     rib_wt(i) = Y_wt(end, 3);
%     hk_wt(i) = Y_wt(end, 6);
%     grat_wt(i) = Y_wt(end, 15);
%     
%     omega_r = k.w_r*Y_wt(end, 2)/(k.theta_r+Y_wt(end, 2));
%     omega_e = k.w_e*Y_wt(end, 2)/(k.theta_x+Y_wt(end, 2));
%     omega_q = k.w_q*Y_wt(end, 2)/(k.theta_x+Y_wt(end, 2))*1/(1+(Y_wt(end, 6)/k.K_q)^k.h_q);
%     omega_p = k.w_p*Y_wt(end, 2)/(k.theta_x+Y_wt(end, 2));
%     P_e_wt(i) = omega_r/(omega_r+omega_e+omega_q+omega_p);
%     P_p_wt(i) = omega_p/(omega_r+omega_e+omega_q+omega_p);
% end
% 
% % compute x_end for the enzyme deletion strain
% k.w_e = 4.14/2;
% for i=1:length(s_ext)
%     k.s = s_ext(i);
%     [T, Y_ed] = ode15s(@diffequs_grat, [0 Tend], Y0_p, opts, k);
%     enz_ed(i) = Y_ed(end, 4) + Y_ed(end, 5);
%     rib_ed(i) = Y_ed(end, 3);
%     hk_ed(i) = Y_ed(end, 6);
%     grat_ed(i) = Y_ed(end, 15);
%     
%     omega_r = k.w_r*Y_ed(end, 2)/(k.theta_r+Y_ed(end, 2));
%     omega_e = k.w_e*Y_ed(end, 2)/(k.theta_x+Y_ed(end, 2));
%     omega_q = k.w_q*Y_ed(end, 2)/(k.theta_x+Y_ed(end, 2))*1/(1+(Y_ed(end, 6)/k.K_q)^k.h_q);
%     omega_p = k.w_p*Y_ed(end, 2)/(k.theta_x+Y_ed(end, 2));
%     P_e_ed(i) = omega_r/(omega_r+omega_e+omega_q+omega_p);
% end
% 
% % compute x_end for the gratuitous protein deletion strain
% k.w_e = 4.14;
% k.w_p = k.w_e/2;
% for i=1:length(s_ext)
%     k.s = s_ext(i);
%     [T, Y_pd] = ode15s(@diffequs_grat, [0 Tend], Y0_p, opts, k);
%     enz_pd(i) = Y_pd(end, 4) + Y_pd(end, 5);
%     rib_pd(i) = Y_pd(end, 3);
%     hk_pd(i) = Y_pd(end, 6);
%     grat_pd(i) = Y_pd(end, 15);
%     omega_r = k.w_r*Y_pd(end, 2)/(k.theta_r+Y_pd(end, 2));
%     omega_e = k.w_e*Y_pd(end, 2)/(k.theta_x+Y_pd(end, 2));
%     omega_q = k.w_q*Y_pd(end, 2)/(k.theta_x+Y_pd(end, 2))*1/(1+(Y_pd(end, 6)/k.K_q)^k.h_q);
%     omega_p = k.w_p*Y_pd(end, 2)/(k.theta_x+Y_pd(end, 2));
%     P_p_pd(i) = omega_p/(omega_r+omega_e+omega_q+omega_p);
% end
% 
% % responsiveness
% for i=1:length(s_ext)
%     resp_enz_e(i) = log2(enz_ed(i)/(enz_wt(i)/2));
%     resp_rib_e(i) = log2(rib_ed(i)/rib_wt(i));
%     resp_hk_e(i) = log2(hk_ed(i)/hk_wt(i));
%     
%     resp_enz_p(i) = log2(enz_pd(i)/enz_wt(i));
%     resp_rib_p(i) = log2(rib_pd(i)/rib_wt(i));
%     resp_hk_p(i) = log2(hk_pd(i)/hk_wt(i));
%     resp_grat_p(i) = log2(grat_pd(i)/(grat_wt(i)/2));
% end
% 
% figure();
% bar([resp_enz_e; resp_rib_e; resp_hk_e]');
% set(gca, 'xticklabel', s_ext);
% legend('enzyme', 'ribosome', 'housekeeping');
% xlabel('External nutrient level');
% ylabel('Responsiveness');
% %title('Enzyme deletion responsiveness');
% 
% 
% figure();
% bar([resp_enz_p; resp_rib_p; resp_hk_p; resp_grat_p]');
% set(gca, 'xticklabel', s_ext);
% legend('enzyme', 'ribosome', 'housekeeping', 'gratuitous');
% xlabel('External nutrient level');
% ylabel('Responsiveness');
% %title('Gratuitous protein deletion responsiveness');
% 
% figure();
% bar([P_e_ed./P_e_wt; P_p_pd./P_p_wt]');
% set(gca, 'xticklabel', s_ext);
% legend('enzyme', 'gratuitous');
% xlabel('External nutrient level');
% ylabel('Relative transcription rate');
% %title('Relative transcription rates');

%% Toggle switch - varying induction rate
% initial conditions: final state of the basic version

% Y0_base = [0 10 10 10 10 0 0 0 0 0 0 0 0 0];
% [T, Y_base] = ode15s(@diffequs_base, [0 Tend], Y0_base, [], k);
% Y0_tog = [Y_base(end,:) 1000 0 0 0 0 0];
% 
% w_TL = logspace(log10(1), log10(1e4), 40); 
% 
% for i=1:length(w_TL)
%     k.w_TL = w_TL(i);
%     [T, Y_tog] = ode15s(@diffequs_toggle, [0 Tend], Y0_tog, [], k);
%     gamma = k.gamma_max*Y_tog(end, 2)/(k.K_gamma+Y_tog(end, 2));
%     gammatog(i) = gamma;
%     lambdatog(i) = gamma*(Y_tog(end, 11)+Y_tog(end, 12)+Y_tog(end, 13)...
%         +Y_tog(end, 14)+Y_tog(end, 19)+Y_tog(end, 20))/k.M;
%     nnu_r(i) = Y_tog(end, 11)*gamma; %multiply by n_x to be homogeneous
%     nnu_et(i) = Y_tog(end, 12)*gamma;
%     nnu_em(i) = Y_tog(end, 13)*gamma;
%     nnu_q(i) = Y_tog(end, 14)*gamma;
%     nnu_L(i) = Y_tog(end, 19)*gamma;
%     nnu_T(i) = Y_tog(end, 20)*gamma;
%     atog(i) = Y_tog(end, 2);
%     laci(i) = Y_tog(end, 15);
%     tetr(i) = Y_tog(end, 16);
%     mlaci(i) = Y_tog(end, 17);
%     mtetr(i) = Y_tog(end, 18);
%     claci(i) = Y_tog(end, 19);
%     ctetr(i) = Y_tog(end, 20);
%     olaci(i) = k.w_TL*Y_tog(end, 2)/(k.theta_x+Y_tog(end, 2))*1/(1+(Y_tog(end, 16)/k.K_L)^k.h_L);
%     otetr(i) = k.w_TL*Y_tog(end, 2)/(k.theta_x+Y_tog(end, 2))*1/(1+(Y_tog(end, 15)/k.K_T)^k.h_T);
% end
% figure();
% bar(log10(w_TL), [nnu_L+nnu_T; nnu_r; nnu_et+nnu_em; nnu_q]', 'stacked'); hold on;
% plot(log10(w_TL), lambdatog*k.M, 'LineWidth', 2);
% set(gca, 'xtick', log10(w_TL));
% set(gca,'xticklabel', w_TL);
% ylabel('Translation rate (aa/min)');
% xlabel('Toggle switch induction level (log)');
% axis([-0.05, 4.05, 0, 2.5e6]);
% legend('toggleswitch', 'ribosome', 'enzymes', 'housekeeping', 'growth rate');

%% Influence of induction of the gratuitous protein

% w_p = logspace(log10(1), log10(1e4), 40);
% Y0_p = [0 10 10 10 10 10 0 0 0 0 0 0 0 0 0 0 0];
% 
% for i=1:length(w_p)
%     k.w_p = w_p(i);
%     [T, Y_p] = ode15s(@diffequs_grat, [0 Tend], Y0_p, [], k);
%     gamma = k.gamma_max*Y_p(end, 2)/(k.K_gamma+Y_p(end, 2));
%     gammap(i)=gamma;
%     lambdap(i) = gamma*(Y_p(end, 11)+Y_p(end, 12)+Y_p(end, 13)...
%         +Y_p(end, 14)+Y_p(end, 17))/k.M;
%     nnu_r(i) = Y_p(end, 11)*gamma; %multiply by n_x to be homogeneous
%     nnu_et(i) = Y_p(end, 12)*gamma;
%     nnu_em(i) = Y_p(end, 13)*gamma;
%     nnu_q(i) = Y_p(end, 14)*gamma;
%     nnu_p(i) = Y_p(end, 17)*gamma;
%     agrat(i) = Y_p(end, 2);
%     grat(i) = Y_p(end, 15);
%     mgrat(i) = Y_p(end, 16);
%     cgrat(i) = Y_p(end, 17);
%     ograt(i) = k.w_p*Y_p(end, 2)/(k.theta_x+Y_p(end, 2));
% end
% figure();
% bar(log10(w_p), [nnu_p; nnu_r; nnu_et+nnu_em; nnu_q]', 'stacked'); hold on;
% plot(log10(w_p), lambdap*k.M, 'LineWidth', 2);
% legend('gratuitous', 'ribosome', 'enzymes', 'housekeeping', 'growth rate');

%% Comparison of protein, mRNA and complex levels in LacI/TetR vs gratuitous

% figure();
% subplot(2,3,1);
% plot(log10(w_p), laci, 'b', 'LineWidth', 2); hold on;
% plot(log10(w_p), tetr, 'r', 'LineWidth', 3); hold on;
% plot(log10(w_p), grat, 'g', 'LineWidth', 2); hold on;
% legend('LacI', 'TetR', 'Grat');
% subplot(2,3,2);
% plot(log10(w_p), mlaci, '--b', 'LineWidth', 2); hold on;
% plot(log10(w_p), mtetr, '--r', 'LineWidth', 3); hold on;
% plot(log10(w_p), mgrat, '--g', 'LineWidth', 2); hold on;
% legend('mLacI', 'mTetR', 'mGrat');
% subplot(2,3,3);
% plot(log10(w_p), claci, ':b', 'LineWidth', 2); hold on;
% plot(log10(w_p), ctetr, ':r', 'LineWidth', 3); hold on;
% plot(log10(w_p), cgrat, ':g', 'LineWidth', 2); hold on;
% legend('cLacI', 'cTetR', 'cGrat');
% subplot(2,3,4);
% plot(log10(w_p), olaci, '-.b', 'LineWidth', 2); hold on;
% plot(log10(w_p), otetr, '-.r', 'LineWidth', 3); hold on;
% plot(log10(w_p), ograt, '-.g', 'LineWidth', 2); hold on;
% legend('oLacI', 'oTetR', 'oGrat');
% subplot(2,3,5);
% plot(log10(w_p), gammatog, '-.b', 'LineWidth', 3); hold on;
% plot(log10(w_p), gammap, '-.g', 'LineWidth', 2); hold on;
% legend('gammatog', 'gammagrat');
% % a is similar to gamma
% % nu is proportional to c
% subplot(2,3,6);
% plot(log10(w_p), lambdatog, 'b', 'LineWidth', 3); hold on;
% plot(log10(w_p), lambdap, 'g', 'LineWidth', 2); hold on;
% legend('lambdaatog', 'lambdagrat');