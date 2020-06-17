clear all; clc; close all; fc = 0;

%% load the data
import_secim2004;
import_secim2009;
import_secim2014;
import_gsyh; % 2004 to 2017
import_nufus_data; % 2008 to 2018
import_girisim; % 2009 to 2017
import_goc_data; % 2008 to 2018
import_issizlik_data; % 2004 to 2019
import_AKP; % 2007 to 2019
import_Kurd; % 2008
import_betweenness; % 2008 to 2019
import_correlations; % 2008 to 2019

% NAN cells are converted to 0.
Goc(isnan(Goc)) = 0;

load('ilmeseafe_data.mat');

N = 81;
%%
algorithms_to_run = [0 0 1];
M_vec = [100000 500000 500000];

%% U, V, Z, and Y are formed

Years_all = 2004:2018;
L_y = length(Years_all);
U_all = cell(1, L_y);
V_all = cell(1, L_y);
Y_all = cell(1, L_y);
Z_all = cell(1, L_y);

lu_ind = 0;
lv_ind = 0;
lz_ind = 0;

one_way_factors_u = {'GDP -->', 'Unemployment -->', 'AKP -->', ...
    'Popularity -->', 'Log-population -->', 'Betweenness -->', 'Kurd -->'};
one_way_factors_v = {'GDP <--', 'Unemployment <--', 'AKP <--', ...
    'Popularity <--', 'Log-population <--', 'Betweenness <--', 'Kurd <--'};
two_way_factors = {'Political distance <-->', 'Spatial distance <-->', ...
    'Previous year <-->', 'Reciprocity <-->', 'In-flow similarity <-->'};

% select_vec_U = [1 1 1 0 1 1 1];
% select_vec_V = [1 1 1 1 1 1 1];
% select_vec_Z = [1 1 1 1 1];

% run this for the empty model (for minimal computation)
select_vec_U = [1 0 0 0 0 0 0];
select_vec_V = [1 0 0 0 0 0 0];
select_vec_Z = [1 0 0 0 0];
selected_factors_U = find(select_vec_U);
selected_factors_V = find(select_vec_V);

% construct legends
legends_U = one_way_factors_u(select_vec_U == 1);
legends_V = one_way_factors_v(select_vec_V == 1);
legends_Z = two_way_factors(select_vec_Z == 1);

for t = 1:L_y
    year = Years_all(t);
    U_all{t} = [];
    V_all{t} = [];
    
    %%%%%%%%%%%%%%%%%%% GSYH  %%%%%%%%%%%%%%%%%%%
    ind = find(GSYH_years == year-1);
    if select_vec_U(1) == 1
        if isempty(ind) == 0
            U_all{t} = [U_all{t} GSYH(:, ind)];
        end
    end
    if select_vec_V(1) == 1
        if isempty(ind) == 0
            V_all{t} = [V_all{t} GSYH(:, ind)];
        end
    end
    
%     %%%%%%%%%%%%%%%%%%% Issizlik %%%%%%%%%%%%%%%%%%%
%     ind1 = find(Issizlik_years == year-1);
%     ind2 = find(Issizlik_years == year-2);
%     if select_vec_U(2) == 1
%         if (isempty(ind1) == 0)
%             if (isempty(ind2) == 0)
%                 U_all{t} = [U_all{t} Issizlik_iller(:, ind1) - Issizlik_iller(:, ind2)];
%             else
%                 U_all{t} = [U_all{t} zeros(N, 1)];
%             end
%         end
%     end
%     if select_vec_V(2) == 1
%         if (isempty(ind1) == 0)
%             if (isempty(ind2) == 0)
%                 V_all{t} = [V_all{t} Issizlik_iller(:, ind1) - Issizlik_iller(:, ind2)];
%             else
%                 V_all{t} = [V_all{t} zeros(N, 1)];
%             end
%         end
%     end
%     
    %%%%%%%%%%%%%%%%%%% Issizlik %%%%%%%%%%%%%%%%%%%
    ind = find(Issizlik_years == year-1);
    if select_vec_U(2) == 1
        if (isempty(ind) == 0)
            U_all{t} = [U_all{t} Issizlik_iller(:, ind)];
        end
    end
    if select_vec_V(2) == 1
        if (isempty(ind) == 0)
            V_all{t} = [V_all{t} Issizlik_iller(:, ind)];
        end
    end
    
    %%%%%%%%%%%%%%%%%%% AKP  %%%%%%%%%%%%%%%%%%%
    ind = find(AKP_years == year-1);
    if select_vec_U(3) == 1
        if (isempty(ind) == 0)
            U_all{t} = [U_all{t} AKP(:, ind)];
        end
    end
    if select_vec_V(3) == 1
        if (isempty(ind) == 0)
            V_all{t} = [V_all{t} AKP(:, ind)];
        end
    end
    
    %%%%%%%%%%%%%%%%%%% Popularity %%%%%%%%%%%%%%%%%%%
    ind = find(Nufus_years == year-1);
    if select_vec_U(4) == 1
        if (isempty(ind) == 0)
            U_all{t} = [U_all{t} sum(Goc(:, :, ind))'];
        end
    end
    if select_vec_V(4) == 1
        if (isempty(ind) == 0)
            V_all{t} = [V_all{t} sum(Goc(:, :, ind))'];
        end
    end
    
    %%%%%%%%%%%%%%%%%%% Nufus %%%%%%%%%%%%%%%%%%%
    ind = find(Nufus_years == year-1);
    if select_vec_U(5) == 1
        if (isempty(ind) == 0)
            U_all{t} = [U_all{t} log(Nufus(:, ind))];
        end
    end
    if select_vec_V(5) == 1
        if (isempty(ind) == 0)
            V_all{t} = [V_all{t} log(Nufus(:, ind))];
        end
    end
    
    %%%%%%%%%%%%%%%%%%% Betweenness %%%%%%%%%%%%%%%%%%%
    ind = find(Betweenness_years == year-1);
    if select_vec_U(6) == 1
        if (isempty(ind) == 0)
            U_all{t} = [U_all{t} Betweenness(:, ind)];
        end
    end
    if select_vec_V(6) == 1
        if (isempty(ind) == 0)
            V_all{t} = [V_all{t} Betweenness(:, ind)];
        end
    end
    % Construct Z{t}
    z_ind = 0;
    
    %%%%%%%%%%%%%%%%%%% Kurd  %%%%%%%%%%%%%%%%%%%
    if select_vec_U(7) == 1
        U_all{t} = [U_all{t} Kurd];
    end
    if select_vec_V(7) == 1
        V_all{t} = [V_all{t} Kurd];
    end
    
    %%%%%%%%%%%%%%%%%%% Political distance %%%%%%%%%%%%%%%%%%%
    if select_vec_Z(1) == 1
        z_ind = z_ind + 1;
        if year >= 2004 && year < 2009
            secim_mtx = Secim2004;
        elseif year >= 2009 && year < 2014
            secim_mtx = Secim2009;
        elseif year >= 2014
            secim_mtx = Secim2014;
        end
        if year >= 2004
            D = zeros(N);
            for i = 1:N
                for j = 1:(i-1)
                    D(i, j) = 0.5*sum(abs(secim_mtx(i, :) - secim_mtx(j, :)));
                end
            end
            D = (D + D')/100;
            Z_all{t}(:, :, z_ind) = D;
        end
    end
    
    %%%%%%%%%%%%%%%%%%% Geographical distance %%%%%%%%%%%%%%%%%%%
    if select_vec_Z(2) == 1
        z_ind = z_ind + 1;
        Z_all{t}(:, :, z_ind) = ilmesafe_abc_abc;
    end
    
    %%%%%%%%%%%%%%%%%%% Autocorrelation %%%%%%%%%%%%%%%%%%%
    if select_vec_Z(3) == 1
        ind = find(Nufus_years == year-1);
        if isempty(ind) == 0
            z_ind = z_ind + 1;
            Z_all{t}(:, :, z_ind) = Goc(:, :, ind);
        end
    end
    
    %%%%%%%%%%%%%%%%%%% Reciprocity %%%%%%%%%%%%%%%%%%%
    if select_vec_Z(4) == 1
        ind = find(Nufus_years == year-1);
        if isempty(ind) == 0
            z_ind = z_ind + 1;
            Z_all{t}(:, :, z_ind) = Goc(:, :, ind)'/max(max(Goc(:, :, ind)));
        end
    end
    
    %%%%%%%%%%%%%%%%%%% Correaltions %%%%%%%%%%%%%%%%%%%
    if select_vec_Z(4) == 1
        ind = find(Correlations_years == year-1);
        if isempty(ind) == 0
            z_ind = z_ind + 1;
            Z_all{t}(:, :, z_ind) = Correlations(:, :, ind);
        end
    end
    
    
    %%%%%%%%%%%%%%%%%%% Construct Y %%%%%%%%%%%%%%%%%%%
    ind1 = find(Nufus_years == year);
    ind2 = find(Nufus_years == year-1);
    if (isempty(ind1) == 0) && (isempty(ind2) == 0)
        Y_all{t} = Goc(:, :, ind1);
        % diagonals
        Y_all{t}((0:N-1)*N + (1:N)) = Nufus(:, ind2) - sum(Y_all{t}, 2);
    end
end

% Common parameters
init_year = 6;
last_year = 15;
T = last_year - init_year + 1;
N = 81;
K1 = sum(select_vec_U);
K2 = sum(select_vec_V);
L = sum(select_vec_Z);

%% Get the inputs
Y = Y_all(init_year:last_year);
U = U_all(init_year:last_year);
V = V_all(init_year:last_year);
Z = Z_all(init_year:last_year);

%% normalise U's, V's, and Z's
max_abs_U = zeros(1, K1);
mean_U = zeros(1, K1);
max_abs_V = zeros(1, K2);
mean_V = zeros(1, K2);
max_abs_Z = zeros(1, L);
mean_Z = zeros(1, L);

for k = 1:K1
    % first, find the mean:
    temp_mean_vec = zeros(1, T);
    for t = 1:T
        temp_mean_vec(t) = mean(U{t}(:, k));
    end
    temp_mean = mean(temp_mean_vec);
    % then, find the max of (absolute value of) centralised values
    temp_max = 0;
    temp_min = inf;
    for t = 1:T
        temp_max = max(temp_max, max(U{t}(:, k)-temp_mean));
        temp_min = min(temp_min, min(U{t}(:, k)-temp_mean));
    end
    max_abs = max(temp_max, abs(temp_min));
    for t = 1:T
        U{t}(:, k) = (U{t}(:, k) - temp_mean)/max_abs;
    end
    max_abs_U(k) = max_abs;
    mean_U(k) = temp_mean;
end

% normalise V's
for k = 1:K2
    % first, find the mean:
    temp_mean_vec = zeros(1, T);
    for t = 1:T
        temp_mean_vec(t) = mean(V{t}(:, k));
    end
    temp_mean = mean(temp_mean_vec);
    % then, find the max of (absolute value of) centralised values
    temp_max = 0;
    temp_min = inf;
    for t = 1:T
        temp_max = max(temp_max, max(V{t}(:, k)-temp_mean));
        temp_min = min(temp_min, min(V{t}(:, k)-temp_mean));
    end
    max_abs = max(temp_max, abs(temp_min));
    for t = 1:T
        V{t}(:, k) = (V{t}(:, k) - temp_mean)/max_abs;
    end
    max_abs_V(k) = max_abs;
    mean_V(k) = temp_mean;
end

% normalise Z's
for l = 1:L
    % first, find the mean:
    temp_mean_vec = zeros(1, T);
    for t = 1:T
        temp_mean_vec(t) = mean(mean(Z{t}(:, :, l)));
    end
    temp_mean = mean(temp_mean_vec);

    temp_max = 0;
    temp_min = inf;
    for t = 1:T
        temp_max = max(temp_max, max(max(Z{t}(:, :, l)-temp_mean)));
        temp_min = min(temp_min, min(min(Z{t}(:, :, l)-temp_mean)));
    end
    max_abs = max(temp_max, abs(temp_min));
    for t = 1:T
        Z{t}(:, :, l) = (Z{t}(:, :, l) - temp_mean)/max_abs;
    end
    max_abs_Z(l) = max_abs;
    mean_Z(l) = temp_mean;
end


%% MCMC for the first model (Multinomial)
if algorithms_to_run(1) == 1
    % initial values
    Theta_init = {ones(1, 1), ones(K1, 1), ones(K2, 1), ones(L, 1)};
    
    % algorithm parameters
    sigma_prop = 0.01; % proposal std
    M = 10000; % number of iterations
    [Theta_samp] = MCMC_goc_model_1(Y, U, V, Z, M, Theta_init, sigma_prop);
    
    fc = fc + 1; figure(fc);
    plot(Theta_samp');
    legend('base', 'GSYH goc verim', 'AKP goc verim',...
        'issizlik goc verim', 'GSYH goc alim', 'AKP goc alim', 'popularity',...
        'issizlik goc alim', 'political distance', 'geographical distance', 'last year', 'reciprocity');
end

%% MCMC for the second model (Dirichlet-Multinomial)
if algorithms_to_run(2) == 1
    Theta_init = {zeros(K1, 1), zeros(K2, 1), zeros(L, 1), ones(1, 1), -5*ones(1, 1), zeros(1, 1)};
    M = M_vec(2);
    prop_params.sigma_prop = 0.0;
    prop_params.sigma_sc_prop = 0.01;
    prop_params.sigma0_a_prop = 0.01;
    prop_params.sigma0_b_prop = 0.01;

    % [Theta_samp_2a] = MCMC_goc_model_2a(Y, U, V, Z, M, Theta_init, prop_params);
    
    legends = [legends_U, legends_V, legends_Z, 'scale parameter', 'intercept', 'time slope'];
    
    range_conv = (M/2+1):M;
    
    D = K1+K2+L+3;
    
    fc = fc + 1; figure(fc);
    K_max = max(K1, K2);
    for i = 1:(K1+K2+L+1+2)
        if i <= K1
            j = selected_factors_U(i);
            J = K_max;
        elseif i > K1 && i <= (K1 + K2)
            j = K_max + selected_factors_V(i - K1);
            J = K_max;
        elseif i > (K1+K2) && i <= (K1 + K2 + L)
            J = L;
            j = 2*L + (i - K1 - K2);
        else
            J = 3;
            j = 3*3 + (i - K1 - K2 - L);
        end
        subplot(4, J, j);
        histogram(Theta_samp_2a(i, range_conv), 20 ,'Normalization', 'probability');
        set(gca, 'ytick', []);
        title(legends{i});
        grid on;
    end
    print(gcf,'full_model_histograms_from_30_May', '-depsc');
    
    fc = fc + 1; figure(fc);
    range_conv2 = M/10+1:M;
    K_max = max(K1, K2);
    for i = 1:(K1+K2+L+1+2)
        if i <= K1
            j = selected_factors_U(i);
            J = K_max;
        elseif i > K1 && i <= (K1 + K2)
            j = K_max + selected_factors_V(i - K1);
            J = K_max;
        elseif i > (K1+K2) && i <= (K1 + K2 + L)
            J = L;
            j = 2*L + (i - K1 - K2);
        else
            J = 3;
            j = 3*3 + (i - K1 - K2 - L);
        end
        subplot(4, J, j);
        plot(range_conv2, Theta_samp_2a(i, range_conv2)');
        title(legends{i});
        set(gca, 'xlim', [range_conv2(1) range_conv2(end)])
        if j <= 4*(J-1)
            set(gca, 'xtick', []);
        end
        grid on;
        
    end
    print(gcf,'full_model_trace_plots_from_30_May', '-depsc');

    
    fc = fc + 1; figure(fc);
    boxplot(Theta_samp_2a(1:(K1+K2+L), range_conv)', 'boxstyle', 'outline',...
        'orientation', 'horizontal', 'labels', legends(1:(K1+K2+L)), ...
        'factordirection', 'list', 'labelorientation', 'inline',...
        'outliersize', 0.1);
    grid on;
    hold on;
    plot((-9:9), (K1+0.5)*ones(1, 19), '-.k');
    plot((-9:9), (K1+K2+0.5)*ones(1, 19), '-.k');
    plot([0 0], [0 K1 + K2 + L+0.5], 'k');
    hold off;
    print(gcf,'full_model_boxplots_from_30_May', '-depsc');
    
    
    % Get the mean, std, and quantiles
    m_2a = zeros(1, D);
    s_2a = zeros(1, D);
    p1_2a = zeros(1, D);
    p2_2a = zeros(1, D);

    
    for i = 1:D
        m_2a(i) = mean(Theta_samp_2a(i, range_conv));
        s_2a(i) = std(Theta_samp_2a(i, range_conv));
        p1_2a(i) = quantile(Theta_samp_2a(i, range_conv), 0.05);
        p2_2a(i) = quantile(Theta_samp_2a(i, range_conv), 0.95);
    end
    

end

%% MCMC for the second model, but with different bases and year effects for provinces
set(0,'DefaultAxesTitleFontWeight','normal');

if algorithms_to_run(3) == 1

    Theta_init = {zeros(K1, 1), zeros(K2, 1), zeros(L, 1), ones(1, 1), ...
        -5*ones(N, 1), zeros(N, 1), [-5; 1], [0; 1]};
    
%     Theta_init = {Theta_samp_2a_ext(1:K1, end), Theta_samp_2a_ext(K1+1:K1+K2, end),...
%         Theta_samp_2a_ext(K1+K2+1:K1+K2+L, end), Theta_samp_2a_ext(K1+K2+L+1, end), ...
%         Theta_samp_2a_ext(K1+K2+L+1+(1:N), end), Theta_samp_2a_ext(K1+K2+L+1+N+(1:N), end),...
%         Theta_samp_2a_ext(K1+K2+L+1+2*N+(1:2), end), Theta_samp_2a_ext(K1+K2+L+1+2*N+(3:4), end)};
    
%     load('Theta_init_for_2a');
    
    M = M_vec(3);
    prop_params.sigma_prop = 0;
    prop_params.sigma_sc_prop = 0.01;
    prop_params.sigma0_a_prop = 0.01;
    prop_params.sigma0_b_prop = 0.001;
    
    prior_params.mu0_a = 0;
    prior_params.var0_a = 100;
    prior_params.alpha0_a = 0.01;
    prior_params.beta0_a = 0.01;
    
    prior_params.mu0_b = 0;
    prior_params.var0_b = 100;
    prior_params.alpha0_b = 0.01;
    prior_params.beta0_b = 0.01;
    
    [Theta_samp_2a_ext, Dec_vec_a, Dec_vec_b] ...
        = MCMC_goc_model_2a_diff_base(Y, U, V, Z, M, Theta_init, prop_params, prior_params);
    

    D = (K1+K2+L)+1+2*N+4;
    
    %%% Plot the results
    range_U = 1:K1;
    range_V = (K1+1):(K1+K2);
    range_Z = (K1+K2+1):(K1+K2+L);
    range_sc = (K1+K2+L)+1;
    range_0a = (K1+K2+L)+1+(1:N);
    range_0b = (K1+K2+L)+1+N+(1:N);
    range_ha = (K1+K2+L)+1+2*N+(1:2);
    range_hb = (K1+K2+L)+1+2*N+(3:4);
    range_conv = (round(M/5)+1):M;
    
    legends = [legends_U, legends_V, legends_Z, 'scale parameter'];
    
    fc = fc + 1; figure(fc);
    K_max = max(K1, K2);
    for i = 1:K1+K2+L
        if i <= K1
            j = selected_factors_U(i);
            J = K_max;
        elseif i > K1 && i <= (K1 + K2)
            j = K_max + selected_factors_V(i - K1);
            J = K_max;
        elseif i > (K1+K2) && i <= (K1 + K2 + L)
            j = 2*L + (i - K1 - K2);
            J = L;
        end
        subplot(4, J, j);
        plot(Theta_samp_2a_ext(i, :)');
        title(legends{i});
        set(gca, 'xtick', [], 'xticklabel', [], 'xlim', [0, M]);
        grid on;
    end
    
    J = 7;
    subplot(4, J, 3*J + 1);
    plot(Theta_samp_2a_ext(K1+K2+L+1, :)');
    set(gca, 'xlim', [0, M]);
    title('scale parameter');
    
    subplot(4, J, 3*J + 2);
    plot(Theta_samp_2a_ext(range_0a, :)');
    set(gca, 'xlim', [0, M]);
    title('intercept');
    
    subplot(4, J, 3*J + 3);
    plot(Theta_samp_2a_ext(range_0b, :)');
    set(gca, 'xlim', [0, M]);
    title('time slope');
    
    subplot(4, J, 3*J + 4);
    plot(Theta_samp_2a_ext(range_ha(1), :)');
    set(gca, 'xlim', [0, M]);
    title('intercept: prior mean');
    
    subplot(4, J, 3*J + 5);
    plot(Theta_samp_2a_ext(range_ha(2), :)');
    set(gca, 'xlim', [0, M]);
    title('intercept: prior variance');

    subplot(4, J, 3*J + 6);
    plot(Theta_samp_2a_ext(range_hb(1), :)');
    set(gca, 'xlim', [0, M]);
    title('time slope: prior mean');
    
    subplot(4, J, 3*J + 7);
    plot(Theta_samp_2a_ext(range_hb(2), :)*1000');
    set(gca, 'xlim', [0, M]);
    title('time slope: prior variance (10^{-3})');

    print(gcf,'full_extended_model_trace_plots_from_30_May', '-depsc');
    
    fc = fc + 1; figure(fc);
    J = ceil((K1+K2+L+1)/4);
    for i = 1:K1+K2+L
        if i <= K1
            j = selected_factors_U(i);
            J = K_max;
        elseif i > K1 && i <= (K1 + K2)
            j = K_max + selected_factors_V(i - K1);
            J = K_max;
        elseif i > (K1+K2) && i <= (K1 + K2 + L+1)
            j = 2*L + (i - K1 - K2);
            J = L;
        end
        subplot(4, J, j);
        histogram(Theta_samp_2a_ext(i, range_conv)', 20,...
            'Normalization', 'probability');
        title(legends{i});
        set(gca, 'ytick', []);
        grid on;
    end
    subplot(4, 5, 3*5 + 1);
    histogram(Theta_samp_2a_ext(K1+K2+L+1, range_conv), 20,...
        'Normalization', 'probability');
    set(gca, 'ytick', []);
    title('scale parameter');
    
    subplot(4, 5, 3*5 + 2);
    histogram(Theta_samp_2a_ext(range_ha(1), range_conv), 20,...
        'Normalization', 'probability');
    set(gca, 'ytick', []);
    title('intercept: mean');
    
    subplot(4, 5, 3*5 + 3);
    histogram(Theta_samp_2a_ext(range_ha(2), range_conv),  20,...
        'Normalization', 'probability');
    set(gca, 'ytick', []);
    title('intercept: variance');
    
    subplot(4, 5, 3*5 + 4);
    histogram(Theta_samp_2a_ext(range_hb(1), range_conv),  20, ...
        'Normalization', 'probability');
    set(gca, 'ytick', []);
    title('time slope: mean');
    
    subplot(4, 5, 3*5 + 5);
    histogram(Theta_samp_2a_ext(range_hb(2), :),  20 ,'Normalization', 'probability');
    set(gca, 'ytick', []);
    title('time slope: variance');
    
    print(gcf,'full_extended_model_histograms_from_30_May', '-depsc');
    
    
    fc = fc + 1; figure(fc);
    % Sort the provinces according to the values of their base parameters
    [~, b] = sort(mean(Theta_samp_2a_ext(range_0a, range_conv), 2));
    subplot(2, 1, 1);
    boxplot(Theta_samp_2a_ext(range_0a(b), range_conv)', 'outliersize', 0.1);
    set(gca, 'xtick', 1:81, 'xTicklabel', Iller(b), 'xlim', [0 82], ...
        'XTickLabelRotation', 45, 'fontsize', 7);
    ylabel('intercept', 'fontsize', 10);
    grid on;
    
    % Sort the provinces according to the values of their base parameters
    subplot(2, 1, 2);
    [~, b] = sort(mean(Theta_samp_2a_ext(range_0b, range_conv), 2));
    boxplot(Theta_samp_2a_ext(range_0b(b), range_conv)', 'outliersize', 0.1);
    set(gca, 'xtick', 1:81, 'xTicklabel', Iller(b), 'xlim', [0 82], ...
        'XTickLabelRotation', 45, 'fontsize', 7);
    ylabel('time slope', 'fontsize', 10);
    grid on;
    print(gcf,'full_extended_intercept_and_time_slope_from_30_May', '-depsc');
    
    [~, b1] = sort(mean(Theta_samp_2a_ext(range_0a, range_conv), 2));
    [~, b2] = sort(mean(Theta_samp_2a_ext(range_0a, range_conv), 2) ...
        + mean(Theta_samp_2a_ext(range_0b, range_conv), 2)*(T-1));
    [~, b3] = sort(mean(Theta_samp_2a_ext(range_0a(b2), range_conv), 2));
    
    b4 = b2(b3);
    
    fc = fc + 1; figure(fc);
    plot([(1:1:81)' b3]');
    yyaxis left
    ylabel('initial year');
    set(gca, 'ytick', 1:81, 'yticklabel', Iller(b4), ...
        'ytickLabelRotation', 45, 'ylim', [0, 82]);
    
    yyaxis right
    ylabel('last year');
    set(gca, 'ytick', 1:81, 'yticklabel', Iller(b2), ...
        'ytickLabelRotation', 45, 'ylim', [0, 82], 'fontsize', 7);
    
    fc = fc + 1; figure(fc);
    boxplot(Theta_samp_2a_ext(1:(K1+K2+L), range_conv)', 'boxstyle', 'outline',...
        'orientation', 'horizontal', 'labels', legends(1:(K1+K2+L)), ...
        'factordirection', 'list', 'labelorientation', 'inline',...
        'outliersize', 0.1);
    grid on;
    hold on;
    plot((-9:9), (K1+0.5)*ones(1, 19), '-.k');
    plot((-9:9), (K1+K2+0.5)*ones(1, 19), '-.k');
    plot([0 0], [0 K1 + K2 + L+0.5], 'k');
    hold off;
    
    print(gcf,'full_extended_model_boxplots_from_30_May', '-depsc');
    
    % Get the mean, std, and quantiles
    m_2a_ext = zeros(1, D);
    s_2a_ext = zeros(1, D);
    p1_2a_ext = zeros(1, D);
    p2_2a_ext = zeros(1, D);
    for i = 1:D
        m_2a_ext(i) = mean(Theta_samp_2a_ext(i, range_conv));
        s_2a_ext(i) = std(Theta_samp_2a_ext(i, range_conv));
        p1_2a_ext(i) = quantile(Theta_samp_2a_ext(i, range_conv), 0.05);
        p2_2a_ext(i) = quantile(Theta_samp_2a_ext(i, range_conv), 0.95);
    end
end

%%
mig_rate  = zeros(N, T);
mig_prob_est_2a = zeros(N, T);
mig_prob_est_2a_ext = zeros(N, T);

for t = 1:T
    Theta_base_2a = Theta_samp_2a(K1+K2+L+2, range_conv) ...
        + (t-1)*Theta_samp_2a(K1+K2+L+3, range_conv);
    mig_prob_est_2a(:, t) = mean(exp(Theta_base_2a)*(N-1)...
        ./(1 + exp(Theta_base_2a)*(N-1)))*ones(N, 1);
    
    Theta_base_2a_ext = Theta_samp_2a_ext(range_0a, range_conv) ...
        + (t-1)*Theta_samp_2a_ext(range_0b, range_conv);
    mig_prob_est_2a_ext(:, t) = mean(exp(Theta_base_2a_ext)*(N-1)...
        ./(1 + exp(Theta_base_2a_ext)*(N-1)), 2);
    
    mig_rate(:, t) = (sum(Y{t}, 2) - diag(Y{t}))./sum(Y{t}, 2);
    
end

fc = fc + 1; figure(fc);
plot(mig_prob_est_2a, 'k')
hold on;
plot(mig_prob_est_2a_ext, 'r')
plot(mig_rate, 'b')
hold off;
    
