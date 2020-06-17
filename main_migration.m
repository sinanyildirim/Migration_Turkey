% This is the main code for the experiments for the paper entitled
% "A model of dynamic migration networks: Explaining Turkey's 
% inter-provincial migration flows"
% 
% Three different models are fitted to the available data:
% 1) Multinomial Logistic regression model
% 2) Dirichlet-multinomial model with a non-hierarchical specification for
% its intercept and slope parameters
% 3) Dirichlet-multinomial model with a hierarchical specification for
% its intercept and slope parameters
% 
% All the models aim to explain migration probabilities by fitting a
% certain regression function whose explanatory variables are various 
% social, economical, political, and geographical factors.
%
% The section of this file named "Algorithm specifications..." has all the
% essential options to be played with by the user, such as
% 
% - The subset of the models to be fitted to the model,
% - Numbers of iterations, and burn-in time
% - Factors (explanatory variables) to be included in the model
% - The year range for during the migration data are to used.
% 
% If you want to run an empty version of the Dirichlet-Multinomial model 
% with the hierarchical setting, set empty_run to 1.
%
% Other critical, but less managerial, choices are the proposal variances,
% and prior distributions, which are further down the script.
%
% The user is further guided through the script by comments around the
% code.
% 
% Last update: 18 June 2020

%%
clear; clc; close all; fc = 0;

seed_no = 1; rng(seed_no); % use this (change seed_no) if you want to 
% control the seed

Outputdirname = 'Outputfiles/'; % This is where the figures are saved
if ~exist(Outputdirname, 'dir')
   mkdir(Outputdirname)
end
set(0,'DefaultAxesTitleFontWeight','normal');

%% load the data from the data files
load('Datafiles/external_data', 'ilmesafe_abc_abc', 'Iller');
[Election2004, Election2009, Election2014] = import_election;
[GDP, GDP_years] = import_gdp;
[Population, Population_years] = import_population;
[Migration, Migration_years] = import_migration;
% NAN cells in the migration data are converted to 0.
Migration(isnan(Migration)) = 0;
[Unemployment, Unemployment_years] = import_unemployment;
[AKP, AKP_years] = import_AKP;
[Kurd] = import_Kurd;
[Betweenness, Betweenness_years] = import_betweenness;
[Correlations, Correlations_years] = import_correlations;

N = 81; % number of provinces

%% Algorithm specifications: Choose factors, algorithms to run, algorithm parameters, etc
algorithms_to_run = [1 1 1]; % make the corresponding index one for the 
% algorithms you want to run
M_vec = [500000 500000 500000]; % numbers of iterations
m_burn = round(M_vec/5);
run_empty = 0; % make this 1 if you want to run the empty model

Years_all = 2004:2018;
L_y = length(Years_all);
U_all = cell(1, L_y);
V_all = cell(1, L_y);
Y_all = cell(1, L_y);
Z_all = cell(1, L_y);

% Select the year range
actual_year_range = 2009:2018;
year_range = actual_year_range - (Years_all(1)-1);
T = year_range(end) - year_range(1) + 1;

one_way_factors_u = {'GDP -->', 'Unemployment -->', 'AKP -->', ...
    'Popularity -->', 'Log-population -->', 'Betweenness -->', 'Kurd -->'};
one_way_factors_v = {'GDP <--', 'Unemployment <--', 'AKP <--', ...
    'Popularity <--', 'Log-population <--', 'Betweenness <--', 'Kurd <--'};
two_way_factors = {'Political distance <-->', 'Spatial distance <-->', ...
    'Previous year <-->', 'Reciprocity <-->', 'In-flow similarity <-->'};

if run_empty == 0 % select the factors to be included here
    select_vec_U = [1 1 1 0 1 1 1];
    select_vec_V = [1 1 1 1 1 1 1];
    select_vec_Z = [1 1 1 1 1];
else % activate this for the empty model (for minimal computation)
    select_vec_U = [1 0 0 0 0 0 0];
    select_vec_V = [1 0 0 0 0 0 0];
    select_vec_Z = [1 0 0 0 0];
end

selected_factors_U = find(select_vec_U);
selected_factors_V = find(select_vec_V);

% construct legends
legends_U = one_way_factors_u(select_vec_U == 1);
legends_V = one_way_factors_v(select_vec_V == 1);
legends_Z = two_way_factors(select_vec_Z == 1);

K1 = sum(select_vec_U);
K2 = sum(select_vec_V);
L = sum(select_vec_Z);

legends{1} = [legends_U, legends_V, legends_Z, 'intercept', 'time slope'];
legends{2} = [legends_U, legends_V, legends_Z, 'scale parameter', 'intercept', 'time slope'];
legends{3} = [legends_U, legends_V, legends_Z, 'scale parameter'];

%% Factors U, V, Z, and observed migration Y are formed
for t = 1:L_y
    year = Years_all(t);
    U_all{t} = [];
    V_all{t} = [];
    
    %%%%%%%%%%%%%%%%%%% GDP  %%%%%%%%%%%%%%%%%%%
    ind = find(GDP_years == year-1);
    if select_vec_U(1) == 1
        if isempty(ind) == 0
            U_all{t} = [U_all{t} GDP(:, ind)];
        end
    end
    if select_vec_V(1) == 1
        if isempty(ind) == 0
            V_all{t} = [V_all{t} GDP(:, ind)];
        end
    end
    
    %%%%%%%%%%%%%%%%%%% Unemployment %%%%%%%%%%%%%%%%%%%
    ind = find(Unemployment_years == year-1);
    if select_vec_U(2) == 1
        if (isempty(ind) == 0)
            U_all{t} = [U_all{t} Unemployment(:, ind)];
        end
    end
    if select_vec_V(2) == 1
        if (isempty(ind) == 0)
            V_all{t} = [V_all{t} Unemployment(:, ind)];
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
    ind = find(Population_years == year-1);
    if select_vec_U(4) == 1
        if (isempty(ind) == 0)
            U_all{t} = [U_all{t} sum(Migration(:, :, ind))'];
        end
    end
    if select_vec_V(4) == 1
        if (isempty(ind) == 0)
            V_all{t} = [V_all{t} sum(Migration(:, :, ind))'];
        end
    end
    
    %%%%%%%%%%%%%%%%%%% Population %%%%%%%%%%%%%%%%%%%
    ind = find(Population_years == year-1);
    if select_vec_U(5) == 1
        if (isempty(ind) == 0)
            U_all{t} = [U_all{t} log(Population(:, ind))];
        end
    end
    if select_vec_V(5) == 1
        if (isempty(ind) == 0)
            V_all{t} = [V_all{t} log(Population(:, ind))];
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
            election_mtx = Election2004;
        elseif year >= 2009 && year < 2014
            election_mtx = Election2009;
        elseif year >= 2014
            election_mtx = Election2014;
        end
        if year >= 2004
            D = zeros(N);
            for i = 1:N
                for j = 1:(i-1)
                    D(i, j) = 0.5*sum(abs(election_mtx(i, :) - election_mtx(j, :)));
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
        ind = find(Population_years == year-1);
        if isempty(ind) == 0
            z_ind = z_ind + 1;
            Z_all{t}(:, :, z_ind) = Migration(:, :, ind);
        end
    end
    
    %%%%%%%%%%%%%%%%%%% Reciprocity %%%%%%%%%%%%%%%%%%%
    if select_vec_Z(4) == 1
        ind = find(Population_years == year-1);
        if isempty(ind) == 0
            z_ind = z_ind + 1;
            Z_all{t}(:, :, z_ind) = Migration(:, :, ind)'/max(max(Migration(:, :, ind)));
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
    ind1 = find(Population_years == year);
    ind2 = find(Population_years == year-1);
    if (isempty(ind1) == 0) && (isempty(ind2) == 0)
        Y_all{t} = Migration(:, :, ind1);
        % diagonals
        Y_all{t}((0:N-1)*N + (1:N)) = Population(:, ind2) - sum(Y_all{t}, 2);
    end
end

%% Get the inputs
Y = Y_all(year_range);
U = U_all(year_range);
V = V_all(year_range);
Z = Z_all(year_range);

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

%% MCMC for the Multinomial logistic regression model
if algorithms_to_run(1) == 1
    % Get the number of iterations
    M = M_vec(1);
    range_conv = m_burn(1)+1:M;
    
    % initial values
    Theta_init = {ones(K1, 1), ones(K2, 1), ones(L, 1), -5*ones(1, 1), 0*ones(1, 1)};
    
    % proposal parameters
    prop_params.sigma_prop = 0.015;
    prop_params.sigma0_a_prop = 0.01;
    prop_params.sigma0_b_prop = 0.01;
    
    % run the algorithm
    [Theta_samp_1] = MCMC_ML_Reg_Migr(Y, U, V, Z, M, Theta_init, prop_params);
    
    %%%%% Plot the results %%%%%
    %%%% Histograms %%%%
    fc = fc + 1; figure(fc);
    K_max = max(K1, K2);
    for i = 1:(K1+K2+L+2)
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
        histogram(Theta_samp_1(i, range_conv), 20 ,'Normalization', 'probability');
        set(gca, 'ytick', []);
        title(legends{1}{i});
        grid on;
    end
    sgtitle('Multinomial regression model - histograms');
    filenametoprint = [Outputdirname 'Multinomial_regression_model_histograms_' date];
    print(gcf, filenametoprint, '-depsc');

    %%%% Trace plots %%%%
    fc = fc + 1; figure(fc);
    K_max = max(K1, K2);
    for i = 1:(K1+K2+L+2)
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
        plot(Theta_samp_1(i, :)');
        title(legends{1}{i});
        if j <= 4*(J-1)
            set(gca, 'xtick', []);
        end
        grid on;
        
    end    
    sgtitle('Multinomial regression model - trace plots');
    filenametoprint = [Outputdirname 'Multinomial_regression_model_trace_plots_' date];
    print(gcf,filenametoprint, '-depsc');
    
    %%% Box plots %%%
    fc = fc + 1; figure(fc);
    boxplot(Theta_samp_1(1:(K1+K2+L), range_conv)', 'boxstyle', 'outline',...
        'orientation', 'horizontal', 'labels', legends{1}(1:(K1+K2+L)), ...
        'factordirection', 'list', 'labelorientation', 'inline',...
        'outliersize', 0.1);
    grid on;
    hold on;
    plot((-9:9), (K1+0.5)*ones(1, 19), '-.k');
    plot((-9:9), (K1+K2+0.5)*ones(1, 19), '-.k');
    plot([0 0], [0 K1 + K2 + L+0.5], 'k');
    hold off;
    sgtitle('Multinomial regression model - box plots');
    
    filenametoprint = [Outputdirname 'Multinomial_regression_model_boxplots_' date];
    print(gcf,filenametoprint, '-depsc');

end

%% MCMC for the Dirichlet-multinomial model with a non-hierarchical specification
if algorithms_to_run(2) == 1
    Theta_init = {zeros(K1, 1), zeros(K2, 1), zeros(L, 1), ones(1, 1), -5*ones(1, 1), zeros(1, 1)};
    M = M_vec(2);
    range_conv = m_burn(2):M;
    
    % proposal parameters
    prop_params.sigma_prop = 0.015;
    prop_params.sigma_sc_prop = 0.01;
    prop_params.sigma0_a_prop = 0.01;
    prop_params.sigma0_b_prop = 0.01;

    % run the algorithm
    [Theta_samp_2a] = MCMC_DM_Reg_Migr(Y, U, V, Z, M, Theta_init, prop_params);
    
    %%%% Plot the results %%%% 
    %%% Histograms %%% 
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
        title(legends{2}{i});
        grid on;
    end
    sgtitle('Dirichlet-multinomial model (non-hier) - histograms');
    
    filenametoprint = [Outputdirname 'Dirichlet_Multinomial_model_histograms_' date];
    print(gcf, filenametoprint, '-depsc');
    
    %%% Trace plots %%% 
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
        plot(Theta_samp_2a(i, :)');
        title(legends{2}{i});
        if j <= 4*(J-1)
            set(gca, 'xtick', []);
        end
        grid on;
        
    end
    sgtitle('Dirichlet-multinomial model (non-hier) - trace plots');
    filenametoprint = [Outputdirname 'Dirichlet_Multinomial_model_trace_plots_' date];
    print(gcf,filenametoprint, '-depsc');
    
    %%% Box plots %%% 
    fc = fc + 1; figure(fc);
    boxplot(Theta_samp_2a(1:(K1+K2+L), range_conv)', 'boxstyle', 'outline',...
        'orientation', 'horizontal', 'labels', legends{2}(1:(K1+K2+L)), ...
        'factordirection', 'list', 'labelorientation', 'inline',...
        'outliersize', 0.1);
    grid on;
    hold on;
    plot((-9:9), (K1+0.5)*ones(1, 19), '-.k');
    plot((-9:9), (K1+K2+0.5)*ones(1, 19), '-.k');
    plot([0 0], [0 K1 + K2 + L+0.5], 'k');
    hold off;
    sgtitle('Dirichlet-multinomial model (non-hier) - box plots');
    
    filenametoprint = [Outputdirname 'Dirichlet_Multinomial_model_boxplots_' date];
    print(gcf,filenametoprint, '-depsc');
    
    %%%% Collect posterior mean, std, and quantiles %%%%
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

%% MCMC for the Dirichlet-multinomial model with a hierarchical specification
if algorithms_to_run(3) == 1

    Theta_init = {zeros(K1, 1), zeros(K2, 1), zeros(L, 1), ones(1, 1), ...
        -5*ones(N, 1), zeros(N, 1), [-5; 1], [0; 1]};
        
    M = M_vec(3);
    range_conv = m_burn(3):M;
    
    % Proposal parameters
    if run_empty == 0
        prop_params.sigma_prop = 0.015;
    else
        prop_params.sigma_prop = 0;
        Theta_init{1} = zeros(K1, 1);
        Theta_init{2} = zeros(K2, 1);
        Theta_init{3} = zeros(L, 1);
    end 
    prop_params.sigma_sc_prop = 0.01;
    prop_params.sigma0_a_prop = 0.01;
    prop_params.sigma0_b_prop = 0.001;
    
    % Prior parameters
    prior_params.mu0_a = 0;
    prior_params.var0_a = 100;
    prior_params.alpha0_a = 0.01;
    prior_params.beta0_a = 0.01;
    prior_params.mu0_b = 0;
    prior_params.var0_b = 100;
    prior_params.alpha0_b = 0.01;
    prior_params.beta0_b = 0.01;
    
    % run the algorithm
    [Theta_samp_2a_ext, Dec_vec_a, Dec_vec_b] ...
        = MCMC_DM_Reg_Migr_hier(Y, U, V, Z, M, Theta_init, prop_params, prior_params);
    
    %%%% Plot the results %%%%
    D = (K1+K2+L)+1+2*N+4;
    range_U = 1:K1;
    range_V = (K1+1):(K1+K2);
    range_Z = (K1+K2+1):(K1+K2+L);
    range_sc = (K1+K2+L)+1;
    range_0a = (K1+K2+L)+1+(1:N);
    range_0b = (K1+K2+L)+1+N+(1:N);
    range_ha = (K1+K2+L)+1+2*N+(1:2);
    range_hb = (K1+K2+L)+1+2*N+(3:4);
    
    %%% Histograms %%%
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
        title(legends{3}{i});
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
    sgtitle('Dirichlet-multinomial model (hier) - histograms');
    
    filenametoprint = [Outputdirname 'Dirichlet_Multinomial_model_hier_histograms_' date];
    print(gcf,filenametoprint, '-depsc');
    
    %%% Trace plots %%%
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
        title(legends{3}{i});
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
    sgtitle('Dirichlet-multinomial model (hier) - trace plots');

    filenametoprint = [Outputdirname 'Dirichlet_Multinomial_model_hier_trace_plots_' date];
    print(gcf,filenametoprint, '-depsc');
    
    %%% Box plots %%%
    fc = fc + 1; figure(fc);
    boxplot(Theta_samp_2a_ext(1:(K1+K2+L), range_conv)', 'boxstyle', 'outline',...
        'orientation', 'horizontal', 'labels', legends{3}(1:(K1+K2+L)), ...
        'factordirection', 'list', 'labelorientation', 'inline',...
        'outliersize', 0.1);
    grid on;
    hold on;
    plot((-9:9), (K1+0.5)*ones(1, 19), '-.k');
    plot((-9:9), (K1+K2+0.5)*ones(1, 19), '-.k');
    plot([0 0], [0 K1 + K2 + L+0.5], 'k');
    hold off;
    
    sgtitle('Dirichlet-multinomial model (hier) - box plots');
    filenametoprint = [Outputdirname 'Dirichlet_Multinomial_model_hier_boxplots_' date];
    print(gcf,filenametoprint, '-depsc');
    
    %%% Box plots for the intercept and slope parameters %%%
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

    filenametoprint = [Outputdirname 'Dirichlet_Multinomial_model_hier_intercept_and_time_slope_' date];
    print(gcf,filenametoprint, '-depsc');

    %%% Intercept + (t-1)*slope at the initial and last years %%%
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
    
    filenametoprint = [Outputdirname 'Dirichlet_Multinomial_model_hier_init_and_last_years_' date];
    print(gcf,filenametoprint, '-depsc');
    
    %%% Get the posterior mean, std, and quantiles %%%
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
