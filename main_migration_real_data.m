% This is the main code for the real-data experiments for the paper entitled
% "A model of dynamic migration networks: Explaining Turkey's
% inter-provincial migration flows"
%
% A Dirichlet-multinomial model and a gravity model is fitted to the data 
%
% The models aim to explain migration probabilities (or counts) by fitting 
% a certain regression function whose explanatory variables are various
% social, economical, political, and geographical factors.
%
% The section of this file named "Algorithm specifications..." has all the
% essential options to be played with by the user, such as
%
% - Numbers of iterations, and burn-in time
% - Factors (explanatory variables) to be included in the model
% - The year range for during the migration data are to used.
%
% If you want to run an empty version of the Dirichlet-Multinomial model
% with the hierarchical setting, set empty_run to 1.
%
% The user is further guided through the script by comments around the
% code.
%
% Sinan Yildirim
% Last update: 22 April 2023

%% Preamble
clear; clc; close all; fc = 0;

seed_no = 1; rng(seed_no); % use this (change seed_no) if you want to
% control the seed

Outputdirname = 'Outputfiles/'; % This is where the output files are saved
if ~exist(Outputdirname, 'dir')
    mkdir(Outputdirname)
end
Outputfigsdirname = 'Outputfigures/'; % This is where the figures are saved
if ~exist(Outputfigsdirname, 'dir')
    mkdir(Outputfigsdirname)
end

set(0,'DefaultAxesTitleFontWeight','normal');

full_or_partial = {'wopast_model', 'full_model'};

%% load the data from the data files
load('Datafiles/migration_data_Turkey');

%% Show posterior results?
show_post_results = 1;

%% Algorithm specifications: Choose factors, algorithms to run, algorithm parameters, etc
N = 81; % number of provinces
M = 500000; % numbers of iterations
m_burn = round(M/2);
range_conv = m_burn:M;
range_conv_pred = (9*M/10+1):M;
M_burn_pred = length(range_conv_pred);

run_empty = 0; % make this 1 if you want to run the empty model
regress_on_past = 1; % make this 1 if you want to regress on the previous observations
years_all = 2004:2018;
K0 = 1; % (order+1) of the polynomial for the sender base parameters of the provinces
theta0_common = 0;
d0 = (theta0_common == 0)*N + (theta0_common == 1)*1;

% number of sender and receiver base parameters (1 for common, 81 for province-based)

% Select the year range
years_data = 2009:2018;
% select the prediction years
test_years_vec = 2009:2018; Years_test_vec = num2cell(test_years_vec);
% choose the one below to run MCMC on the whole dataset
% test_years_vec = []; Years_test_vec = {[]};

T = length(years_data);

one_way_factors = {'GDP', 'Unemployment', 'AKP', ...
    'Popularity', 'Population', 'Betweenness', 'Kurd'};
two_way_factors = {'Political distance', 'Spatial distance', ...
    'Previous year', 'Reciprocity', 'In-flow similarity'};

one_way_factors_u = {'GDP -->', 'Unemployment -->', 'AKP -->', ...
    'Popularity -->', 'Population -->', 'Betweenness -->', 'Kurd -->'};
one_way_factors_v = {'GDP <--', 'Unemployment <--', 'AKP <--', ...
    'Popularity <--', 'Population <--', 'Betweenness <--', 'Kurd <--'};
two_way_factors_z = {'Political distance <-->', 'Spatial distance <-->', ...
    'Previous year <-->', 'Reciprocity <-->', 'In-flow similarity <-->'};

if regress_on_past == 1
    % select the factors to be included here
    select_vec_U = [1 1 1 1 1 1 1];
    select_vec_V = [1 1 1 1 1 1 1];
    select_vec_Z = [1 1 1 1 1];
else
    select_vec_U = [1 1 1 0 1 0 1];
    select_vec_V = [1 1 1 0 1 0 1];
    select_vec_Z = [1 1 0 0 0];
end

filename_extension = ['_U_' sprintf('%d', select_vec_U) ...
    '_V_' sprintf('%d', select_vec_V) ...
    '_Z_' sprintf('%d', select_vec_Z)...
    sprintf('_K0_%d_b_%d_M_%d', K0, d0, M)...
    '_test' sprintf('_%d', test_years_vec)];

K1 = sum(select_vec_U);
K2 = sum(select_vec_V);
L = sum(select_vec_Z);


%% Construct variables according to the choices above
L_y = length(years_all);
years_all_str = cell(1, L_y);
for i = 1:L_y
    years_all_str{i} = num2str(years_all(i));
end

data_year_inds = zeros(1, T);
for t = 1:T
    data_year_inds(t) = find(years_all == years_data(t));
end
data_year_names = years_all_str(data_year_inds);

selected_factors_U = find(select_vec_U);
selected_factors_V = find(select_vec_V);

% construct legends
legends_U0 = one_way_factors(select_vec_U == 1);
legends_V0 = one_way_factors(select_vec_V == 1);
legends_Z0 = two_way_factors(select_vec_Z == 1);
legends_U = one_way_factors_u(select_vec_U == 1);
legends_V = one_way_factors_v(select_vec_V == 1);
legends_Z = two_way_factors_z(select_vec_Z == 1);

legends{1} = [legends_U, legends_V, legends_Z, 'intercept', 'time slope'];
legends{2} = [legends_U, legends_V, legends_Z, 'scale parameter', 'intercept', 'time slope'];
legends{3} = [legends_U, legends_V, legends_Z, 'scale parameter'];

%% Factors U, V, Z, and observed migration Y are formed
U_all = cell(1, L_y);
V_all = cell(1, L_y);
Y_all = cell(1, L_y);
Z_all = cell(1, L_y);

for t = 1:L_y
    year = years_all(t);
    U_all{t} = [];
    V_all{t} = [];

    %%%%%%%%%%%%%%%%%%% GDP  %%%%%%%%%%%%%%%%%%%
    ind = find(GDP_years == year-1);
    if select_vec_U(1) == 1  && isempty(ind) == 0
        U_all{t} = [U_all{t} GDP(:, ind)];
    end
    if select_vec_V(1) == 1 && isempty(ind) == 0
        V_all{t} = [V_all{t} GDP(:, ind)];
    end

    %%%%%%%%%%%%%%%%%%% Unemployment %%%%%%%%%%%%%%%%%%%
    ind = find(Unemployment_years == year-1);
    if select_vec_U(2) == 1 && (isempty(ind) == 0)
        U_all{t} = [U_all{t} Unemployment(:, ind)];
    end
    if select_vec_V(2) == 1 && isempty(ind) == 0
        V_all{t} = [V_all{t} Unemployment(:, ind)];
    end

    %%%%%%%%%%%%%%%%%%% AKP  %%%%%%%%%%%%%%%%%%%
    ind = find(AKP_years == year-1);
    if select_vec_U(3) == 1 && isempty(ind) == 0
        U_all{t} = [U_all{t} AKP(:, ind)];
    end
    if select_vec_V(3) == 1 && isempty(ind) == 0
        V_all{t} = [V_all{t} AKP(:, ind)];
    end

    %%%%%%%%%%%%%%%%%%% Popularity %%%%%%%%%%%%%%%%%%%
    ind = find(Population_years == year-1);
    if select_vec_U(4) == 1 && isempty(ind) == 0
        U_all{t} = [U_all{t} log(max(1, sum(Migration(:, :, ind))))'];
    end
    if select_vec_V(4) == 1 && isempty(ind) == 0
        V_all{t} = [V_all{t} log(max(1, sum(Migration(:, :, ind))))'];
    end

    %%%%%%%%%%%%%%%%%%% Population %%%%%%%%%%%%%%%%%%%
    ind = find(Population_years == year-1);
    if select_vec_U(5) == 1 && isempty(ind) == 0
        U_all{t} = [U_all{t} log(Population(:, ind))];
    end
    if select_vec_V(5) == 1 && isempty(ind) == 0
        V_all{t} = [V_all{t} log(Population(:, ind))];
    end

    %%%%%%%%%%%%%%%%%%% Betweenness %%%%%%%%%%%%%%%%%%%
    ind = find(Betweenness_years == year-1);
    if select_vec_U(6) == 1 && isempty(ind) == 0
        U_all{t} = [U_all{t} Betweenness(:, ind)];
    end
    if select_vec_V(6) == 1 && isempty(ind) == 0
        V_all{t} = [V_all{t} Betweenness(:, ind)];
    end

    % Construct Z{t}
    z_ind = 0;

    %%%%%%%%%%%%%%%%%%% Kurd  %%%%%%%%%%%%%%%%%%%
    if select_vec_U(7) == 1
        if year < 2013
            U_all{t} = [U_all{t} Kurds(:, 2)];
        else
            U_all{t} = [U_all{t} Kurds(:, 1)];
        end
    end
    if select_vec_V(7) == 1
        if year < 2013
            V_all{t} = [V_all{t} Kurds(:, 2)];
        else
            V_all{t} = [V_all{t} Kurds(:, 1)];
        end
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
            DD = zeros(N);
            for i = 1:N
                for j = 1:(i-1)
                    DD(i, j) = 0.5*sum(abs(election_mtx(i, :) - election_mtx(j, :)));
                end
            end
            DD = (DD + DD')/100;
            Z_all{t}(:, :, z_ind) = DD;
        end
    end

    %%%%%%%%%%%%%%%%%%% Geographical distance %%%%%%%%%%%%%%%%%%%
    if select_vec_Z(2) == 1
        z_ind = z_ind + 1;
        Z_all{t}(:, :, z_ind) = ilmesafe_abc_abc;
    end

    %%%%%%%%%%%%%%%%%%% Previous year %%%%%%%%%%%%%%%%%%%
    ind = find(Population_years == year-1);
    if select_vec_Z(3) == 1 && isempty(ind) == 0 
        z_ind = z_ind + 1;
        Z_all{t}(:, :, z_ind) = log(max(1, Migration(:, :, ind)));
    end

    %%%%%%%%%%%%%%%%%%% Reciprocity %%%%%%%%%%%%%%%%%%%
    ind = find(Population_years == year-1);
    if select_vec_Z(4) == 1 && isempty(ind) == 0
        z_ind = z_ind + 1;
        Z_all{t}(:, :, z_ind) = log(max(1, Migration(:, :, ind)))' ...
            - log(max(max(Migration(:, :, ind))));
    end

    %%%%%%%%%%%%%%%%%%% Correaltions %%%%%%%%%%%%%%%%%%%
    ind = find(Correlations_years == year-1);
    if select_vec_Z(4) == 1  && isempty(ind) == 0
        z_ind = z_ind + 1;
        Z_all{t}(:, :, z_ind) = Correlations(:, :, ind);
    end

    %%%%%%%%%%%%%%%%%%% Construct Y %%%%%%%%%%%%%%%%%%%
    ind1 = find(Population_years == year);
    ind2 = find(Population_years == year-1);
    if (isempty(ind1) == 0) && (isempty(ind2) == 0)
        Y_all{t} = Migration(:, :, ind1);
        % diagonals
        Y_all{t}((0:N-1)*N + (1:N)) = Population(:, ind2) - sum(Y_all{t}, 2);
        % correction due to birth and death
        birth_minus_death = Population(:, ind1) - sum(Migration(:, :, ind1), 2) ...
            + sum(Migration(:, :, ind1), 1)' - Population(:, ind2);
        Y_all{t}((0:N-1)*N + (1:N)) = Y_all{t}((0:N-1)*N + (1:N)) + birth_minus_death';
    end
end

%% Get the inputs
Y = Y_all(data_year_inds);
U = U_all(data_year_inds);
V = V_all(data_year_inds);
Z = Z_all(data_year_inds);

%% normalise U's, V's, and Z's
max_abs_U = zeros(1, K1); mean_U = zeros(1, K1);
max_abs_V = zeros(1, K2); mean_V = zeros(1, K2);
max_abs_Z = zeros(1, L);  mean_Z = zeros(1, L);

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

%% Prepare X and Y for the gravity model
[X_mod, Y_mod] = convert_data_DM2grav(U, V, Z, Y, K0, theta0_common);

%% Initalize the error metrics for the algorithms
Num_of_iters_pred = length(Years_test_vec);
% Get the dimensions of the arrays of the raw errors
T_tests_vec = zeros(1, Num_of_iters_pred);
for iter_pred = 1:Num_of_iters_pred
    T_tests_vec(iter_pred) = length(Years_test_vec{iter_pred});
end
T_test_cumsum = [0 cumsum(T_tests_vec)];
T_test_sum = sum(T_tests_vec);

diag_ind = (0:N-1)*N + (1:N);
all_ind = 1:N^2;
non_diag_ind = setdiff(all_ind, diag_ind);
inds_to_compare_cell = {non_diag_ind, all_ind};
L_inds = [length(inds_to_compare_cell{1}) length(inds_to_compare_cell{2})];

% Initialize the error metrics
temp_cell = repmat({repmat({zeros(1, Num_of_iters_pred)}, 3, 2)}, 1, 2);
MSE = repmat({repmat({zeros(1, Num_of_iters_pred)}, 3, 2)}, 1, 2);

% Initialise the arrays of all errors
Temp_mtx1 = repmat({zeros(L_inds(1), T_test_sum)}, 3, 1);
Temp_mtx2 = repmat({zeros(L_inds(2), T_test_sum)}, 3, 1);
Temp_mtx_SE = repmat({[Temp_mtx1 Temp_mtx2]}, 1, 2);

% initialize the arrays of true and predicted values
All_true_vals = cell(1, Num_of_iters_pred);
All_pred_vals = repmat({cell(1, Num_of_iters_pred)}, 1, 2);

%% Training and estimation for every training-testing data combination
for iter_pred = 1:Num_of_iters_pred

    disp(iter_pred);
    years_test = Years_test_vec{iter_pred};

    % get the data (mostly indices) for the training and test years
    years_train = setdiff(years_data, years_test);
    T_train = length(years_train);
    T_test = length(years_test);

    train_year_inds = zeros(1, T_train);
    test_year_inds = zeros(1, T_test);
    train_sub_t = zeros(1, T_train);
    test_sub_t = zeros(1, T_test);

    for t = 1:T_train
        train_year_inds(t) = find(years_all == years_train(t));
        train_sub_t(t) = find(years_data == years_train(t));
    end
    for t = 1:T_test
        test_year_inds(t) = find(years_all == years_test(t));
        test_sub_t(t) = find(years_data == years_test(t));
    end
    train_year_names = years_all_str(train_year_inds);
    test_year_names = years_all_str(test_year_inds);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%% Method 1: MCMC %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Train with MCMC
    Y_train = Y(train_sub_t);
    U_train = U(train_sub_t);
    V_train = V(train_sub_t);
    Z_train = Z(train_sub_t);

    load('theta_init_migTR'); % this loads theta_init
    [Theta_samp, range_theta] = MCMC_Migration(Y_train, U_train, V_train, Z_train, ...
        train_sub_t, K0, theta0_common, M, theta_init);
    
    %%% Predict with MCMC
    U_test = U(test_sub_t);
    V_test = V(test_sub_t);
    Z_test = Z(test_sub_t);
    Y_test = Y(test_sub_t);

    Prop_test = cell(1, T_test);
    log_Y_test = cell(1, T_test);
    P_test = cell(1, T_test);
    
    for t = 1:T_test
        log_Y_test{t} = log(max(Y_test{t}, 1));
        Prop_test{t} = Y_test{t}./sum(Y_test{t}, 2);
        P_test{t} = sum(Y_test{t}, 2);
    end

    Theta_for_pred = Theta_samp(:, range_conv_pred);
    if T_test > 0
        [Y_pred_MCMC, log_Y_pred_MCMC, Prop_pred_MCMC] = ...
            pred_migration(P_test, U_test, V_test, Z_test, test_sub_t, Theta_for_pred, ...
            range_theta);
    end

    %% %%%%%%%%%%%%%%%%%%% Method 2: Gravity model %%%%%%%%%%%%%%%%%%%%%%
    % determine the indices of test and training data
    rows_per_year = N*(N-1);
    n = T*rows_per_year;
    test_inds = zeros(1, T_test*rows_per_year);
    for t = 1:T_test
        t_test = test_sub_t(t);
        test_ind_temp = (t-1)*rows_per_year+1:(rows_per_year*t);
        test_inds(test_ind_temp) = ((t_test-1)*rows_per_year+1):t_test*rows_per_year;
    end
    train_inds = setdiff(1:n, test_inds);

    % get the X, y data
    X_mod_train = X_mod(train_inds, :);
    X_mod_test = X_mod(test_inds, :);
    Y_mod_train = Y_mod(train_inds);
    Y_mod_test = Y_mod(test_inds);

    % perform estimation
    beta_est = (X_mod_train'*X_mod_train)\(X_mod_train'*Y_mod_train);

    % perform prediction
    log_y_pred = X_mod_test*beta_est;
    y_pred = exp(log_y_pred);

    % store the predictions
    Prop_pred_grav = repmat({zeros(N)}, 1, T_test);
    Y_pred_grav = repmat({zeros(N)}, 1, T_test);
    log_Y_pred_grav = repmat({zeros(N)}, 1, T_test);

    for t = 1:T_test
        test_ind_temp = (t-1)*rows_per_year+1:(rows_per_year*t);
        log_Y_pred_grav{t}(inds_to_compare_cell{1}) = log_y_pred(test_ind_temp);
        Y_pred_grav{t}(inds_to_compare_cell{1}) = y_pred(test_ind_temp);

        % take the transposes to put in the right form
        log_Y_pred_grav{t} = log_Y_pred_grav{t}';
        Y_pred_grav{t} = Y_pred_grav{t}';
        
        % fill in the diagonal terms
        y_pred_diag_temp = sum(Y_test{t}, 2) - sum(Y_pred_grav{t}, 2);
        Y_pred_grav{t}(diag_ind) = y_pred_diag_temp;
        log_Y_pred_grav{t}(diag_ind) = log(max(1, y_pred_diag_temp));

        % proportions
        Prop_pred_grav{t} = Y_pred_grav{t}./sum(Y_test{t}, 2);
    end
    
    % Store the predictions and the true values
    All_pred_vals{1}{iter_pred} = {Y_pred_MCMC, log_Y_pred_MCMC, Prop_pred_MCMC};
    All_pred_vals{2}{iter_pred} = {Y_pred_grav, log_Y_pred_grav, Prop_pred_grav};
    All_true_vals{iter_pred} = {Y_test, log_Y_test, Prop_test};

    %% %%%%%%%%%%%%%%%%%%%%%%% Error metrics %%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1:3 % type of prediction (log, linear, proportion)
        for j = 1:2 % non-diagonal and all
            inds_to_compare = inds_to_compare_cell{j};
            for algo_type = [1 2] % MCMC and gravity
                Temp_vec_SE = zeros(L_inds(j), T_test);
                for t = 1:T_test
                    pred_vals = All_pred_vals{algo_type}{iter_pred}{i}{t}(inds_to_compare);
                    true_vals = All_true_vals{iter_pred}{i}{t}(inds_to_compare);
                    temp_err = pred_vals - true_vals;
                    Temp_vec_SE(:, t) = temp_err.^2;
                end
   
                % calculate the mean and the median of the errors
                MSE{algo_type}{i, j}(iter_pred) = mean(Temp_vec_SE(:));

                % These are the arrays of raw errrors and they necessary for the median
                t_indices = T_test_cumsum(iter_pred)+1:T_test_cumsum(iter_pred+1);
                Temp_mtx_SE{algo_type}{i, j}(:, t_indices) = Temp_vec_SE;
            end
        end
    end
end

%% Overall testing results
MSE_total = repmat({zeros(3, 2)}, 1, 2);

for i = 1:3 % type of prediction (log, linear, proportion)
    for j = 1:2 % non-diagonal and all
        for algo_type = [1 2] % MCMC and gravity
            % Calculate the mean and median of all errors
            MSE_total{algo_type}(i, j) = mean(Temp_mtx_SE{algo_type}{i, j}(:));
        end
    end
end

%% Show and plot results
if show_post_results == 1
    % prediction results
    i_plot = 1; j_plot = 1;
    fc = fc + 1; figure(fc); % new figure for each algorithm

    for algo_type = [1 2]
        subplot(1, 2, algo_type);
        hold on;
        for iter_pred = 1:Num_of_iters_pred
            for t = 1:T_test
                inds_to_compare = inds_to_compare_cell{j_plot};
                pred_vals = All_pred_vals{algo_type}{iter_pred}{i_plot}{t}(inds_to_compare);
                true_vals = All_true_vals{iter_pred}{i_plot}{t}(inds_to_compare);
                plot(pred_vals, true_vals, '.k');
            end
        end
        hold off;
        % set(gca, 'xlim', [0, 20], 'ylim', [0, 20]);
        xlabel('predicted values');
        ylabel('true values');
        if algo_type == 1
            title('Dirichlet-Multinomial');
        elseif algo_type == 2
            title('Gravity');
        end

    end
    
    %% Correlation between the factors
    fc = fc + 1; figure(fc);
    UU = zeros(T*N, K1);
    VV = zeros(T*N, K2);
    ZZ = zeros(T*N^2, L);
    
    for t = 1:T
        for k = 1:K1, UU((t-1)*N+1:t*N, k) = U{t}(:, k); end
        for k = 1:K2, VV((t-1)*N+1:t*N, k) = V{t}(:, k); end
        for l = 1:L
            temp_mtx = Z{t}(:, :, l);
            ZZ((t-1)*N^2+1:t*N^2, l) = temp_mtx(:);
        end
    end
    C = corr(UU);
    imagesc(C); colormap(gray); colorbar; caxis([-1, 1]);
    % title('correlation matrix between one-way factors');
    set(gca, 'xtick', 1:K1, 'xticklabel', legends_U0);
    set(gca, 'ytick', 1:K1, 'yticklabel', legends_U0);
    str_t = cell(K1);

    for i = 1:K1
        for j = 1:K1
            str_t{i, j} = sprintf('%.2f', C(i, j));
        end
    end
    x = repmat(1:K1,K1,1); % generate x-coordinates
    y = x'; % generate y-coordinates
    text(x(:), y(:), str_t, 'HorizontalAlignment', 'Center', 'color','k');
    filenametoprint = [Outputfigsdirname 'factor_corr_oneway_U' filename_extension];
    print(gcf,'-depsc', filenametoprint);

    C = corr(VV);
    fc = fc + 1; figure(fc);
    imagesc(C); colormap(gray); colorbar; caxis([-1, 1]);
    % title('correlation matrix between one-way factors');
    set(gca, 'xtick', 1:K2, 'xticklabel', legends_V0);
    set(gca, 'ytick', 1:K2, 'yticklabel', legends_V0);
    str_t = cell(K2);
    for i = 1:K2
        for j = 1:K2
            str_t{i, j} = sprintf('%.2f', C(i, j));
        end
    end
    x = repmat(1:K2,K2,1); % generate x-coordinates
    y = x'; % generate y-coordinates
    text(x(:), y(:), str_t, 'HorizontalAlignment', 'Center', 'color','k');
    filenametoprint = [Outputfigsdirname 'factor_corr_oneway_V' filename_extension];
    print(gcf,'-depsc', filenametoprint);

    fc = fc + 1; figure(fc);
    C = corr(ZZ);
    imagesc(C); colormap(gray); colorbar; caxis([-1, 1]);
    % title('correlation matrix between two-way factors');
    set(gca, 'xtick', 1:L, 'xticklabel', legends_Z0);
    set(gca, 'ytick', 1:L, 'yticklabel', legends_Z0);
    str_t = cell(L);

    for i = 1:L
        for j = 1:L
            str_t{i, j} = sprintf('%.2f', C(i, j));
        end
    end
    x = repmat(1:L,L,1); % generate x-coordinates
    y = x'; % generate y-coordinates
    text(x(:), y(:), str_t, 'HorizontalAlignment', 'Center', 'color','k');
    filenametoprint = [Outputfigsdirname 'factor_corr_twoways' filename_extension];
    print(gcf,'-depsc', filenametoprint);

    clear C x y t;

    %% Posterior correlation
    fc = fc + 1; figure(fc);
    D_f = K1 + K2 + L;
    C = corr(Theta_samp(1:D_f, range_conv)');
    imagesc(C); colormap(gray); colorbar; caxis([-1, 1]);
    % title('posterior correlation matrix of $$\theta_{1} \, \theta_{2} \, \theta_{3}$$', 'Interpreter', 'Latex', 'fontsize', 18);
    set(gca, 'xtick', 1:D_f, 'xticklabel', [legends_U legends_V legends_Z]);
    set(gca, 'ytick', 1:D_f, 'yticklabel', [legends_U legends_V legends_Z]);
    str_t = cell(D_f);

    for i = 1:D_f
        for j = 1:D_f
            str_t{i, j} = sprintf('%.2f', C(i, j));
        end
    end
    x = repmat(1:D_f,D_f,1); % generate x-coordinates
    y = x'; % generate y-coordinates
    text(x(:), y(:), str_t, 'HorizontalAlignment', 'Center', 'color', [0.7, 0, 0]);
    hold on;
    plot([K1+0.5 K1+0.5], [0 D_f+0.5], 'b', 'linewidth', 2);
    plot([K1+K2+0.5 K1+K2+0.5], [0 D_f+0.5], 'b', 'linewidth', 2);
    plot([0 D_f+0.5], [K1+0.5 K1+0.5], 'b', 'linewidth', 2);
    plot([0 D_f+0.5], [K1+K2+0.5 K1+K2+0.5], 'b', 'linewidth', 2);

    hold off;
    filenametoprint = [Outputfigsdirname 'DM_post_corr' filename_extension];
    print(gcf,'-depsc', filenametoprint);

    %% %%%%%%%%%%%%%%%%%%%%%%%% Histograms and Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%% Histogram  of -->, <--, and <--> parameters %%%%%%%%%%%%
    fc = fc + 1; figure(fc);
    for i = 1:K1+K2+L
        if i <= K1
            J = K1;
            j = i;
        elseif i > K1 && i <= (K1 + K2)
            J = K2;
            j = K2 + (i - K1);
        elseif i > (K1+K2) && i <= (K1 + K2 + L+1)
            J = L;
            j = 2*L + (i - K1 - K2);
        end
        subplot(3, J, j);
        histogram(Theta_samp(i, range_conv)', 20, 'Normalization', 'probability');
        title(legends{3}{i});
        set(gca, 'ytick', []);
        grid on;
    end
    % sgtitle('Dirichlet-multinomial model - histograms');
    filenametoprint = [Outputfigsdirname 'DM_histograms' filename_extension];
    print(gcf,filenametoprint, '-depsc');

    %%%%%%%%%%%% Histogram  of sender base hyparparameters %%%%%%%%%%%%
    ylabels = {'intercept', 'time slope', 'quad coeff'};
    fc = fc + 1; figure(fc);
    for i = 1:K0
        subplot(K0+1, K0, i);
        histogram(Theta_samp(range_theta{6}(i), range_conv)', 20, 'Normalization', 'probability');
        set(gca, 'ytick', []);
        title(sprintf(' mean (%d)', i));
        for j = i:K0
            subplot(K0+1, K0, K0+(i-1)*K0+j);
            histogram(Theta_samp(range_theta{7}((i-1)*K0+j), range_conv)', 20,...
                'Normalization', 'probability');
            set(gca, 'ytick', []);
            title(sprintf('cov. (%d, %d)', i, j));
        end
    end
    % sgtitle('Dirichlet-multinomial model - base mean and covariance histograms');
    filenametoprint = [Outputfigsdirname 'DM_base_mean_var_hist' filename_extension];
    print(gcf,filenametoprint, '-depsc');

    %%%%%%%%%%%% Histogram  of scale parameter %%%%%%%%%%%%
    fc = fc + 1; figure(fc);

    histogram(Theta_samp(range_theta{4}, range_conv)', 20, 'Normalization', 'probability');
    % title('Dirichlet-multinomial model - scale parameter histogram');
    set(gca, 'ytick', []);
    filenametoprint = [Outputfigsdirname 'DM_scale_hist' ...
        filename_extension];
    print(gcf,filenametoprint, '-depsc');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%% Trace plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Trace plots of -->, <--, and <--> parameters %%%%%%%%%%%%
    fc = fc + 1; figure(fc);
    K_max = max(K1, K2);
    for i = 1:K1+K2+L
        if i <= K1
            J = K1;
            j = i;
        elseif i > K1 && i <= (K1 + K2)
            J = K2;
            j = K2 + (i - K1);
        elseif i > (K1+K2) && i <= (K1 + K2 + L+1)
            J = L;
            j = 2*L + (i - K1 - K2);
        end
        subplot(3, J, j);
        plot(Theta_samp(i, :)');
        title(legends{3}{i});
        set(gca, 'xtick', [], 'xticklabel', [], 'xlim', [0, M]);
        grid on;
    end
    % sgtitle('Dirichlet-multinomial model - trace plots');
    filenametoprint = [Outputfigsdirname 'DM_trace_plots' filename_extension];
    print(gcf,filenametoprint, '-depsc');

    %%% Trace plots of sender base hyperparameters %%%%%%%%%%%%
    fc = fc + 1; figure(fc);
    for i = 1:K0
        subplot(K0+1, K0, i);
        plot(Theta_samp(range_theta{6}(i), range_conv)');
        title(sprintf(' mean (%d)', i));
        for j = i:K0
            subplot(K0+1, K0, K0+(i-1)*K0+j);
            plot(Theta_samp(range_theta{7}((i-1)*K0+j), range_conv)');
            % set(gca, 'ytick', []);
            title(sprintf('cov. (%d, %d)', i, j));
        end
    end
    % sgtitle('Dirichlet-multinomial model - base sender mean and covariance trace plots');
    filenametoprint = [Outputfigsdirname 'DM_hier_base_sender_mean_var_trace_plots' filename_extension];
    print(gcf,filenametoprint, '-depsc');

    %%% Trace plots of scale parameters %%%%%%%%%%%%
    fc = fc + 1; figure(fc);
    plot(Theta_samp(range_theta{4}, range_conv));
    set(gca, 'xticklabel', []);
    % title('Dirichlet-multinomial model - scale parameter trace plot');
    filenametoprint = [Outputfigsdirname 'DM_scale_trace_plots' ...
        filename_extension];
    print(gcf,filenametoprint, '-depsc');

    %% %%%%%%%%%%%%%%%%%%%%%%%%%% Box plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Boxplots of -->, <--, and <--> parameters %%%%%%%%%%%%
    fc = fc + 1; figure(fc);
    boxplot(Theta_samp(1:(K1+K2+L), range_conv)', 'boxstyle', 'outline',...
        'orientation', 'horizontal', 'labels', legends{3}(1:(K1+K2+L)), ...
        'factordirection', 'list', 'labelorientation', 'inline',...
        'outliersize', 0.1);
    grid on;
    hold on;
    plot((-9:9), (K1+0.5)*ones(1, 19), '-.k');
    plot((-9:9), (K1+K2+0.5)*ones(1, 19), '-.k');
    plot([0 0], [0 K1 + K2 + L+0.5], 'k');
    hold off;

    % sgtitle('Dirichlet-Multinomial model - box plots');
    filenametoprint = [Outputfigsdirname 'DM_hier_boxplots' filename_extension];
    print(gcf,filenametoprint, '-depsc');

    %% Boxplots of -->, <--, and <--> parameters separately %%%%%%%%%%%%
    % fc = fc + 1; figure(fc);
    subplot(1, 3, 1);
    boxplot(Theta_samp(1:(K1), range_conv)', 'boxstyle', 'outline',...
        'orientation', 'horizontal', 'labels', legends{3}(1:K1), ...
        'factordirection', 'list', 'outliersize', 0.1);
    hold on;
    plot([0 0], [0 K1+0.5], 'k');
    grid on;
    hold off;

    subplot(1, 3, 2);
    boxplot(Theta_samp(K1+1:(K1+K2), range_conv)', 'boxstyle', 'outline',...
        'orientation', 'horizontal', 'labels', legends{3}(K1+1:(K1+K2)), ...
        'factordirection', 'list', 'outliersize', 0.1);    
    hold on;
    plot([0 0], [0 K2+0.5], 'k');
    grid on;
    hold off;

    subplot(1, 3, 3);
    boxplot(Theta_samp((K1+K2+1):(K1+K2+L), range_conv)', 'boxstyle', 'outline',...
        'orientation', 'horizontal', 'labels', legends{3}(K1+K2+1:(K1+K2+L)), ...
        'factordirection', 'list', 'outliersize', 0.1);
    hold on;
    plot([0 0], [0 L+0.5], 'k');
    grid on;
    hold off;

    % sgtitle('Dirichlet-Multinomial model - box plots');
    filenametoprint = [Outputfigsdirname 'DM_hier_boxplots_separate' filename_extension];
    print(gcf,filenametoprint, '-depsc');

    %% Box plots for the intercept and slope parameters - sender %%%
    if theta0_common == 0
        Range_poly = reshape(range_theta{5}, K0, N);
        temp_vec = mean(Theta_samp(range_theta{5}, range_conv), 2);
        Poly_coeff_mtx = reshape(temp_vec, K0, N);

        Base_mtx = zeros(T, N);
        for t = 1:T
            t_vec = (t-1).^(0:(K0-1));
            Base_mtx(t, :) = t_vec*Poly_coeff_mtx;
        end

        fc = fc + 1; figure(fc);
        ylabels = {'intercept', 'time slope', 'quad coeff'};
        for k = 1:K0
            % Sort the provinces according to the values of their base parameters
            [~, b] = sort(Poly_coeff_mtx(k, :));
            subplot(K0, 1, k);
            boxplot(Theta_samp(Range_poly(k, b), range_conv)', 'outliersize', 0.1);
            set(gca, 'xtick', 1:81, 'xTicklabel', Iller(b), 'xlim', [0 82], ...
                'XTickLabelRotation', 45, 'fontsize', 8);
            ylabel(ylabels{k}, 'fontsize', 10);
            grid on;
        end

        filenametoprint = [Outputfigsdirname ...
            'DM_intercept_and_time_slope_sender_'...
            filename_extension];
        print(gcf,filenametoprint, '-depsc');

        %%% Intercept + (t-1)*slope at the initial and last years %%%
        % base parameters - beginning:

        [~, b1] = sort(Base_mtx(1, :));
        [~, b2] = sort(Base_mtx(T, :));
        [~, b3] = sort(Base_mtx(1, b2));
        b4 = b2(b3);

        fc = fc + 1; figure(fc);
        plot([(1:1:81)' b3']');
        yyaxis left
        ylabel('initial year');
        set(gca, 'ytick', 1:81, 'yticklabel', Iller(b4), ...
            'ytickLabelRotation', 45, 'ylim', [0, 82]);

        yyaxis right
        ylabel('last year');
        set(gca, 'ytick', 1:81, 'yticklabel', Iller(b2), ...
            'ytickLabelRotation', 45, 'ylim', [0, 82], 'fontsize', 7);

        filenametoprint = [Outputfigsdirname ...
            'DM_init_and_last_years_sender_' filename_extension];
        print(gcf,filenametoprint, '-depsc');
    end

    % Get the posterior mean, std, and quantiles %%%
    D = K1 + K2 + L + 1 + K0*d0 + K0 + K0^2;
    m_2a_ext = zeros(1, D);
    s_2a_ext = zeros(1, D);
    p1_2a_ext = zeros(1, D);
    p2_2a_ext = zeros(1, D);
    for i = 1:D
        m_2a_ext(i) = mean(Theta_samp(i, range_conv));
        s_2a_ext(i) = std(Theta_samp(i, range_conv));
        p1_2a_ext(i) = quantile(Theta_samp(i, range_conv), 0.05);
        p2_2a_ext(i) = quantile(Theta_samp(i, range_conv), 0.95);
    end

end % show post results

%% Save file
save([Outputfigsdirname 'Migration_Turkey' filename_extension '_' date]);
