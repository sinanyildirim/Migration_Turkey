% This is the main code for the simulated-data experiments for the paper 
% entitled "A model of dynamic migration networks: Explaining Turkey's
% inter-provincial migration flows"
%
% A Dirichlet-multinomial model and a gravity model is fitted to the 
% generated data.
% 
% The user is further guided through the script by comments around the
% code.
%
% Last update: 22 April 2023

clear; clc; close all; fc = 0;
rng(1);

Outputdirname = 'Artificial_outputfiles/'; % This is where the output files are saved
if ~exist(Outputdirname, 'dir')
    mkdir(Outputdirname)
end

%% Generate data
N = 10;
T = 5;
lambda_p = 1000;
phi = 0.95;

writetofile = 1;
saveresults = 1;

model_type_gen = 1; %1 for dirichlet-multinomial, 2 for correlated normal 

K1 = 3; K2 = 3; L = 2;
legends_U = cell(1, K1); for i = 1:K1, legends_U{i} = sprintf('U%d', i); end
legends_V = cell(1, K2); for i = 1:K2, legends_V{i} = sprintf('V%d', i); end
legends_Z = cell(1, L); for i = 1:L, legends_Z{i} = sprintf('Z%d', i); end

% algorithms you want to run
M = 100000; % numbers of iterations
m_burn = round(M/2);
range_conv = m_burn:M;
range_conv_pred = (9*M/10+1):M;
M_burn_pred = length(range_conv_pred);

K0 = 1; % (order+1) of the polynomial for the sender base parameters of the provinces
theta0_common = 1; % sender base parameters (1 for common, 0 for province-based)

%% Run the experiments
theta0_vec = 0:-0.5:-5; L_t0 = length(theta0_vec);
Num_MC_runs = 20;

MSE_Total = repmat({repmat({zeros(3, 2)}, 1, 2)}, L_t0, Num_MC_runs);

for it = 1:L_t0
    theta0 = theta0_vec(it);
    theta_true = {ones(K1, 1), ones(K2, 1), ones(L, 1),  0, theta0};

    for mc = 1:Num_MC_runs
        disp([it, mc]);

        [U, V, Z, Y, P] = generate_migration_data(N, T, theta_true, lambda_p, phi, model_type_gen);

        %% Prepare data for the log-linear (gravity) model
        [X_mod, Y_mod] = convert_data_DM2grav(U, V, Z, Y, K0, theta0_common);
    
        %% Initalize the error metrics for the algorithms
        % Get the dimensions of the arrays of the raw errors
        T_tests_vec = zeros(1, T);

        diag_ind = (0:N-1)*N + (1:N);
        all_ind = 1:N^2;
        non_diag_ind = setdiff(all_ind, diag_ind);
        inds_to_compare_cell = {non_diag_ind, all_ind};
        L_inds = [length(inds_to_compare_cell{1}) length(inds_to_compare_cell{2})];

        % Initialise the arrays of all errors
        Temp_mtx1 = repmat({zeros(L_inds(1), T)}, 3, 1);
        Temp_mtx2 = repmat({zeros(L_inds(2), T)}, 3, 1);
        SE_mtx = repmat({[Temp_mtx1 Temp_mtx2]}, 1, 2);

        % initialize the arrays of true and predicted values
        All_true_vals = cell(1, T);
        All_pred_vals = repmat({cell(1, T)}, 1, 2);

        %% Run the methods
        for t = 1:T
            disp([it, mc, t]);
            % get the data (mostly indices) for the training and test years
            train_sub_t = setdiff(1:T, t);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%% Method 1: MCMC %%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Train with MCMC
            Y_train = Y(train_sub_t); U_train = U(train_sub_t);
            V_train = V(train_sub_t); Z_train = Z(train_sub_t);

            %%% Predict with MCMC
            U_test = U{t}; V_test = V{t}; Z_test = Z{t};
            Y_test = Y{t}; P_test = sum(Y_test, 2);
            log_Y_test = log(max(Y_test, 1));
            Prop_test = Y_test./sum(Y_test, 2);

            [Theta_samp, range_theta] = MCMC_Migration(Y_train, U_train, V_train, ...
                Z_train, train_sub_t, K0, theta0_common, M);

            % perform prediction
            Theta_for_pred = Theta_samp(:, range_conv_pred);
            [Y_pred_MCMC, log_Y_pred_MCMC, Prop_pred_MCMC] = ...
                pred_migration({P_test}, {U_test}, {V_test}, {Z_test}, t, Theta_for_pred, ...
                range_theta);

            %%%%%%%%%%%%%%%%%%%%% Method 2: Gravity model %%%%%%%%%%%%%%%%%%%%%%
            % determine the indices of test and training data
            rows_per_year = N*(N-1);
            n = rows_per_year*T;
            test_inds = ((t-1)*rows_per_year+1):t*rows_per_year;
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
            Y_pred_grav = zeros(N);
            log_Y_pred_grav = zeros(N);

            log_Y_pred_grav(inds_to_compare_cell{1}) = log_y_pred;
            Y_pred_grav(inds_to_compare_cell{1}) = y_pred;

            % take the transposes to put in the right form
            log_Y_pred_grav = log_Y_pred_grav';
            Y_pred_grav = Y_pred_grav';

            % fill in the diagonal terms
            y_pred_diag_temp = sum(Y_test, 2) - sum(Y_pred_grav, 2);
            Y_pred_grav(diag_ind) = y_pred_diag_temp;
            log_Y_pred_grav(diag_ind) = log(max(1, y_pred_diag_temp));

            % proportions
            Prop_pred_grav = Y_pred_grav./sum(Y_test, 2);

            %%%%%%%%% Store the predictions and the true values %%%%%%%%%
            All_pred_vals{1}{t} = {Y_pred_MCMC{1}, log_Y_pred_MCMC{1}, Prop_pred_MCMC{1}};
            All_pred_vals{2}{t} = {Y_pred_grav, log_Y_pred_grav, Prop_pred_grav};
            All_true_vals{t} = {Y_test, log_Y_test, Prop_test};

            %%%%%%%%%%%%%%%%%%%%%%%%% Error metrics %%%%%%%%%%%%%%%%%%%%%%%%%
            for i = 1:3 % type of prediction (log, linear, proportion)
                for j = 1:2 % non-diagonal and all
                    inds_to_compare = inds_to_compare_cell{j};
                    for algo_type = [1 2] % MCMC and gravity

                        % linear errors
                        % non-diagonal errors
                        pred_vals = All_pred_vals{algo_type}{t}{i}(inds_to_compare);
                        true_vals = All_true_vals{t}{i}(inds_to_compare);

                        temp_err = pred_vals - true_vals;

                        Temp_vec_SE = temp_err.^2;

                        % These are the arrays of raw errrors and they necessary for the median
                        SE_mtx{algo_type}{i, j}(:, t) = Temp_vec_SE;
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
                    MSE_total{algo_type}(i, j) = mean(SE_mtx{algo_type}{i, j}(:));
                end
            end
        end
        if saveresults == 1
            save([Outputdirname 'Migration_artificial_data' ...
                sprintf('genmodel_%d_N_%d_T_%d_K1_%d_K2_%d_L_%d_lambda_pop_%d_theta0_minus%d_MCrun_%d', ...
                model_type_gen, N, T, K1, K2, L, lambda_p, -10*theta0, mc)]);
        end
        MSE_Total{it, mc} = MSE_total;
    end
end

%% plot results
theta_0_vec = 0:-0.5:-5; L_t0 = length(theta_0_vec);
Num_MC_runs = 10; N = 10; T = 5; K1 = 3; K2 = 3; L = 2; lambda_p = 1000;
for itit = 1:L_t0
    theta_0 = theta_0_vec(itit);
    for mcmc = 1:Num_MC_runs
        load([Outputdirname 'Migration_artificial_data' ...
                sprintf('genmodel_%d_N_%d_T_%d_K1_%d_K2_%d_L_%d_lambda_pop_%d_theta0_minus%d_MCrun_%d', ...
                model_type_gen, N, T, K1, K2, L, lambda_p, -10*theta_0, mcmc)]);

        MSE_lin{1}(itit, mcmc, 1) = MSE_total{1}(1, 1);
        MSE_lin{2}(itit, mcmc, 1) = MSE_total{2}(1, 1);

        MSE_lin{1}(itit, mcmc, 2) = MSE_total{1}(1, 2);
        MSE_lin{2}(itit, mcmc, 2) = MSE_total{2}(1, 2);
    end
end

%%
fc = fc + 1; figure(fc);
subplot(1, 2, 1);
plot(theta_0_vec, mean(MSE_lin{1}(:, :, 1)./MSE_lin{2}(:,:, 1), 2), '.-');
xlabel('$\theta_{0}$', 'Interpreter', 'Latex');
ylabel('ratio of MSEs (DM/grav)');
set(gca, 'xtick', fliplr(theta_0_vec));
title('Diagonals excluded');
subplot(1, 2, 2);
plot(theta_0_vec, mean(MSE_lin{1}(:, :, 2)./MSE_lin{2}(:,:, 2), 2), '.-');
xlabel('$\theta_{0}$', 'Interpreter', 'Latex');
ylabel('ratio of MSEs (DM/grav)');
set(gca, 'xtick', fliplr(theta_0_vec));
title('Diagonals included');