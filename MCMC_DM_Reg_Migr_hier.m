function [Theta_samp, Dec_vec_a, Dec_vec_b] = MCMC_DM_Reg_Migr_hier(Y, U, V, Z, M, Theta_init, prop_params, prior_params)

% [output] = MCMC_DM_Reg_Migr_hier(Y, U, V, Z, M, Theta_init, prop_params, prior_params)
%
% Y is a (T x 1) cell array, each cell contains a matrix of size (N x N)
% Y{t}(i, j): # people migrating from city i to city j.
%
% U is a (T x 1) cell array, each cell matrix of size (N x K)
% U(t)(i, j): the value of the i'th feature of the j'th city
%
% Z is a (T x 1) cell array, each cell contains a matrix of size (N x N x L)
% Z{t}(i, j, l): the value for the l'th feature of cities (i, j)
%
% M: number of iterations
% theta0 is the initial value for the parameter vector, whose size must be
% (K + L x 2)

%% Get the dimensions
T = length(Y);
[N, K1] = size(U{1});
K2 = size(V{1}, 2);
L = size(Z{1}, 3);
D = K1 + K2 + L + 1 + 2*N + 2*2; % extra parameter to scale the alpha's


%% Prior and proposal parameters

% get the prior parameters
mu0_a = prior_params.var0_a;
var0_a = prior_params.var0_a;
alpha0_a = prior_params.alpha0_a;
beta0_a = prior_params.beta0_a;

mu0_b = prior_params.var0_b;
var0_b = prior_params.var0_b;
alpha0_b = prior_params.alpha0_b;
beta0_b = prior_params.beta0_b;

% get the proposal parameters
sigma_prop = prop_params.sigma_prop;
sigma_sc_prop = prop_params.sigma_sc_prop;
sigma0_a_prop = prop_params.sigma0_a_prop;
sigma0_b_prop = prop_params.sigma0_b_prop;

%% Initialisation
theta1 = Theta_init{1};
theta2 = Theta_init{2};
theta3 = Theta_init{3};
theta4 = Theta_init{4};
theta0a = Theta_init{5};
theta0b = Theta_init{6};
thetaha = Theta_init{7};
thetahb = Theta_init{8};

% initialise the log-likelihood
[log_lkl, log_lkl_vec] = calculate_log_lkl_migration(Y, U, V, Z, Theta_init,...
    T, N, L, 2);

% initialise the log-priors
log_prior_vec_a = -0.5*(theta0a - thetaha(1)).^2/thetaha(2)...
    - 0.5*log(2*pi*thetaha(2));

log_prior_vec_b = -0.5*(theta0b - thetahb(1)).^2/thetahb(2)...
    - 0.5*log(2*pi*thetahb(2));

%% Iterations
% Initialise array of theta samples
Theta_samp = zeros(D, M);
Dec_vec_a = zeros(1, M);
Dec_vec_b = zeros(1, M);
j = -1;
for m = 1:M
    if mod(m, 1000) == 0
        disp(m);
    end
    
    %%%%%%%%%%%%% Step 1: Propose factor-related parameters %%%%%%%%%%%%%
    % propose new theta
    j = mod(j+1, 4);
    theta1_prop = theta1 + (j==0)*sigma_prop*randn(K1, 1);
    theta2_prop = theta2 + (j==1)*sigma_prop*randn(K2, 1);
    theta3_prop = theta3 + (j==2)*sigma_prop*randn(L, 1);
    theta4_prop = theta4 + (j==3)*sigma_sc_prop*randn(1, 1);
    
    Theta_prop = {theta1_prop, theta2_prop, theta3_prop, theta4_prop, theta0a, theta0b};
    
    % calculate the log-likelihood
    [log_lkl_prop, log_lkl_vec_prop] = calculate_log_lkl_migration(Y, U, V, Z, ...
        Theta_prop, T, N, L, 2);
        
    log_r = log_lkl_prop - log_lkl;

    if rand < exp(log_r)
        theta1 = theta1_prop;
        theta2 = theta2_prop;
        theta3 = theta3_prop;
        theta4 = theta4_prop;
        log_lkl = log_lkl_prop;
        log_lkl_vec = log_lkl_vec_prop;
    end
    
    %%%%%%%%%%%%% Step 2. update base parameters %%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%% First, the base parameter %%%%%%%%%%%%%%%%%%%%%%%
    theta0a_prop = theta0a + sigma0_a_prop*randn(N, 1);
    
    % calculate the log-likelihood
    Theta_prop = {theta1, theta2, theta3, theta4, theta0a_prop, theta0b};
    [~, log_lkl_vec_prop] = calculate_log_lkl_migration(Y, U, V, Z, Theta_prop,...
        T, N, L, 2);
    
    % calculate the log-prior vector
    log_prior_vec_prop = -0.5*(theta0a_prop - thetaha(1)).^2/thetaha(2)...
        - 0.5*log(2*pi*thetaha(2));
    
    % acceptance ratios
    log_r_vec = log_lkl_vec_prop + log_prior_vec_prop - log_lkl_vec - log_prior_vec_a;
    
    % decision and update
    decision_vec = rand(N, 1) < log_r_vec;
    Dec_vec_a(m) = sum(decision_vec);
    theta0a(decision_vec) = theta0a_prop(decision_vec);
    log_lkl_vec(decision_vec) = log_lkl_vec_prop(decision_vec);
    % log_prior_vec_a(decision_vec) = log_prior_vec_prop(decision_vec);
    
    %%%%%%%%%%%%%%%%%%%%% Secondly, the year effect %%%%%%%%%%%%%%%%%%%%%%%
    theta0b_prop = theta0b + sigma0_b_prop*randn(N, 1);
    
    % calculate the log-likelihood
    Theta_prop = {theta1, theta2, theta3, theta4, theta0a, theta0b_prop};
    [~, log_lkl_vec_prop] = calculate_log_lkl_migration(Y, U, V, Z, Theta_prop,...
        T, N, L, 2);
    
    % calculate the log-prior vector
    log_prior_vec_prop = -0.5*(theta0b_prop - thetahb(1)).^2/thetahb(2)...
        - 0.5*log(2*pi*thetahb(2));
    
    % acceptance ratios
    log_r_vec = log_lkl_vec_prop + log_prior_vec_prop - log_lkl_vec - log_prior_vec_b;
    
    % decision and update
    decision_vec = rand(N, 1) < log_r_vec;
    Dec_vec_b(m) = sum(decision_vec);
    theta0b(decision_vec) = theta0b_prop(decision_vec);
    log_lkl_vec(decision_vec) = log_lkl_vec_prop(decision_vec);
    % log_prior_vec_b(decision_vec) = log_prior_vec_prop(decision_vec);
    
    %%%%%%%%%%%%% Step 3. update the hyperparameters %%%%%%%%%%%%%
    %%% first the base parameter
    % first the mean conditional on the variance
    var_post = 1/(1/var0_a + N/thetaha(2));
    mu_post = var_post*(mu0_a/var0_a + sum(theta0a)/thetaha(2));
    thetaha(1) = sqrt(var_post)*randn + mu_post;
    
    % secondly, the variance conditional on the mean
    alpha_post = alpha0_a + N/2;
    beta_post = beta0_a + sum((theta0a - thetaha(1)).^2)/2;
    thetaha(2) = 1/gamrnd(alpha_post, 1/beta_post);
    
    log_prior_vec_a = -0.5*(theta0a - thetaha(1)).^2/thetaha(2)...
        - 0.5*log(2*pi*thetaha(2));
    
    %%% secondly, the year effect
    % first the mean conditional on the variance
    var_post = 1/(1/var0_b + N/thetahb(2));
    mu_post = var_post*(mu0_b/var0_b + sum(theta0b)/thetahb(2));
    thetahb(1) = sqrt(var_post)*randn + mu_post;
    
    % secondly, the variance conditional on the mean
    alpha_post = alpha0_b + N/2;
    beta_post = beta0_b + sum((theta0b - thetahb(1)).^2)/2;
    thetahb(2) = 1/gamrnd(alpha_post, 1/beta_post);
    
    log_prior_vec_b = -0.5*(theta0b - thetahb(1)).^2/thetahb(2)...
        - 0.5*log(2*pi*thetahb(2));
    
    % store theta's
    Theta_samp(:, m) = [theta1; theta2; theta3; theta4; theta0a; theta0b; thetaha; thetahb];
end