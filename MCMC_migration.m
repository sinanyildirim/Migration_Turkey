function [Theta_samp, range_theta] = MCMC_migration(Y, U, V, Z, year_vec, K0, theta0_common, M, theta_init)

% [Theta_samp, range_theta] = MCMC_migration(Y, U, V, Z, year_vec, K0, theta0_common, M, theta_init)
%
% Input data contains the migration counts and the one-way and two-way
% external factors regarding N provinces for T years.
% 
% Y is a (T x 1) cell array, each cell contains a matrix of size (N x N)
% Y{t}(i, j): # people migrating from city i to city j.
%
% U is a (T x 1) cell array, each cell matrix of size (K1 x N)
% U(t)(k, i): the value of the i'th feature of the i'th migration-sending
% province
%
% V is a (T x 1) cell array, each cell matrix of size (K2 x N)
% V(t)(k, j): the value of the k'th feature of the i'th migration-receiving
% province.
%
% Z is a (T x 1) cell array, each cell contains a matrix of size (N x N x L)
% Z{t}(i, j, l): the value for the l'th feature of cities (i, j)
%
% K0: polynomial order of the base parameter
% theta0_common: set 1 for a common base parameter, 0 for province-based base parameter
% M: number of iterations
% 
% theta_init: initial theta value. If not specified, it is set to default 
% value inside the function
% 
% Those two moments are also random.
% 
% year_vec: This is the vector of years that correspond to time steps 1:T
%
% Output is the a M x D matrix of samples from the MCMC. Each column 
% (sample) is formed as 
% 
% [theta1; theta2; theta3; theta4; theta0(:), theta0_mu, theta0_Cov(:)]
% 
% Sinan Yildirim, 22 April 2023

%% Get the dimensions
T = length(Y);
[N, K1] = size(U{1});
K2 = size(V{1}, 2);
L = size(Z{1}, 3);

d0 = (theta0_common == 0)*N + (theta0_common == 1)*1;

% dimension of theta
D = K1 + K2 + L + 1 + K0*d0 + K0 + K0^2;


%% the ranges %%%%
range_theta{1} = 1:K1;
range_theta{2} = range_theta{1}(end)+(1:K2);
range_theta{3} = range_theta{2}(end)+(1:L);
range_theta{4} = range_theta{3}(end) + 1;
range_theta{5} = range_theta{4}(end) +(1:(K0*d0));
if K0 > 0, range_last = range_theta{5}; end
range_theta{6} = range_last(end)+(1:K0);
if K0 > 0, range_last = range_theta{6}; end
range_theta{7} = range_last(end)+(1:K0^2);

%% Proposal parameters
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

%%% Proposal parameters
S1 = eye(K1)/chol(cov(UU));
S2 = eye(K2)/chol(cov(VV));
S3 = eye(L)/chol(cov(ZZ));
sigma_prop = {2/sqrt(N^2*T*K1)*S1, 2/sqrt(N^2*T*K2)*S2, ...
    2/sqrt(N^2*T*L)*S3, 0.02, 0.1*[1; 0.1*ones(K0-1, 1)]};

%% Prior parameters
% prior covariances for the factor related parameters
prior_var1 = 10^4*ones(K1, 1);
prior_var2 = 10^4*ones(K2, 1);
prior_var3 = 10^4*ones(L, 1);
prior_var4 = 10^4;

%%% Prior hyperparameters
mu0_mean = zeros(K0, 1);
mu0_Cov = 100*eye(K0);
nu0 = K0;
Psi0 = eye(K0);

%% Initial values
if nargin < 9 % if an initial value is not provided
    theta1 = zeros(K1, 1);
    theta2 = zeros(K2, 1);
    theta3 = zeros(L, 1);
    theta4 = 0;
    theta0 = [-10*ones(min(1, K0), d0); zeros(K0-1, d0)];
    % These variables become hyperparameters if fe_re = 0
    theta0_mu = [0; zeros(K0-1, 1)];
    theta0_Cov = diag([1 0.01*ones(1, K0-1)]);
else % if an initial value is provided
    theta1 = theta_init(range_theta{1});
    theta2 = theta_init(range_theta{2});
    theta3 = theta_init(range_theta{3});
    theta4 = theta_init(range_theta{4});
    theta0 = reshape(theta_init(range_theta{5}), K0, d0);
    theta0_mu = theta_init(range_theta{6});
    theta0_Cov = reshape(theta_init(range_theta{7}), K0, K0);
end

%% Initialisation of the probability densities
% initialise the log-likelihood
theta_init = {theta1, theta2, theta3, theta4, theta0};
[~, log_lkl_vec] = log_lkl_migration(Y, U, V, Z, theta_init, year_vec);

% Calculate the prior probability vector of theta0 (must be 1 x d0)
temp_mtx = chol(theta0_Cov);
log_prior_0_vec = -0.5*sum((temp_mtx'\(theta0 - theta0_mu)).^2, 1)...
    - 0.5*K0*log(2*pi) - sum(log(diag(temp_mtx)));

%% Iterations
% Initialise array of theta samples
Theta_samp = zeros(D, M);
for m = 1:M
    if mod(m, 1000) == 0
        fprintf('%d ', m);
        if mod(m, 10000) == 0
            fprintf('\n');
        end
    end

    %%%%%%%%%%%%% Step 1: Update factor-related parameters %%%%%%%%%%%%%
    % propose new theta    
    j = mod(m-1, 3) + 1;
    theta1_prop = theta1 + (j==1)*sigma_prop{1}*randn(K1, 1);
    theta2_prop = theta2 + (j==2)*sigma_prop{2}*randn(K2, 1);
    theta3_prop = theta3 + (j==3)*sigma_prop{3}*randn(L, 1);

    % combine the proposed variables
    theta_prop = {theta1_prop, theta2_prop, theta3_prop, theta4, theta0};

    % ration of the priors
    log_prior_ratio = ...
        sum(-(theta1_prop.^2 - theta1.^2)./(2*prior_var1))...
        + sum(-(theta2_prop.^2 - theta2.^2)./(2*prior_var2))...
        + sum(-(theta3_prop.^2 - theta3.^2)./(2*prior_var3));
    
    % calculate the log-likelihood of the proposal
    [~, log_lkl_vec_prop] = log_lkl_migration(Y, U, V, Z, theta_prop, year_vec);
    
    % calculate the log acceptance ratio
    log_r = sum(log_lkl_vec_prop) - sum(log_lkl_vec) + log_prior_ratio;

    % decision
    decision = rand < exp(log_r);
    if decision == 1
        theta1 = theta1_prop;
        theta2 = theta2_prop;
        theta3 = theta3_prop;
        log_lkl_vec = log_lkl_vec_prop;
    end

    %%%%%%%%%%%%% Step 2a. update sender base parameters (polynomial coefficients) %%%%%
    theta0_prop = theta0 + sigma_prop{5}*randn(K0, d0);

    % calculate the log-likelihood of the proposal
    theta_prop = {theta1, theta2, theta3, theta4, theta0_prop};
    [~, log_lkl_vec_prop] = log_lkl_migration(Y, U, V, Z, theta_prop, year_vec);

    % calculate the log-prior (vector) - this is a 1 x d_b vector
    temp_mtx = chol(theta0_Cov);
    log_prior_0_vec_prop = -0.5*sum((temp_mtx'\(theta0_prop - theta0_mu)).^2, 1)...
        - 0.5*K0*log(2*pi) - sum(log(diag(temp_mtx)));

    if d0 == 1
        % acceptance ratio
        log_r = sum(log_lkl_vec_prop) - sum(log_lkl_vec) ...
            + log_prior_0_vec_prop - log_prior_0_vec;
        
        % decision and update
        decision = rand < exp(log_r);
        if decision == 1
            theta0 = theta0_prop;
            log_lkl_vec = log_lkl_vec_prop;
            log_prior_0_vec = log_prior_0_vec_prop;
        end
    elseif d0 == N
        % acceptance ratios
        log_r_vec = log_lkl_vec_prop' + log_prior_0_vec_prop ...
            - log_lkl_vec' - log_prior_0_vec;
        
        % decision and update
        decision_vec = rand(1, d0) < exp(log_r_vec);
        theta0(:, decision_vec) = theta0_prop(:, decision_vec);
        log_lkl_vec(decision_vec) = log_lkl_vec_prop(decision_vec);
        log_prior_0_vec(decision_vec) = log_prior_0_vec_prop(decision_vec);
    end

    %%%%%%%%%%%%% Step 2b. Update the base hyperparameters %%%%%%%%%%%%%
    if d0 == N
        % first, sample the mean
        Cov_post = mu0_Cov - mu0_Cov*((theta0_Cov +  mu0_Cov*N)\(N*mu0_Cov));
        mu_post = Cov_post*(mu0_Cov\mu0_mean + theta0_Cov\(sum(theta0, 2)));
        theta0_mu = mvnrnd(mu_post, Cov_post)';
    
        % secondly, the variance conditional on the mean
        nu_post = nu0 + N;
        Psi_post = Psi0 + (theta0 - theta0_mu)*(theta0 - theta0_mu)';
        theta0_Cov = iwishrnd(Psi_post, nu_post);  
    
        % refresh the log-prior vector
        temp_mtx = chol(theta0_Cov);
        log_prior_0_vec = -0.5*sum((temp_mtx'\(theta0 - theta0_mu)).^2, 1)...
            - 0.5*K0*log(2*pi) - sum(log(diag(temp_mtx)));
    end
    
    %%%%%%%%%%%%% Step 3: Update the scale parameter %%%%%%%%%%%%%
    theta4_prop = theta4 + sigma_prop{4}*randn;
    theta_prop = {theta1, theta2, theta3, theta4_prop, theta0};
    [~, log_lkl_vec_prop] = log_lkl_migration(Y, U, V, Z, theta_prop,...
        year_vec);

    log_prior_ratio = -(theta4_prop.^2 - theta4.^2)./(2*prior_var4);

    % acceptance ratios
    log_r = sum(log_lkl_vec_prop) - sum(log_lkl_vec) + log_prior_ratio;
    decision = rand < exp(log_r);

    if decision == 1
        theta4 = theta4_prop;
        log_lkl_vec = log_lkl_vec_prop;
    end

    % store theta's
    Theta_samp(:, m) = [theta1; theta2; theta3; theta4; theta0(:); ...
        theta0_mu; theta0_Cov(:)];
end


