function [Theta_samp] = MCMC_DM_Reg_Migr(Y, U, V, Z, M, Theta_init, prop_params)

% [output] = MCMC_DM_Reg_Migr(Y, U, V, Z, M, Theta_init, prop_params)
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

T = length(Y);
[N, K1] = size(U{1});
K2 = size(V{1}, 2);
L = size(Z{1}, 3);
D = K1 + K2 + L + 3;

% calculate the likelihood
theta1 = Theta_init{1};
theta2 = Theta_init{2};
theta3 = Theta_init{3};
theta4 = Theta_init{4};
theta0a = Theta_init{5};
theta0b = Theta_init{6};

% get the proposal parameters
sigma_prop = prop_params.sigma_prop;
sigma_sc_prop = prop_params.sigma_sc_prop;
sigma0_a_prop = prop_params.sigma0_a_prop;
sigma0_b_prop = prop_params.sigma0_b_prop;

% initialise the log-likelihood
[log_lkl, ~] = calculate_log_lkl_migration(Y, U, V, Z, Theta_init, T, N, L, 2);

%%
j = -1;
Theta_samp = zeros(D, M);
for m = 1:M
    if mod(m, 1000) == 0
        disp(m);
    end
    
    % propose new theta
    j = mod(j+1, 5);
    
    theta1_prop = theta1 + (j==0)*sigma_prop*randn(K1, 1);
    theta2_prop = theta2 + (j==1)*sigma_prop*randn(K2, 1);
    theta3_prop = theta3 + (j==2)*sigma_prop*randn(L, 1);
    theta4_prop = theta4 + (j==3)*sigma_sc_prop*randn(1, 1);
    theta0a_prop = theta0a + (j==4)*sigma0_a_prop*randn(1, 1);
    theta0b_prop = theta0b + (j==4)*sigma0_b_prop*randn(1, 1);    
       
    % calculate the log-likelihood
    Theta_prop = {theta1_prop, theta2_prop, theta3_prop, theta4_prop, theta0a_prop, theta0b_prop};    
    
    [log_lkl_prop, ~] = calculate_log_lkl_migration(Y, U, V, Z, Theta_prop, T, N, L, 2);

    % acceptance ratio
    log_r = log_lkl_prop - log_lkl;

    % decision
    if rand < exp(log_r)
        theta1 = theta1_prop;
        theta2 = theta2_prop;
        theta3 = theta3_prop;
        theta4 = theta4_prop;
        theta0a = theta0a_prop;
        theta0b = theta0b_prop;
        
        log_lkl = log_lkl_prop;
    end
    
    % store the theta parameters
    Theta_samp(:, m) = [theta1; theta2; theta3; theta4; theta0a; theta0b];
    
end