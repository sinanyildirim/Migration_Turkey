function [Theta_samp] = MCMC_goc_model_1(Y, U, V, Z, M, Theta_init, sigma_prop)

% [output] = MCMC_goc_model1(Y, U, Z, M, theta_init, sigma_prop)

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
% (1 + 2*K + L)

T = length(Y);
[N, K1] = size(U{1});
K2 = size(V{1}, 2);
L = size(Z{1}, 3);
D = 1 + K1 + K2 + L;

% calculate the likelihood
theta0 = Theta_init{1};
theta1 = Theta_init{2};
theta2 = Theta_init{3};
theta3 = Theta_init{4};

log_lkl = calculate_log_lkl_goc(Y, U, V, Z, Theta_init, T, N, L, 1);

%%

Theta_samp = zeros(D, M);
for m = 1:M
    if mod(m, 1000) == 0
        disp(m);
    end
    
    % propose new theta
    theta0_prop = theta0 + sigma_prop*randn(1, 1);
    theta1_prop = theta1 + sigma_prop*randn(K1, 1);
    theta2_prop = theta2 + sigma_prop*randn(K2, 1);
    theta3_prop = theta3 + sigma_prop*randn(L, 1);
    
    % calculate the log-likelihood
    Theta_prop = {theta0_prop, theta1_prop, theta2_prop, theta3_prop};
    
    [log_lkl_prop] = calculate_log_lkl_goc(Y, U, V, Z, Theta_prop,...
        T, N, L, 1);
        
    % acceptance ratio    
    log_r = log_lkl_prop - log_lkl;
    
    % decision
    if rand < exp(log_r)
        theta0 = theta0_prop;
        theta1 = theta1_prop;
        theta2 = theta2_prop;
        theta3 = theta3_prop;
        log_lkl = log_lkl_prop;
    end
    
    % store the theta parameters
    Theta_samp(:, m) = [theta0; theta1; theta2; theta3];
    
end