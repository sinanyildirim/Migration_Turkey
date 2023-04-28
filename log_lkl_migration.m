function [log_lkl, log_lkl_vec] = log_lkl_migration(Y, U, V, Z, theta, year_vec)

% [log_lkl, log_lkl_vec, Y_pred, p_pred] = log_lkl_migration(Y, ...
%    U, V, Z, Thetas, year_vec, year_init)
% 
% This function calculates the log-likelihood of the observations according
% to the Dirichlet-multinomial model for migration.
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
% theta is a cell of parameters, in the form 
% 
% {theta1, theta2, theta3, theta4, theta0}
% 
% where 
% - theta0 is a K0 x N matrix, whose i'th column is the baseline
% probability parameter of the i'th origin
% Each column of theta0 is assumed to have a normal distribution with mean
% 
% year_vec: This is the vector of years that correspond to time steps 1:T


% number of years
T = length(year_vec);

[N, L] = size(Z{1}, [2, 3]);

theta1 = theta{1}; % params for sender factors
theta2 = theta{2}; % params for receiver factors
theta3 = theta{3}; % params for pair factors
theta4 = theta{4}; % scale parameter
theta0 = theta{5}; % receiver base, may be K0r x 1 or K0 x N
K0 = size(theta0, 1);

% initialize
log_lkl_vec = zeros(N, 1);

% indices on the diagonal
diag_inds = (0:N-1)'*N + (1:N)';

for t = 1:T
    u = U{t}*theta1;
    v = V{t}*theta2;

    % get the polynomial bases
    t_vec_s = (year_vec(t) - 1).^(0:(K0-1));

    Z_t_2D = reshape(Z{t}, N*N, L);
    
    % prepare the parameters of the model
    log_prob_mtx = theta4 + (t_vec_s*theta0)' + u + v' + reshape(Z_t_2D*theta3, N, N);
    log_prob_mtx(diag_inds) = theta4;

    % calculate the alpha parameters
    alpha_mtx = exp(log_prob_mtx);
        
    % calculate the log-likelihood
    y = Y{t};
    log_lkl_inc = gammaln(sum(alpha_mtx, 2)) - gammaln(sum(y + alpha_mtx, 2))...
        + sum(gammaln(y + alpha_mtx), 2) - sum(gammaln(alpha_mtx), 2);

    % update the log-lkl vector
    log_lkl_vec = log_lkl_vec + log_lkl_inc;
end

log_lkl = sum(log_lkl_vec);