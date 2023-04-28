function [Y_pred, log_Y_pred, Prop_pred] = pred_migration(P, U, V, Z, year_vec, Thetas, range_theta)

% [Y_pred, log_Y_pred, Prop_pred] = pred_migration(P, U, V, Z, year_vec, Thetas, range_theta)
% 
% This function predicts migrations by calculating their expected value
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
% Thetas is a matrix of theta values, in the form 
% 
% [theta1; theta2; theta3; theta4; theta0(:), theta0_mu, theta0_Cov(:)]
% 
% where 
% - theta0 is a K0 x N matrix, whose i'th column is the baseline
% probability parameter of the i'th origin
% Each column of theta0 is assumed to have a normal distribution with mean
% 
% year_vec: This is the vector of years that correspond to time steps 1:T
% 
% Sinan Yildirim
% Last update: 22 April 2023

% number of years
T = length(year_vec);

[N, L] = size(Z{1}, [2, 3]);
% indices on the diagonal
diag_inds = (0:N-1)'*N + (1:N)';

Prop_pred = repmat({zeros(N)}, 1, T);
Y_pred = repmat({zeros(N)}, 1, T);
log_Y_pred= repmat({zeros(N)}, 1, T);

K0 = length(range_theta{6});
d0 = length(range_theta{5})/K0;

M = size(Thetas, 2);

for m = 1:M
           
    theta_curr = Thetas(:, m);
    theta1 = theta_curr(range_theta{1}); % params for sender factors
    theta2 = theta_curr(range_theta{2}); % params for receiver factors
    theta3 = theta_curr(range_theta{3}); % params for pair factors
    theta4 = theta_curr(range_theta{4}); % scale parameter
    theta0 = reshape(theta_curr(range_theta{5}), K0, d0); % receiver base, may be K0r x 1 or K0 x N

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
        p_pred = alpha_mtx./sum(alpha_mtx, 2);
        y_pred = P{t}.*p_pred;

        % update the predictions
        Y_pred{t} = Y_pred{t} + y_pred/M;
        log_Y_pred{t} = log_Y_pred{t} + log(y_pred)/M;
        Prop_pred{t} = Prop_pred{t} + p_pred/M;
    end
end