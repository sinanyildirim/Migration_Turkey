function [U, V, Z, Y, P] = generate_migration_data(N, T, theta, lambda_p, phi, model_type)

% [U, V, Z, Y, P] = generate_migration_data(N, T, theta, lambda_p, phi, model_type)
% 
% This function generates example migration data for T time periods and for
% N provinces. 
% 
% Sinan Yildirim, 22.04.2023

theta1 = theta{1}; 
theta2 = theta{2};
theta3 = theta{3};
theta4 = theta{4};
theta0 = theta{5};

K1 = length(theta1);
K2 = length(theta2);
L = length(theta3);

sigma_x = sqrt(1 - phi^2);

U = cell(1, T);
V = cell(1, T);
Z = cell(1, T);
Y = cell(1, T);
P = cell(1, T);

diag_inds = (0:N-1)*N + (1:N);

for t = 1:T
    if t == 1
        U{t} = randn(N, K1);
        V{t} = randn(N, K2);
        Z{t} = randn(N, N, L);
        P{t} = poissrnd(lambda_p, N, 1);
    else
       U{t} = phi*U{t-1} + sigma_x*randn(N, K1);
       V{t} = phi*V{t-1} + sigma_x*randn(N, K2);
       Z{t} = phi*Z{t-1} + sigma_x*randn(N, N, L);
       P{t} = sum(Y{t-1}, 1)';
    end
    
    u = U{t}*theta1;
    v = V{t}*theta2;
    Z_t_2D = reshape(Z{t}, N*N, L);
    
    alpha_mtx = theta0 + theta4 + u + v' + reshape(Z_t_2D*theta3, N, N);
    alpha_mtx(diag_inds) = theta4;
    
    for n = 1:N 
        if model_type == 1
            prob_vec = dirichlet_rnd(exp(alpha_mtx(n, :)), 1);
            Y{t}(n, :) = mnrnd(P{t}(n), prob_vec);
        elseif model_type == 2
            % calculate the covariance matrix
            a_v = exp(alpha_mtx(n, :));
            p_v = a_v/sum(a_v);
            a0 = sum(a_v);

            D = -(p_v'*p_v)*P{t}(n)*(P{t}(n) + a0)/(1 + a0);
            D(diag_inds) = p_v.*(1 - p_v)*P{t}(n)*(P{t}(n) + a0)/(1 + a0);
            D = (D + D')/2;
            D = D + eye(N);
            mu_vec = P{t}(n)*p_v;
            % randomize the covariance matrix
            Cov_mtx = wishrnd(D/N, N);
            
            % sample the migration counts
            Y{t}(n, :) = mvnrnd(mu_vec, Cov_mtx);
            Y{t}(n, :) = round(max(0, Y{t}(n, :)));

            % renew the population
            P{t}(n) = sum(Y{t}(n, :));
        end
    end
end