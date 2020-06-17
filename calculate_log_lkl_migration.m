function [log_lkl, log_lkl_vec] = calculate_log_lkl_migration(Y, U, V, Z, Thetas,...
    T, N, L, model_type)

if model_type == 1 % This is the Multinomial model
    
    theta1 = Thetas{1};
    theta2 = Thetas{2};
    theta3 = Thetas{3};
    theta0a = Thetas{4};
    theta0b = Thetas{4};
    
    log_lkl = 0;
    log_lkl_vec = zeros(N, 1);
    
    for t = 1:T
        u = U{t}*theta1;
        v = V{t}*theta2;
        
        Z_t_2D = reshape(Z{t}, N*N, L);
        
        log_prob_mtx = theta0a + (t-1)*theta0b + u*ones(1, N) + ones(N, 1)*v' + reshape(Z_t_2D*theta3, N, N);
        log_prob_mtx((0:N-1)*N + (1:N)) = 0;
        
        % normalise the log-probabilities
        log_prob_mtx = log_prob_mtx - repmat(max(log_prob_mtx, [], 2), 1, N);
        
        log_p_max = max(log_prob_mtx, [], 2);
        R = repmat(log_p_max, 1, N);
        log_sum_vec = log(sum(exp(log_prob_mtx - R), 2)) + log_p_max;
        log_prob_mtx_norm = log_prob_mtx - repmat(log_sum_vec, 1, N);
        
        log_lkl_inc = Y{t}.*log_prob_mtx_norm;
        log_lkl = log_lkl + sum(log_lkl_inc(:));
        log_lkl_vec = log_lkl_vec + sum(log_lkl_inc, 2);
    end
    
    
elseif model_type == 2 % This covers both types of Dirichlet-Multinomial model
    
    theta1 = Thetas{1};
    theta2 = Thetas{2};
    theta3 = Thetas{3};
    theta4 = Thetas{4};
    theta0a = Thetas{5};
    theta0b = Thetas{6};
        
    % if theta0 has dimension N, it has to be replicated
    theta0a = repmat(theta0a, 1, length(theta0a));
    theta0b = repmat(theta0b, 1, length(theta0b));
    
    log_lkl = 0;
    log_lkl_vec = zeros(N, 1);
    
    for t = 1:T
        u = U{t}*theta1;
        v = V{t}*theta2;
        
        Z_t_2D = reshape(Z{t}, N*N, L);
        
        log_alpha_mtx = theta4 + theta0a + theta0b*(t-1) + u*ones(1, N) + ones(N, 1)*v' ...
            + reshape(Z_t_2D*theta3, N, N);
        
        log_alpha_mtx((0:N-1)*N + (1:N)) = theta4;
        alpha_mtx = exp(log_alpha_mtx);
        
        log_lkl_inc = gammaln(sum(alpha_mtx, 2)) - gammaln(sum(Y{t} + alpha_mtx, 2))...
            + sum(gammaln(Y{t}+alpha_mtx), 2) - sum(gammaln(alpha_mtx), 2);
        
        log_lkl = log_lkl + sum(log_lkl_inc);
        log_lkl_vec = log_lkl_vec + log_lkl_inc;
    end
end