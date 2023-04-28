function [X_g, Y_g] = convert_data_DM2grav(U, V, Z, Y, K0, theta0_common)

% [X_g, Y_g] = convert_data_DM2grav(U, V, Z, Y, K0, theta0_common)
% 
% This function converts the data suitable to the DM model to the classical
% table format that is suitable to the gravity model that assumes that the
% the log-migration counts regress on the factors.
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
% K0 is the (polynomial order + 1) for the baseline probability parameter
% theta0_common: set to 1 for a common baseline probability parameter, set to 0
% for a distinct baseline probability parameter per province.
% 
% Sinan Yildirim
% 22.04.2023

T = length(Y);
[N, K1] = size(U{1});
K2 = size(V{1}, 2);
L = size(Z{1}, 3);
d0 = (theta0_common == 0)*N + (theta0_common == 1)*1;

rows_per_year = N*(N - 1);
n = T*rows_per_year;
Y_g = zeros(n, 1);

X_g = zeros(n, max(K0*d0, 1) + K1+K2+L+1);

r = 0;
for t = 1:T
    P = sum(Y{t}, 2);
    for i = 1:N
        for j = 1:N
            if i ~= j
                r = r + 1;
                Y_g(r) = log(max(1, Y{t}(i, j)));

                % prepare the one-hot parts of the X row:
                if d0 == N
                    sender_part = [zeros(i-1, 1); 1; zeros(N-i, 1)]*(t-1).^(0:K0-1);
                    sender_part = sender_part(:)';
                else
                    sender_part = 1;
                end

                % if no sender/receiver base parameter, add an intercept
                intercept_part = ones(1, K0 == 0);

                X_g(r, :) = [intercept_part sender_part ...
                    U{t}(i, :) V{t}(j, :) squeeze(Z{t}(i, j, :))' log(max(1, P(i)))];
            end
        end
    end
end