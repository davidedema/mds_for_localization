function X_hat = emds(D_prev, D_curr, X_prev, V_hat, delta_t)
    n = size(D_curr, 1);
    dim = size(X_prev, 2);
    
    % Initialize optimization variables (2n variables)
    % First n variables are positions at t-1, last n are positions at t
    X0 = zeros(2*n, dim);
    X0(1:n, :) = X_prev;  
    X0(n+1:end, :) = X_prev + V_hat * delta_t;  
    
    stress_fun = @(x) emds_stress_function(x, D_prev, D_curr, n, dim);
    
    % Define constraints (positions at t are related to positions at t-1 by velocity)
    % A*x <= b + epsilon_v and -A*x <= -b + epsilon_v
    A = zeros(n*dim, 2*n*dim);
    b = zeros(n*dim, 1);
    
    for i = 1:n
        for d = 1:dim
            row_idx = (i-1)*dim + d;
            col_idx1 = (i-1)*dim + d;
            col_idx2 = n*dim + (i-1)*dim + d;
            
            A(row_idx, col_idx1) = 1;
            A(row_idx, col_idx2) = -1;
            b(row_idx) = V_hat(i, d) * delta_t;
        end
    end
    % Combined constraints
    epsilon_v = 0.05;
    A_combined = [A; -A];
    b_combined = [b + epsilon_v; -b + epsilon_v];
    
    options = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'interior-point', ...
        'MaxIterations', 2000, 'MaxFunctionEvaluations', 10000, ...
        'OptimalityTolerance', 1e-6, 'StepTolerance', 1e-6);
    
    % Optimization
    x = fmincon(stress_fun, reshape(X0', [], 1), A_combined, b_combined, [], [], [], [], [], options);
    
    X_reshaped = reshape(x, dim, 2*n)';
    
    X_hat = X_reshaped(n+1:end, :);
end

% eMDS stress function
function stress = emds_stress_function(x, D_prev, D_curr, n, dim)
    X_reshaped = reshape(x, dim, 2*n)';
    X_prev = X_reshaped(1:n, :);
    X_curr = X_reshaped(n+1:end, :);
    
    % Calculate stress at t-1
    stress_prev = 0;
    for i = 1:n
        for j = (i+1):n
            d_ij = norm(X_prev(i, :) - X_prev(j, :));
            stress_prev = stress_prev + (d_ij - D_prev(i, j))^2;
        end
    end
    
    % Calculate stress at t
    stress_curr = 0;
    for i = 1:n
        for j = (i+1):n
            d_ij = norm(X_curr(i, :) - X_curr(j, :));
            stress_curr = stress_curr + (d_ij - D_curr(i, j))^2;
        end
    end
    
    stress = stress_prev + stress_curr;
end