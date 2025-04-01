%% MAIN EMDS 
function main_emds_fully()
    % Simulation parameters
    n = 10;              % Number of nodes
    sim_time = 100;     % Simulation time in (s)
    delta_t = 1;        % Time step in (s)
    sigma_d = 0.3;      % Standard deviation for distance measurements (m) 
    sigma_v = 0.03;     % Standard deviation for velocity measurements (m/s)
    arena_size = 20;    % Size of the arena (m)
    
    [X, V, A] = NodeUtils.initialize_nodes(n, arena_size);
    
    X_hat = zeros(sim_time, n, 2);  % Estimated positions
    X_real = zeros(sim_time, n, 2); % Real positions
    V_real = zeros(sim_time, n, 2); % Real velocities
    errors = zeros(sim_time, 1);    % Position estimation errors
    
    X_real(1, :, :) = X;
    V_real(1, :, :) = V;
    
    % Run MDS to get initial positions
    D_hat = NodeUtils.generate_distance_matrix(X, sigma_d);
    X_hat(1, :, :) = mds(D_hat, 2);
    
    kf = KalmanF.initialize_kalman_filters(n);
    
    for t = 2:sim_time
        [X, V, A] = NodeUtils.update_node_positions(X, V, A, delta_t, arena_size);
        
        X_real(t, :, :) = X;
        V_real(t, :, :) = V;
        
        D_hat_prev = D_hat;
        D_hat = NodeUtils.generate_distance_matrix(X, sigma_d);
        
        V_hat = NodeUtils.generate_velocity_measurements(V, sigma_v);
        
        X_hat_t = emds(D_hat_prev, D_hat, squeeze(X_hat(t-1, :, :)), V_hat, delta_t);
        
        % Apply Kalman Filter for each node
        for i = 1:n
            pos_est = squeeze(X_hat_t(i, :))';
            vel_est = V_hat(i, :);
            
            kf{i} = KalmanF.update_kalman_filter(kf{i}, pos_est, vel_est, delta_t);
            
            X_hat(t, i, :) = kf{i}.x(1:2);
        end
        
        errors(t) = NodeUtils.calculate_position_error(squeeze(X_real(t, :, :)), squeeze(X_hat(t, :, :)));
        
        if mod(t, 10) == 0
            fprintf('Time step: %d/%d, Average error: %.2f meters\n', t, sim_time, mean(errors(1:t)));
        end
    end
    
    VisualizationUtils.visualize_results(X_real, X_hat, errors);
end
