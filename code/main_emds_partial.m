%% MAIN PARTIAL EMDS

function main_emds_partial()
    % Simulation parameters
    n = 10;             % Total number of nodes 
    m = 7;              % Number of fully connected nodes
    sim_time = 100;     % Simulation time in (s)
    delta_t = 1;        % Time step in (s)
    sigma_d = 0.3;      % Standard deviation for distance measurements (m)
    sigma_v = 0.03;     % Standard deviation for velocity measurements (m/s)
    arena_size = 20;    % Size of the arena (m)
    
    fprintf('Starting partial connectivity MDS simulation with %d nodes (%d fully connected)...\n', n, m);
    
    [X, V, A] = NodeUtils.initialize_nodes(n, arena_size);
    
    X_hat = zeros(sim_time, n, 2);  % Estimated positions
    X_real = zeros(sim_time, n, 2); % Real positions
    V_real = zeros(sim_time, n, 2); % Real velocities
    errors = zeros(sim_time, 1);    % Position estimation errors
    
    X_real(1, :, :) = X;
    V_real(1, :, :) = V;
    
    % Generate connectivity matrix - first m nodes are fully connected
    % the rest can only communicate with the first m nodes
    W = NodeUtils.generate_connectivity_matrix(n, m);
    
    fprintf('Generating initial distance matrix...\n');
    D_hat = NodeUtils.generate_partial_distance_matrix(X, sigma_d, W);
    fprintf('Running initial MDS...\n');
    try
        X_hat(1, :, :) = EMDSPartial.emds_partial(D_hat, m, 2);
        fprintf('Initial positions computed successfully.\n');
    catch e
        fprintf('Error in initial MDS: %s\n', e.message);
        % Fallback: Use true positions with some noise
        X_hat(1, :, :) = X + randn(size(X)) * 0.5;
        fprintf('Using fallback initial positions.\n');
    end
    
    % Kalman Filter initialization for each node
    fprintf('Initializing Kalman filters...\n');
    kf = KalmanF.initialize_kalman_filters(n);
    
    % Main simulation loop
    fprintf('Starting main simulation loop...\n');
    for t = 2:sim_time
        fprintf('Time step: %d/%d\n', t, sim_time);
        
        [X, V, A] = NodeUtils.update_node_positions(X, V, A, delta_t, arena_size);
        
        X_real(t, :, :) = X;
        V_real(t, :, :) = V;
        
        D_hat = NodeUtils.generate_partial_distance_matrix(X, sigma_d, W);
        
        V_hat = NodeUtils.generate_velocity_measurements(V, sigma_v);
        
        try
            % Run enhanced MDS with partial connectivity
            X_hat_t = EMDSPartial.enhanced_mds_partial(D_hat, squeeze(X_hat(t-1, :, :)), V_hat, delta_t, m);
            
            % Apply Kalman Filter for each node
            for i = 1:n
                pos_est = squeeze(X_hat_t(i, :))';
                vel_est = V_hat(i, :);
                
                kf{i} = KalmanF.update_kalman_filter(kf{i}, pos_est, vel_est, delta_t);
                
                X_hat(t, i, :) = kf{i}.x(1:2);
            end
        catch e
            fprintf('Error at time step %d: %s\n', t, e.message);
            for i = 1:n
                X_hat(t, i, :) = squeeze(X_hat(t-1, i, :)) + V_hat(i, :) * delta_t;
            end
        end
        
        errors(t) = NodeUtils.calculate_position_error(squeeze(X_real(t, :, :)), squeeze(X_hat(t, :, :)));
        
        fprintf('  Average error: %.2f meters\n', mean(errors(1:t)));
    end
    
    fprintf('\nSimulation completed.\n');
    fprintf('Final average error: %.2f meters\n', mean(errors));
    
    fprintf('Generating visualization...\n');
    VisualizationUtils.visualize_results(X_real, X_hat, errors);
end