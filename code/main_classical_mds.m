%% MAIN CLASSICAL MDS SIMULATION

function main_classical_mds()
    % Simulation parameters
    n = 6;              % Number of nodes
    sim_time = 100;     % Simulation time in (s)
    delta_t = 1;        % Time step in (s)
    sigma_d = 0.3;      % Standard deviation for distance measurements (meters)
    arena_size = 20;    % Size of the arena (meters)
    
    fprintf('Starting Classical MDS simulation with %d nodes...\n', n);
    
    [X, V, A] = NodeUtils.initialize_nodes(n, arena_size);
    
    X_hat = zeros(sim_time, n, 2);  % Estimated positions
    X_real = zeros(sim_time, n, 2); % Real positions
    V_real = zeros(sim_time, n, 2); % Real velocities
    errors = zeros(sim_time, 1);    % Position estimation errors
    
    X_real(1, :, :) = X;
    V_real(1, :, :) = V;
    
    D_hat = NodeUtils.generate_distance_matrix(X, sigma_d);
    X_hat(1, :, :) = mds(D_hat, 2);
    
    kf = KalmanF.initialize_kalman_filters(n);
    
    % Main simulation loop
    fprintf('Starting main simulation loop...\n');
    for t = 2:sim_time
        [X, V, A] = NodeUtils.update_node_positions(X, V, A, delta_t, arena_size);

        X_real(t, :, :) = X;
        V_real(t, :, :) = V;
        
        D_hat = NodeUtils.generate_distance_matrix(X, sigma_d);
        
        % Run classical MDS 
        X_mds = mds(D_hat, 2);
        
        if t > 1
            V_hat = (squeeze(X_mds) - squeeze(X_hat(t-1, :, :))) / delta_t;
        else
            V_hat = zeros(n, 2);
        end
        
        for i = 1:n
            % Get position estimates from MDS
            pos_est = squeeze(X_mds(i, :))';
            vel_est = V_hat(i, :);
            
            % Update Kalman Filter
            kf{i} = KalmanF.update_kalman_filter(kf{i}, pos_est, vel_est, delta_t);
            
            % Store the filtered position estimates
            X_hat(t, i, :) = kf{i}.x(1:2);
        end
        
        errors(t) = NodeUtils.calculate_position_error(squeeze(X_real(t, :, :)), squeeze(X_hat(t, :, :)));
        
        if mod(t, 10) == 0
            fprintf('Time step: %d/%d, Average error: %.2f meters\n', t, sim_time, mean(errors(1:t)));
        end
    end
    
    fprintf('\nSimulation completed.\n');
    fprintf('Final average error: %.2f meters\n', mean(errors));
    
    fprintf('Generating visualization...\n');
    VisualizationUtils.visualize_results(X_real, X_hat, errors);
end