classdef NodeUtils
    methods (Static)
        % Initialize node positions, velocities, and accelerations
        function [X, V, A] = initialize_nodes(n, arena_size)
            % Initialize positions randomly within the arena
            X = (rand(n, 2) - 0.5) * 2 * arena_size;
            
            % Initialize velocities randomly
            speed = rand(n, 1) * 2;         % Random speed between 0 and 2 m/s
            angle = rand(n, 1) * 2 * pi;    % Random angle between 0 and 2π
            V = [speed .* cos(angle), speed .* sin(angle)];
            
            % Initialize accelerations as zeros
            A = zeros(n, 2);
        end
        
        % Update node positions based on velocities and accelerations
        function [X_new, V_new, A_new] = update_node_positions(X, V, A, delta_t, arena_size)
            % Update velocities
            V_new = V + A * delta_t;
            
            % Update positions
            X_new = X + V_new * delta_t;
            
            % Update accelerations 
            A_new = 0.3 * A + 0.7 * (rand(size(A)) - 0.5) * 0.2;
            
            % Randomly change direction for some nodes
            for i = 1:size(X, 1)
                if rand() < 0.1  
                    angle_change = (rand() - 0.5) * pi;  
                    speed = norm(V_new(i, :));
                    current_angle = atan2(V_new(i, 2), V_new(i, 1));
                    new_angle = current_angle + angle_change;
                    V_new(i, :) = speed * [cos(new_angle), sin(new_angle)];
                end
            end
            
            % Bounce off arena boundaries
            for i = 1:size(X, 1)
                for j = 1:2
                    if abs(X_new(i, j)) > arena_size
                        % Bounce off the boundary
                        X_new(i, j) = sign(X_new(i, j)) * arena_size;
                        V_new(i, j) = -V_new(i, j);
                    end
                end
            end
        end
        
        % Generate distance matrix with noise
        function D = generate_distance_matrix(X, sigma_d)
            n = size(X, 1);
            D = zeros(n, n);
            
            for i = 1:n
                for j = i:n
                    if i == j
                        D(i, j) = 0;
                    else
                        % True distance
                        true_dist = norm(X(i, :) - X(j, :));
                        
                        % Add noise
                        D(i, j) = true_dist + normrnd(0, sigma_d);
                        D(j, i) = D(i, j);  
                    end
                end
            end
        end
        
        % Generate velocity measurements with noise
        function V_hat = generate_velocity_measurements(V, sigma_v)
            n = size(V, 1);
            V_hat = zeros(n, 2);
            
            for i = 1:n
                % Add noise to velocity magnitude
                speed = norm(V(i, :));
                noisy_speed = max(0, speed + normrnd(0, sigma_v));
                
                % Add minimal noise to direction (assuming compass with high accuracy)
                angle = atan2(V(i, 2), V(i, 1));
                noisy_angle = angle + normrnd(0, 0.017); 
                
                % Reconstruct velocity vector
                V_hat(i, :) = noisy_speed * [cos(noisy_angle), sin(noisy_angle)];
            end
        end
        
        % Calculate position estimation error
        function error = calculate_position_error(X_real, X_hat)
            
            % Calculate centroids
            centroid_real = mean(X_real, 1);
            centroid_hat = mean(X_hat, 1);
            
            % Center both point sets
            X_real_centered = X_real - repmat(centroid_real, size(X_real, 1), 1);
            X_hat_centered = X_hat - repmat(centroid_hat, size(X_hat, 1), 1);
            
            % Find optimal rotation using Procrustes analysis
            [~, ~, transform] = procrustes(X_real_centered, X_hat_centered, 'Scaling', false);
            
            % Apply the transformation to the estimated positions
            X_hat_aligned = transform.b * X_hat_centered * transform.T + repmat(centroid_real, size(X_hat, 1), 1);
            
            % Calculate the average Euclidean distance
            error = mean(sqrt(sum((X_real - X_hat_aligned).^2, 2)));
        end

        % Generate connectivity matrix for partially connected network
        function W = generate_connectivity_matrix(n, m)
            W = zeros(n, n);
            
            % First m nodes are fully connected with each other
            W(1:m, 1:m) = 1;
            
            % Other n-m nodes only communicate with the first m nodes
            W(m+1:n, 1:m) = 1;
            W(1:m, m+1:n) = 1;
            
            % Set diagonal elements to 0 (no self-connections)
            for i = 1:n
                W(i, i) = 0;
            end
        end
        
        % Generate distance matrix with noise, respecting connectivity constraints
        function D = generate_partial_distance_matrix(X, sigma_d, W)
            n = size(X, 1);
            D = zeros(n, n);
            
            for i = 1:n
                for j = 1:n
                    if i == j
                        % Set diagonal to zero
                        D(i, j) = 0;
                    elseif W(i, j) == 1
                        % True distance
                        true_dist = norm(X(i, :) - X(j, :));
                        
                        % Add noise while ensuring distance is non-negative
                        D(i, j) = max(0.1, true_dist + normrnd(0, sigma_d));
                        D(j, i) = D(i, j); % Ensure matrix is symmetric
                    else
                        % Set disconnected nodes' distances to NaN
                        D(i, j) = NaN;
                        D(j, i) = NaN;
                    end
                end
            end
        end
    end
end