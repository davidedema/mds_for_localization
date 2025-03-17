function main_emds_distributed()
            % Simulation parameters
            n = 30;              % Number of nodes (increased for better connectivity)
            sigma_d = 0.3;       % Standard deviation for distance measurements (meters)
            arena_size = 15;     % Size of the arena (meters) (reduced for higher density)
            max_range = 10;      % Maximum communication range (meters) (increased)
            
            fprintf('Starting robust distributed MDS simulation with %d nodes...\n', n);
            
            % Initialize positions, velocities, and accelerations
            [X, ~, ~] = NodeUtils.initialize_nodes(n, arena_size);
            
            % Store true positions for comparison
            X_true = X;
            
            % Generate full distance matrix with noise
            D_full = NodeUtils.generate_distance_matrix(X, sigma_d);
            
            % Apply connectivity constraints based on maximum range
            D = D_full;
            conn_matrix = zeros(n, n);
            for i = 1:n
                for j = 1:n
                    if i ~= j && norm(X(i,:) - X(j,:)) <= max_range
                        conn_matrix(i,j) = 1;
                    else
                        D(i,j) = NaN;
                    end
                end
            end
            
            % Set robustness threshold based on noise
            dmin = 2.5 * sigma_d;  
            
            X_hat = cell(n, 1);             % Estimated positions in local coordinate systems
            localizable = false(n, 1);      % Nodes that can be localized
            cluster_sizes = zeros(n, 1);    % Size of each node's cluster
            
            fprintf('Phase I: Performing cluster localization at each node...\n');
            
            % Phase I: Cluster Localization
            for i = 1:n
                quads = RobustDistributedMDS.find_robust_quads(i, D, dmin);
                
                if isempty(quads)
                    continue;
                end
                
                fprintf('Node %d: Found %d robust quadrilaterals\n', i, size(quads, 1));
                
                overlap_graph = RobustDistributedMDS.create_overlap_graph(quads);
                
                components = RobustDistributedMDS.find_connected_components(overlap_graph);
                
                [~, largest_idx] = max(cellfun(@length, components));
                if ~isempty(components)
                    largest_component = components{largest_idx};
                
                    % Localize nodes in the largest component
                    [positions, localized_nodes] = RobustDistributedMDS.localize_component(i, quads, largest_component, D);
                    
                    localizable(localized_nodes) = true;
                    
                    % Store positions in local coordinate system
                    X_hat{i} = struct('positions', positions, 'nodes', localized_nodes);
                    
                    % Store cluster size
                    cluster_sizes(i) = length(localized_nodes);
                end
            end
            
            % Phase II: Cluster Optimization
            fprintf('\nPhase II: Optimizing cluster localizations...\n');
            
            for i = 1:n
                if cluster_sizes(i) > 0
                    positions = X_hat{i}.positions;
                    cluster_nodes = X_hat{i}.nodes;
                    
                    if length(cluster_nodes) >= 4
                        optimized_positions = RobustDistributedMDS.optimize_cluster(cluster_nodes, positions, D);
                        
                        X_hat{i}.positions = optimized_positions;
                        
                        fprintf('Optimized cluster for node %d with %d nodes\n', i, length(cluster_nodes));
                    end
                end
            end
            
            % Phase III: Cluster Transformation
            fprintf('\nPhase III: Computing transformations between clusters...\n');
            
            % Initialize global coordinate system
            [~, largest_cluster_idx] = max(cluster_sizes);
            
            if cluster_sizes(largest_cluster_idx) > 0
                fprintf('Using cluster of node %d as the global reference (size: %d nodes)\n', ...
                    largest_cluster_idx, cluster_sizes(largest_cluster_idx));
                
                % Global positions start with positions from the largest cluster
                global_positions = X_hat{largest_cluster_idx}.positions;
                global_nodes = X_hat{largest_cluster_idx}.nodes(:);
                nodes_in_global = false(n, 1);
                nodes_in_global(global_nodes) = true;
                
                % Keep track of which clusters have been merged
                processed_clusters = false(n, 1);
                processed_clusters(largest_cluster_idx) = true;
                
                % Repeatedly try to merge unprocessed clusters
                made_progress = true;
                while made_progress
                    made_progress = false;
                    
                    for i = 1:n
                        if cluster_sizes(i) == 0 || processed_clusters(i)
                            continue;
                        end
                        
                        cluster_nodes = X_hat{i}.nodes;
                        
                        common_nodes = find(nodes_in_global & ismember((1:n)', cluster_nodes));
                        
                        % Need at least 3 common nodes
                        if length(common_nodes) >= 3
                            local_positions = X_hat{i}.positions(common_nodes, :);
                            global_common_positions = global_positions(common_nodes, :);
                            
                            % Check if common nodes form a robust triangle for reliable transformation
                            robust_common = false;
                            if length(common_nodes) >= 3
                                for a = 1:length(common_nodes)-2
                                    for b = a+1:length(common_nodes)-1
                                        for c = b+1:length(common_nodes)
                                            node_a = common_nodes(a);
                                            node_b = common_nodes(b);
                                            node_c = common_nodes(c);
                                            if RobustDistributedMDS.is_robust_triangle(node_a, node_b, node_c, D, dmin)
                                                robust_common = true;
                                                break;
                                            end
                                        end
                                        if robust_common, break; end
                                    end
                                    if robust_common, break; end
                                end
                            end
                            
                            if robust_common
                                transform = RobustDistributedMDS.compute_transformation(local_positions, global_common_positions);
                                
                                new_nodes = setdiff(cluster_nodes, common_nodes);
                                
                                % Transform their positions to the global coordinate system
                                for j = 1:length(new_nodes)
                                    node = new_nodes(j);
                                    local_pos = X_hat{i}.positions(node, :);
                                    
                                    % Apply transformation
                                    global_positions(node, :) = transform.scale * local_pos * transform.rotation + transform.translation;
                                    
                                    nodes_in_global(node) = true;
                                    global_nodes = [global_nodes; node];
                                end
                                
                                processed_clusters(i) = true;
                                made_progress = true;
                                
                                fprintf('Merged cluster of node %d into global coordinate system\n', i);
                            end
                        end
                    end
                end
                
                fprintf('\nGlobal localization complete.\n');
                fprintf('Nodes in global coordinate system: %d/%d (%.1f%%)\n', ...
                    sum(nodes_in_global), n, 100*sum(nodes_in_global)/n);
                
                VisualizationUtils.visualize_global_results_distributed(X_true, global_nodes, global_positions, conn_matrix);
                
            else
                fprintf('No clusters were localized.\n');
            end
            
            num_localizable = sum(localizable);
            percent_localizable = 100 * num_localizable / n;
            avg_cluster_size = mean(cluster_sizes(cluster_sizes > 0));
            
            fprintf('\nCluster localization results:\n');
            fprintf('Nodes localizable: %d/%d (%.1f%%)\n', num_localizable, n, percent_localizable);
            fprintf('Average cluster size: %.1f nodes\n', avg_cluster_size);
            
            VisualizationUtils.visualize_results_distributed(X_true, localizable, conn_matrix);
        end