classdef RobustDistributedMDS
    
    methods (Static)
        
        function optimized_positions = optimize_cluster(nodes, initial_positions, D)
            
            n = length(nodes);
            
            % Extract the distance submatrix for this cluster
            D_cluster = zeros(n, n);
            for i = 1:n
                for j = 1:n
                    D_cluster(i, j) = D(nodes(i), nodes(j));
                end
            end
            
            X = initial_positions(nodes, :);
            
            stress_func = @(x) RobustDistributedMDS.compute_stress(x, D_cluster, n);
            
            options = optimset('Display', 'off', 'MaxIter', 100, 'TolFun', 1e-6);
            X_flat = reshape(X', [], 1);  
            
            % Perform optimization
            X_opt = fminsearch(stress_func, X_flat, options);
            
            X_opt_reshaped = reshape(X_opt, 2, n)';
            
            optimized_positions = initial_positions;
            optimized_positions(nodes, :) = X_opt_reshaped;
            
            return;
        end
        
        function stress = compute_stress(x, D, n)
            
            X = reshape(x, 2, n)';
            
            stress = 0;
            for i = 1:n
                for j = i+1:n
                    if ~isnan(D(i, j))
                        d_ij = norm(X(i, :) - X(j, :));
                        stress = stress + (d_ij - D(i, j))^2;
                    end
                end
            end
        end
        
        function transform = compute_transformation(source, target)
            
            % Center both point sets
            source_centroid = mean(source, 1);
            target_centroid = mean(target, 1);
            source_centered = source - repmat(source_centroid, size(source, 1), 1);
            target_centered = target - repmat(target_centroid, size(target, 1), 1);
            
            % Compute covariance matrix
            covariance = source_centered' * target_centered;
            
            % Compute SVD
            [U, ~, V] = svd(covariance);
            
            % Compute rotation matrix
            R = V * U';
            
            % Check for reflection
            if det(R) < 0
                V(:, 2) = -V(:, 2);
                R = V * U';
            end
            
            
            % Compute translation
            t = target_centroid' - R * source_centroid';
            
            transform = struct('rotation', R, 'translation', t', 'scale', scale);
        end
        
        function quads = find_robust_quads(node_idx, D, dmin)
            
            n = size(D, 1);
            quads = [];
            
            % Get neighbors of the central node
            neighbors = [];
            for j = 1:n
                if ~isnan(D(node_idx, j)) && j ~= node_idx
                    neighbors = [neighbors, j];
                end
            end
            
            % Check all possible combinations of 3 neighbors
            for i = 1:length(neighbors)
                node_a = neighbors(i);
                for j = i+1:length(neighbors)
                    node_b = neighbors(j);
                    % Check if nodeA and nodeB are neighbors
                    if isnan(D(node_a, node_b))
                        continue;
                    end
                    
                    for k = j+1:length(neighbors)
                        node_c = neighbors(k);
                        % Check if all nodes are connected to each other
                        if isnan(D(node_a, node_c)) || isnan(D(node_b, node_c))
                            continue;
                        end
                        
                        % Check if all triangles are robust
                        if RobustDistributedMDS.is_robust_triangle(node_idx, node_a, node_b, D, dmin) && ...
                           RobustDistributedMDS.is_robust_triangle(node_idx, node_a, node_c, D, dmin) && ...
                           RobustDistributedMDS.is_robust_triangle(node_idx, node_b, node_c, D, dmin) && ...
                           RobustDistributedMDS.is_robust_triangle(node_a, node_b, node_c, D, dmin)
                            quads = [quads; [node_idx, node_a, node_b, node_c]];
                        end
                    end
                end
            end
        end
        
        function is_robust = is_robust_triangle(node_a, node_b, node_c, D, dmin)
            
            % Get distances between nodes
            d_ab = D(node_a, node_b);
            d_ac = D(node_a, node_c);
            d_bc = D(node_b, node_c);
            
            % Find shortest side
            b = min([d_ab, d_ac, d_bc]);
            
            if d_ab > 0 && d_ac > 0 && d_bc > 0
                angle_a = acos((d_ab^2 + d_ac^2 - d_bc^2) / (2*d_ab*d_ac));
                angle_b = acos((d_ab^2 + d_bc^2 - d_ac^2) / (2*d_ab*d_bc));
                angle_c = acos((d_ac^2 + d_bc^2 - d_ab^2) / (2*d_ac*d_bc));
                
                theta = min([angle_a, angle_b, angle_c]);
                
                is_robust = (b * (sin(theta)^2) > dmin);
            else
                is_robust = false;
            end
        end
        
        function overlap_graph = create_overlap_graph(quads)
            
            num_quads = size(quads, 1);
            overlap_graph = zeros(num_quads, num_quads);
            
            % Check all pairs of quads
            for i = 1:num_quads
                for j = i+1:num_quads
                    quad1 = quads(i, :);
                    quad2 = quads(j, :);
                    
                    % Count shared nodes
                    shared_nodes = length(intersect(quad1, quad2));
                    
                    if shared_nodes == 3
                        overlap_graph(i, j) = 1;
                        overlap_graph(j, i) = 1;
                    end
                end
            end
        end
        
        function components = find_connected_components(graph)
            
            num_vertices = size(graph, 1);
            visited = false(1, num_vertices);
            components = {};
            
            for i = 1:num_vertices
                if ~visited(i)
                    % Start a new component
                    component = RobustDistributedMDS.depth_first_search(graph, i, visited);
                    visited(component) = true;
                    components{end+1} = component;
                end
            end
        end
        
        function component = depth_first_search(graph, start, visited)
            
            % Initialize
            stack = start;
            component = start;
            visited_local = visited;
            visited_local(start) = true;
            
            % DFS
            while ~isempty(stack)
                node = stack(end);
                stack(end) = [];
                
                % Find neighbors
                neighbors = find(graph(node, :));
                
                for i = 1:length(neighbors)
                    neighbor = neighbors(i);
                    if ~visited_local(neighbor)
                        visited_local(neighbor) = true;
                        stack = [stack, neighbor];
                        component = [component, neighbor];
                    end
                end
            end
        end
        
        function [positions, localized_nodes] = localize_component(central_node, quads, component, D)
            
            n = size(D, 1);
            positions = zeros(n, 2);
            localized_nodes = central_node;  % Start with central node
            
            % Start with the first quad in the component
            first_quad_idx = component(1);
            first_quad = quads(first_quad_idx, :);
            
            positions(central_node, :) = [0, 0];  
            
            first_neighbor = first_quad(2);
            positions(first_neighbor, :) = [D(central_node, first_neighbor), 0];
            localized_nodes = [localized_nodes, first_neighbor];
            
            second_neighbor = first_quad(3);
            d_origin_first = D(central_node, first_neighbor);
            d_origin_second = D(central_node, second_neighbor);
            d_first_second = D(first_neighbor, second_neighbor);
            
            % Using law of cosines to find the angle
            cos_angle = (d_origin_first^2 + d_origin_second^2 - d_first_second^2) / (2 * d_origin_first * d_origin_second);
            sin_angle = sqrt(1 - cos_angle^2);
            
            positions(second_neighbor, :) = [d_origin_second * cos_angle, d_origin_second * sin_angle];
            localized_nodes = [localized_nodes, second_neighbor];
            
            third_neighbor = first_quad(4);
            positions(third_neighbor, :) = RobustDistributedMDS.trilaterate(third_neighbor, [central_node, first_neighbor, second_neighbor], positions, D);
            localized_nodes = [localized_nodes, third_neighbor];
            
            % Initialize BFS queue with the first quad
            queue = component(1);
            visited = false(length(component), 1);
            visited(1) = true;
            
            % BFS through the overlap graph
            while ~isempty(queue)
                current_quad_idx = queue(1);
                queue(1) = [];
                
                % Find neighbors in overlap graph that haven't been visited
                for i = 1:length(component)
                    if visited(i)
                        continue;
                    end
                    
                    neighbor_quad_idx = component(i);
                    
                    % Check if these quads are connected (share 3 nodes)
                    curr_quad = quads(current_quad_idx, :);
                    neighbor_quad = quads(neighbor_quad_idx, :);
                    shared_nodes = intersect(curr_quad, neighbor_quad);
                    
                    if length(shared_nodes) == 3
                        unique_node = setdiff(neighbor_quad, shared_nodes);
                        
                        if ~ismember(unique_node, localized_nodes)
                            % Trilaterate its position using the three shared nodes
                            positions(unique_node, :) = RobustDistributedMDS.trilaterate(unique_node, shared_nodes, positions, D);
                            localized_nodes = [localized_nodes, unique_node];
                        end
                        
                        visited(i) = true;
                        queue = [queue, neighbor_quad_idx];
                    end
                end
            end
            
            localized_nodes = unique(localized_nodes);
        end
        
        function pos = trilaterate(node_to_localize, reference_nodes, known_positions, D)
            
            p1 = known_positions(reference_nodes(1), :);
            p2 = known_positions(reference_nodes(2), :);
            p3 = known_positions(reference_nodes(3), :);
            
            r1 = D(node_to_localize, reference_nodes(1));
            r2 = D(node_to_localize, reference_nodes(2));
            r3 = D(node_to_localize, reference_nodes(3));
            
            A = 2 * [p2(1)-p1(1), p2(2)-p1(2);
                     p3(1)-p1(1), p3(2)-p1(2)];
            
            b = [r1^2 - r2^2 - p1(1)^2 + p2(1)^2 - p1(2)^2 + p2(2)^2;
                 r1^2 - r3^2 - p1(1)^2 + p3(1)^2 - p1(2)^2 + p3(2)^2];
            
            pos = (A \ b)';
        end
    end
end