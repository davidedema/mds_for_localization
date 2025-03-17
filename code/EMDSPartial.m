classdef EMDSPartial
    methods (Static)
        function X = emds_partial(D, m, dim)
            n = size(D, 1);
            
            % First m nodes form the fully connected group (X1)
            % Remaining n-m nodes form the disconnected group (X2)
            
            % X1 distances
            D11 = D(1:m, 1:m);
            
            % Clear the matrix 
            for i = 1:m
                for j = 1:m
                    if isnan(D11(i,j))
                        D11(i,j) = 0; 
                    end
                end
            end
            
            % Projection type A for the fully connected group
            Im = eye(m);
            one_m = ones(m, 1);
            P_a = Im - (1/m)*(one_m*one_m');
            
            % Now classic MDS
            % Double centering for the first group 
            B11 = -0.5 * P_a * (D11.^2) * P_a;
            
            % Eigendecomposition with robust handling
            [V, E] = eig((B11+B11')/2);
            
            % Sort eigenvalues in descending order
            [evals, idx] = sort(diag(E), 'descend');
            V = V(:, idx);
            
            % Filter small eigenvalues to improve stability
            tol = max(abs(evals)) * 1e-10;
            valid_idx = abs(evals) > tol;
            
            % If we don't have enough eigenvectors, use the available ones
            if sum(valid_idx) < dim
                warning('Not enough significant eigenvalues. Using %d dimensions.', sum(valid_idx));
                V = V(:, valid_idx);
                E = diag(evals(valid_idx));
            else
                % Take the top 'dim' eigenvalues/vectors
                V = V(:, 1:dim);
                E = diag(evals(1:dim));
            end
            
            % Compute the coordinates for the first group
            X1 = V * sqrt(max(E,0)); 
            
            % Calculate center of gravity for first group
            x1_cg = (1/m) * sum(X1, 1);
            
            % Center the positions of the first group
            X1 = X1 - repmat(x1_cg, m, 1);
            
            % Projection type B for X2
            X2 = zeros(n-m, dim);
            
            for i = 1:n-m
                d_i = D(m+i, 1:m)';
                
                valid_distances = ~isnan(d_i);
                if sum(valid_distances) < 3
                    % Use a default position if not enough distances
                    X2(i,:) = [0 0];
                    continue;
                end
                
                % Use only valid distances
                d_i_valid = d_i(valid_distances);
                X1_valid = X1(valid_distances, :);
                
                % Use least squares to solve the position
                % Formulating as Ax = b
                A = -2 * X1_valid;
                b = d_i_valid.^2 - sum(X1_valid.^2, 2);
                
                % Add a regularization term for stability
                reg_lambda = 0.01;
                I_reg = reg_lambda * eye(size(A,2));
                A_reg = [A; I_reg];
                b_reg = [b; zeros(size(A,2), 1)];
                x = A_reg \ b_reg;
                
                % Store the result
                X2(i,:) = x';
            end
            
            % Combine the two groups
            X = [X1; X2];
        end
        
        % Enhanced MDS for partial connectivity
        function X_hat = enhanced_mds_partial(D_curr, X_prev, V_hat, delta_t, m)
            n = size(D_curr, 1);
            dim = size(X_prev, 2);
            
            % For fully connected first m nodes, use a modified version of the original enhanced_mds
            % Extract sub-matrices for the first group
            D_curr_sub = D_curr(1:m, 1:m);
            X_prev_sub = X_prev(1:m, :);
            V_hat_sub = V_hat(1:m, :);
            
            % Classical MDS for the first group
            X_mds = mds(D_curr_sub, dim);
            
            % Incorporate velocity information 
            X_hat_first = zeros(m, dim);
            for i = 1:m
                % Predicted position based on previous position and velocity
                pred_pos = X_prev_sub(i,:) + V_hat_sub(i,:) * delta_t;
                confidence = 0.7; % Adjust this weight based on relative confidence
                X_hat_first(i,:) = confidence * X_mds(i,:) + (1-confidence) * pred_pos;
            end
            
            % Projection type B
            X_hat_second = zeros(n-m, dim);
            
            for i = 1:n-m
                idx = i + m;
                d_i = D_curr(idx, 1:m);
                
                valid_idx = find(~isnan(d_i));
                if length(valid_idx) < 3 
                    if idx <= size(X_prev, 1)
                        X_hat_second(i,:) = X_prev(idx,:);
                    else
                        X_hat_second(i,:) = [0 0];
                    end
                    continue;
                end
                
                ref_idx = valid_idx(1);
                ref_dist_sq = d_i(ref_idx)^2;
                ref_pos = X_hat_first(ref_idx,:)';
                ref_pos_norm_sq = sum(ref_pos.^2);
                
                % Build equation system using differences
                A = zeros(length(valid_idx)-1, dim);
                b = zeros(length(valid_idx)-1, 1);
                
                for j = 2:length(valid_idx)
                    j_idx = valid_idx(j);
                    pos_j = X_hat_first(j_idx,:)';
                    
                    A(j-1,:) = 2 * (pos_j - ref_pos)';
                    b(j-1) = ref_dist_sq - d_i(j_idx)^2 + sum(pos_j.^2) - ref_pos_norm_sq;
                end
                
                % Regularizattion term
                reg_lambda = 0.01;
                A_reg = [A; reg_lambda * eye(dim)];
                
                if idx <= size(X_prev, 1)
                    b_reg = [b; reg_lambda * X_prev(idx,:)'];
                else
                    b_reg = [b; zeros(dim, 1)];
                end
                
                x = A_reg \ b_reg;
                
                X_hat_second(i,:) = x';
            end
            
            % Combine the two groups
            X_hat = [X_hat_first; X_hat_second];
        end
    end
end
