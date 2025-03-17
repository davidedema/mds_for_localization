%% Classical MDS implementation 
function X = mds(D, dim)
    n = size(D, 1);
    
    % Replace any NaN with the mean of non-NaN entries
    D_mean = mean(D(~isnan(D)));
    D(isnan(D)) = D_mean;
    
    % Ensure matrix is symmetric
    D = (D + D')/2;
    
    % Double centering
    J = eye(n) - ones(n)/n;
    B = -0.5 * J * (D.^2) * J;
    
    % Make sure B is symmetric for numerical stability
    B = (B + B')/2;
    
    % Eigendecomposition with robust handling
    [V, E] = eig(B);
    
    % Sort eigenvalues in descending order
    [evals, idx] = sort(diag(E), 'descend');
    
    % Keep only positive eigenvalues with reasonable magnitude
    tol = max(abs(evals)) * 1e-10;
    valid_idx = evals > tol;
    
    % Check if we have enough dimensions
    if sum(valid_idx) < dim
        warning('Not enough positive eigenvalues for requested dimensions. Using %d dimensions.', sum(valid_idx));
        V = V(:, valid_idx);
        E = diag(evals(valid_idx));
    else
        V = V(:, idx(1:dim));
        E = diag(evals(1:dim));
    end
    
    % Compute the coordinates
    X = V * sqrt(E);
end