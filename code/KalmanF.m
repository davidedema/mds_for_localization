classdef KalmanF
    methods (Static)
        % Initialize Kalman Filters for each node
        function kf = initialize_kalman_filters(n)
            kf = cell(n, 1);
            
            for i = 1:n
                kf{i}.x = zeros(4, 1);  % State: [x, y, vx, vy]'
                kf{i}.P = eye(4);       % State covariance
                
                % Process noise covariance
                kf{i}.Q = diag([0.1, 0.1, 0.5, 0.5]);
                
                % Measurement noise covariance
                kf{i}.R = diag([0.5, 0.5, 0.1, 0.1]);
                
                % State transition matrix
                kf{i}.F = [1 0 0 0;
                          0 1 0 0;
                          0 0 1 0;
                          0 0 0 1];
                
                % Measurement matrix
                kf{i}.H = eye(4);
            end
        end
        
        % Update Kalman Filter
        function kf = update_kalman_filter(kf, pos, vel, dt)
            % Update state transition matrix for the current dt
            kf.F = [1 0 dt 0;
                    0 1 0 dt;
                    0 0 1 0;
                    0 0 0 1];
            
            % Prediction step
            kf.x = kf.F * kf.x;
            kf.P = kf.F * kf.P * kf.F' + kf.Q;
            
            % Measurement
            z = [pos; vel'];
            
            % Update step
            y = z - kf.H * kf.x;            % Innovation
            S = kf.H * kf.P * kf.H' + kf.R; % Innovation covariance
            K = kf.P * kf.H' / S;           % Kalman gain
            
            kf.x = kf.x + K * y;            % Updated state estimate
            kf.P = (eye(4) - K * kf.H) * kf.P; % Updated covariance
        end
    end
end