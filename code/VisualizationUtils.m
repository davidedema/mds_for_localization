classdef VisualizationUtils
    methods (Static)

        % Visualize simulation results with animation
        function visualize_results(X_real, X_hat, errors)
            fig = figure('Position', [100, 100, 1200, 600]);
            
            % Number of nodes and time steps
            [num_steps, n, ~] = size(X_real);
            
            colors = hsv(n);
            
            ax1 = subplot(1, 2, 1);
            hold(ax1, 'on');
            grid(ax1, 'on');
            
            % Get bounds for the plot
            x_min = min(min(min(X_real(:,:,1))), min(min(X_hat(:,:,1)))) - 2;
            x_max = max(max(max(X_real(:,:,1))), max(max(X_hat(:,:,1)))) + 2;
            y_min = min(min(min(X_real(:,:,2))), min(min(X_hat(:,:,2)))) - 2;
            y_max = max(max(max(X_real(:,:,2))), max(max(X_hat(:,:,2)))) + 2;
            
            % Create axes for error plot
            ax2 = subplot(1, 2, 2);
            hold(ax2, 'on');
            grid(ax2, 'on');
            title(ax2, 'Position Estimation Error Over Time');
            xlabel(ax2, 'Time Step');
            ylabel(ax2, 'Average Error (meters)');
            error_line = plot(ax2, 1, errors(1), 'b-', 'LineWidth', 2);
            xlim(ax2, [1, num_steps]);
            ylim(ax2, [0, max(errors)*1.1]);
            
            % Add mean error text that will be updated
            error_text = text(ax2, 0.5, 0.95, 'Mean Error: calculating...', ...
                'Units', 'normalized', 'FontSize', 12, 'HorizontalAlignment', 'center');
            
            real_plots = cell(n, 1);
            est_plots = cell(n, 1);
            connect_lines = cell(n, 1);
            real_traces = cell(n, 1);
            est_traces = cell(n, 1);
            node_labels = cell(n, 1);
            
            % Specify number of past positions to show in trace
            trace_length = 10;
            
            % Initialize all plots
            for i = 1:n
                real_plots{i} = plot(ax1, NaN, NaN, 'o', 'MarkerSize', 10, 'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', 'k');
                est_plots{i} = plot(ax1, NaN, NaN, 's', 'MarkerSize', 8, 'MarkerFaceColor', colors(i,:));
                connect_lines{i} = line(ax1, [NaN, NaN], [NaN, NaN], 'Color', colors(i,:), 'LineStyle', '--');
                real_traces{i} = plot(ax1, NaN, NaN, '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
                est_traces{i} = plot(ax1, NaN, NaN, ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
                node_labels{i} = text(ax1, NaN, NaN, ['  Node ', num2str(i)], 'FontSize', 10);
            end
            
            % Add time step display
            time_text = text(ax1, 0.02, 0.98, 'Time step: 0', 'Units', 'normalized', 'FontSize', 12, ...
                'FontWeight', 'bold', 'BackgroundColor', [1 1 1 0.7]);
            
            % Add play/pause button - fixed implementation
            play_button = uicontrol('Style', 'pushbutton', 'String', 'Pause', ...
                'Position', [20 20 80 30]);
            is_playing = true;
            
            % Add slider for time control
            time_slider = uicontrol('Style', 'slider', 'Min', 1, 'Max', num_steps, ...
                'Value', 1, 'SliderStep', [1/(num_steps-1), 10/(num_steps-1)], ...
                'Position', [120 20 400 30]);
            
            % Add speed control
            speed_text = uicontrol('Style', 'text', 'String', 'Speed:', ...
                'Position', [540 20 50 30], 'HorizontalAlignment', 'left');
            speed_slider = uicontrol('Style', 'slider', 'Min', 0.1, 'Max', 5, ...
                'Value', 1, 'SliderStep', [0.1/4.9, 0.5/4.9], ...
                'Position', [590 20 100 30]);
            
            % Function to align estimated positions with real positions
            function X_aligned = align_positions(X_real_t, X_hat_t)
                % Calculate centroids
                centroid_real = mean(X_real_t, 1);
                centroid_hat = mean(X_hat_t, 1);
                
                % Center both point sets
                X_real_centered = X_real_t - repmat(centroid_real, size(X_real_t, 1), 1);
                X_hat_centered = X_hat_t - repmat(centroid_hat, size(X_hat_t, 1), 1);
                
                % Find optimal rotation using Procrustes analysis
                [~, ~, transform] = procrustes(X_real_centered, X_hat_centered, 'Scaling', false);
                
                % Apply the transformation to the estimated positions
                X_aligned = transform.b * X_hat_centered * transform.T + repmat(centroid_real, size(X_hat_t, 1), 1);
            end
            
            % Animation loop
            t = 1;
            
            % Set axis limits
            axis(ax1, [x_min, x_max, y_min, y_max]);
            title(ax1, 'Node Positions');
            xlabel(ax1, 'X Coordinate (meters)');
            ylabel(ax1, 'Y Coordinate (meters)');
            
            % Function to update the plot for a given time step
            function update_plot(t)
                % Update time step display
                set(time_text, 'String', ['Time step: ', num2str(t)]);
                
                % Get current data
                X_real_t = squeeze(X_real(t, :, :));
                X_hat_t = squeeze(X_hat(t, :, :));
                
                % Align estimated positions with real positions
                X_hat_aligned = align_positions(X_real_t, X_hat_t);
                
                % Update trajectories
                start_trace = max(1, t - trace_length);
                for i = 1:n
                    % Update current positions
                    set(real_plots{i}, 'XData', X_real_t(i, 1), 'YData', X_real_t(i, 2));
                    set(est_plots{i}, 'XData', X_hat_aligned(i, 1), 'YData', X_hat_aligned(i, 2));
                    
                    % Update connection lines
                    set(connect_lines{i}, 'XData', [X_real_t(i, 1), X_hat_aligned(i, 1)], ...
                        'YData', [X_real_t(i, 2), X_hat_aligned(i, 2)]);
                    
                    % Update labels
                    set(node_labels{i}, 'Position', [X_real_t(i, 1), X_real_t(i, 2)]);
                    
                    % Update trajectory traces
                    if t > 1
                        trace_t = start_trace:t;
                        
                        % Real trajectory
                        real_trace_x = squeeze(X_real(trace_t, i, 1));
                        real_trace_y = squeeze(X_real(trace_t, i, 2));
                        set(real_traces{i}, 'XData', real_trace_x, 'YData', real_trace_y);
                        
                        % Estimated trajectory - need to align each point
                        est_trace_x = zeros(length(trace_t), 1);
                        est_trace_y = zeros(length(trace_t), 1);
                        
                        for j = 1:length(trace_t)
                            time_idx = trace_t(j);
                            X_real_tj = squeeze(X_real(time_idx, :, :));
                            X_hat_tj = squeeze(X_hat(time_idx, :, :));
                            X_hat_aligned_tj = align_positions(X_real_tj, X_hat_tj);
                            est_trace_x(j) = X_hat_aligned_tj(i, 1);
                            est_trace_y(j) = X_hat_aligned_tj(i, 2);
                        end
                        
                        set(est_traces{i}, 'XData', est_trace_x, 'YData', est_trace_y);
                    end
                end
                
                % Update error plot
                set(error_line, 'XData', 1:t, 'YData', errors(1:t));
                
                % Update mean error text
                if t > 10
                    meanError = mean(errors(10:t));
                    set(error_text, 'String', ['Mean Error: ', num2str(meanError, '%.2f'), ' meters']);
                end
                
                % Update slider position
                set(time_slider, 'Value', t);
                
                % Refresh figure
                drawnow;
            end
            
            % Function to handle animation
            function animate()
                if ~isvalid(fig)
                    return;
                end
                
                if is_playing && t < num_steps
                    t = t + 1;
                    update_plot(t);
                    
                    % Get speed
                    speed = get(speed_slider, 'Value');
                    pause(0.1/speed);
                    
                    if is_playing && isvalid(fig)
                        % Queue the next iteration
                        drawnow;
                        pause(0.01); % Small pause to allow for button interactions
                        if is_playing && isvalid(fig)
                            animate();
                        end
                    end
                end
            end
            
            % Callback for the slider
            set(time_slider, 'Callback', @(src,~) slider_callback(src));
            function slider_callback(src)
                t = round(get(src, 'Value'));
                update_plot(t);
            end
            
            % Callback for the play/pause button
            set(play_button, 'Callback', @(src,~) play_callback(src));
            function play_callback(src)
                is_playing = ~is_playing;
                if is_playing
                    set(src, 'String', 'Pause');
                    animate();
                else
                    set(src, 'String', 'Play');
                end
            end
            
            % Initial plot
            update_plot(1);
            
            % Start animation
            animate();
        end

        function visualize_results_distributed(X_true, localizable, conn_matrix)
            % Visualize the results of the distributed MDS algorithm
            
            % Create figure
            figure('Position', [100, 100, 900, 500]);
            
            % Plot node connectivity
            subplot(1, 2, 1);
            hold on;
            grid on;
            title('Node Connectivity Graph');
            xlabel('X Coordinate (meters)');
            ylabel('Y Coordinate (meters)');
            
            % Plot connectivity edges
            n = size(X_true, 1);
            for i = 1:n
                for j = i+1:n
                    if conn_matrix(i, j) == 1
                        line([X_true(i,1), X_true(j,1)], ...
                             [X_true(i,2), X_true(j,2)], ...
                             'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);
                    end
                end
            end
            
            % Plot all nodes
            scatter(X_true(:,1), X_true(:,2), 80, 'b', 'filled');
            
            % Add node labels
            for i = 1:n
                text(X_true(i,1)+0.2, X_true(i,2), num2str(i), 'FontSize', 8);
            end
            
            % Plot localization results
            subplot(1, 2, 2);
            hold on;
            grid on;
            title('Distributed MDS Localization Results (Cluster Level)');
            xlabel('X Coordinate (meters)');
            ylabel('Y Coordinate (meters)');
            
            % Plot all nodes (gray for non-localized)
            scatter(X_true(:,1), X_true(:,2), 80, [0.7 0.7 0.7], 'filled');
            
            % Highlight localized nodes
            localized_nodes = find(localizable);
            scatter(X_true(localized_nodes,1), X_true(localized_nodes,2), 80, 'r', 'filled');
            
            % Add node labels
            for i = 1:n
                color = 'k';
                if localizable(i)
                    color = 'r';
                end
                text(X_true(i,1)+0.2, X_true(i,2), num2str(i), 'FontSize', 8, 'Color', color);
            end
            
            % Add legend
            legend_handles = zeros(2, 1);
            legend_handles(1) = plot(NaN, NaN, 'o', 'MarkerSize', 8, 'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerEdgeColor', 'none');
            legend_handles(2) = plot(NaN, NaN, 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'none');
            legend(legend_handles, {'Non-localizable Node', 'Localizable Node'}, 'Location', 'best');
            
            % Set axis properties
            subplot(1, 2, 1); axis equal;
            subplot(1, 2, 2); axis equal;
            
            % Create a second figure for network statistics
            figure('Position', [1000, 100, 500, 400]);
            
            % Compute node degree statistics
            degrees = sum(conn_matrix, 2);
            localized_degrees = degrees(localizable);
            nonlocalized_degrees = degrees(~localizable);
            
            % Create histogram of node degrees
            histogram(degrees, 'FaceColor', [0.7 0.7 0.7]);
            hold on;
            histogram(localized_degrees, 'FaceColor', 'r');
            
            title('Node Degree Distribution');
            xlabel('Node Degree');
            ylabel('Number of Nodes');
            legend('All Nodes', 'Localizable Nodes', 'Location', 'best');
            grid on;
            
            % Add text with statistics
            avg_degree = mean(degrees);
            text(0.05, 0.92, sprintf('Average degree: %.1f', avg_degree), ...
                'Units', 'normalized', 'FontSize', 12);
            
            if ~isempty(localized_degrees)
                avg_localized_degree = mean(localized_degrees);
                text(0.05, 0.85, sprintf('Avg. degree of localizable nodes: %.1f', avg_localized_degree), ...
                    'Units', 'normalized', 'FontSize', 12);
            end
            
            localizable_pct = 100 * sum(localizable) / length(localizable);
            text(0.05, 0.78, sprintf('Localizable nodes: %.1f%%', localizable_pct), ...
                'Units', 'normalized', 'FontSize', 12);
        end
        
        
        function visualize_global_results_distributed(X_true, global_nodes, global_positions, conn_matrix)
            % Create figure
            figure('Position', [100, 600, 900, 500]);
            
            % Plot node connectivity
            subplot(1, 2, 1);
            hold on;
            grid on;
            title('Node Connectivity Graph');
            xlabel('X Coordinate (meters)');
            ylabel('Y Coordinate (meters)');
            
            % Plot connectivity edges
            n = size(X_true, 1);
            for i = 1:n
                for j = i+1:n
                    if conn_matrix(i, j) == 1
                        line([X_true(i,1), X_true(j,1)], ...
                             [X_true(i,2), X_true(j,2)], ...
                             'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);
                    end
                end
            end
            
            % Plot all nodes
            scatter(X_true(:,1), X_true(:,2), 80, 'b', 'filled');
            
            % Add node labels
            for i = 1:n
                text(X_true(i,1)+0.2, X_true(i,2), num2str(i), 'FontSize', 8);
            end
            
            % Plot global localization results
            subplot(1, 2, 2);
            hold on;
            grid on;
            title('Global Coordinate System (Phase III)');
            xlabel('X Coordinate (meters)');
            ylabel('Y Coordinate (meters)');
            
            % Plot all nodes (gray for non-localized)
            scatter(X_true(:,1), X_true(:,2), 80, [0.7 0.7 0.7], 'filled');
            
            % Plot ground truth positions of localized nodes
            scatter(X_true(global_nodes,1), X_true(global_nodes,2), 80, 'r', 'filled');
            
            % Align estimated positions with ground truth for visualization
            X_est = zeros(size(X_true));
            X_est(global_nodes,:) = global_positions(global_nodes,:);
            
            % Find the Procrustes transformation between estimated and true positions
            [~, ~, transform] = procrustes(X_true(global_nodes,:), X_est(global_nodes,:), 'scaling', false);
            
            % Apply transformation to all estimated positions
            X_est_aligned = transform.b * X_est * transform.T + repmat(transform.c(1,:), size(X_est,1), 1);
            
            % Draw lines between true and estimated positions
            for i = 1:length(global_nodes)
                node = global_nodes(i);
                line([X_true(node,1), X_est_aligned(node,1)], ...
                     [X_true(node,2), X_est_aligned(node,2)], ...
                     'Color', 'g', 'LineWidth', 1);
            end
            
            % Plot estimated positions
            scatter(X_est_aligned(global_nodes,1), X_est_aligned(global_nodes,2), 50, 'g', 'filled');
            
            % Add node labels
            for i = 1:n
                color = 'k';
                if ismember(i, global_nodes)
                    color = 'r';
                end
                text(X_true(i,1)+0.2, X_true(i,2), num2str(i), 'FontSize', 8, 'Color', color);
            end
            
            % Add legend
            h1 = plot(NaN, NaN, 'o', 'MarkerSize', 8, 'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerEdgeColor', 'none');
            h2 = plot(NaN, NaN, 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'none');
            h3 = plot(NaN, NaN, 'o', 'MarkerSize', 8, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'none');
            h4 = plot(NaN, NaN, '-', 'Color', 'g', 'LineWidth', 1);
            
            legend([h1, h2, h3, h4], ...
                   {'Non-localized Node', 'True Position (Localized)', 'Estimated Position', 'Position Error'}, ...
                   'Location', 'best');
            
            % Set axis properties
            axis equal;
            
            % Calculate and display error metrics
            errors = sqrt(sum((X_true(global_nodes,:) - X_est_aligned(global_nodes,:)).^2, 2));
            mean_error = mean(errors);
            max_error = max(errors);
            
            % Add text with error statistics
            text(0.05, 0.05, sprintf('Mean error: %.2f meters', mean_error), ...
                'Units', 'normalized', 'FontSize', 10);
            text(0.05, 0.10, sprintf('Max error: %.2f meters', max_error), ...
                'Units', 'normalized', 'FontSize', 10);
            text(0.05, 0.15, sprintf('Nodes in global system: %d/%d', length(global_nodes), n), ...
                'Units', 'normalized', 'FontSize', 10);
                
            % Create a second figure for error distribution
            figure('Position', [1000, 600, 500, 400]);
            
            % Plot histogram of errors
            histogram(errors, 10, 'FaceColor', 'g');
            title('Distribution of Localization Errors');
            xlabel('Error (meters)');
            ylabel('Number of Nodes');
            grid on;
            
            % Add text with error statistics
            text(0.05, 0.92, sprintf('Mean error: %.2f meters', mean_error), ...
                'Units', 'normalized', 'FontSize', 12);
            text(0.05, 0.85, sprintf('Max error: %.2f meters', max_error), ...
                'Units', 'normalized', 'FontSize', 12);
            text(0.05, 0.78, sprintf('Nodes in global system: %d/%d (%.1f%%)', ...
                length(global_nodes), n, 100*length(global_nodes)/n), ...
                'Units', 'normalized', 'FontSize', 12);
        end
    end
end
