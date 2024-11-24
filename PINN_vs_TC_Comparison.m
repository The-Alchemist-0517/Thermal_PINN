% This script performs the following tasks:
% 1. Extracts temperature predictions from the PINN model at specific x-positions and a fixed y-position, corresponding to TC sensor locations.
% 2. Calculates statistical metrics (R² and RMSE) to quantify the agreement between the PINN predictions and TC measurements over time.
% 3. Visualizes the comparison between PINN predictions and TC data, as well as the absolute error, for each x-position over time.
% 4. Outputs overall R² and RMSE metrics for each x-position, summarizing model performance.
%% data import
datafiles = {
    '80_SS_160_2_face.csv';
    '80_SS_RTD2_2_face.csv';
    'RT-60_RTD.csv';
    'RT-60_TC.csv';
    'RT-80.csv';
    'RT-80_TC.csv';
    'RTD2_verify.csv';
    'Heater_exp.csv';};

%% Reading PINN and TC data
% Load results from the PINN model
pinn_data = load('PINN_exp_data_Trasient_input_1007_v.mat');
x_pinn = pinn_data.x;  % x-coordinate data (200x200)
y_pinn = pinn_data.y;  % y-coordinate data (200x200)
times_pinn = pinn_data.times;  % Time data
u_pinn = pinn_data.u;  % Temperature data (200x200xTime)

% Read TC (thermocouple) sensor data
tc_data = readtable(datafiles{8});
time_tc = tc_data.time(5:end) - tc_data.time(4);  % Adjust time data
tc_temps = [tc_data.TC1(5:end), tc_data.TC2(5:end), tc_data.TC3(5:end), tc_data.TC4(5:end)];  % TC temperature data

% Limit the number of time steps in the PINN results to 200
num_time_steps = 200;
time_indices_pinn = round(linspace(1, length(times_pinn), num_time_steps));
times_pinn_limited = times_pinn(time_indices_pinn);  % Restrict PINN time points to 200

% Find the closest TC data points to the PINN time points
[~, tc_time_indices] = arrayfun(@(t) min(abs(time_tc - t)), times_pinn_limited);
time_tc_limited = time_tc(tc_time_indices);  % TC time points closest to PINN time points
tc_temps_limited = tc_temps(tc_time_indices, :);  % Temperatures corresponding to those time points

% Define x positions to compare
x_values = [0.00641, 0.0491, 0.09013, 0.16];  % x positions corresponding to the TC sensors
y_fixed = 0.01;  % Fixed y-coordinate (0.01)

% Find the index in the y-coordinate closest to the fixed y=0.01
[~, y_index_pinn] = min(abs(y_pinn(:,1) - y_fixed));  % Assuming y_pinn values are the same in each column, use the first column

% Find the indices for each x-coordinate closest to the defined x_values
for i = 1:length(x_values)
    [~, x_indices_pinn(i)] = min(abs(x_pinn(1,:) - x_values(i)));  % Find the closest index for each x value
end

% Extract temperature data from the PINN model for the specified t, x, y
pinn_temps = zeros(length(x_values), num_time_steps);
for i = 1:length(x_values)
    for j = 1:num_time_steps
        pinn_temps(i, j) = u_pinn(time_indices_pinn(j), y_index_pinn, x_indices_pinn(i));
    end
end

%% Plot the comparision and absolute error(prediction error)
r2_time = zeros(length(time_tc_limited), length(x_values));  % Store R^2 for each time point
rmse_time = zeros(length(time_tc_limited), length(x_values)); % Store RMSE for each time point

for i = 1:length(x_values)
    % Extract TC and PINN temperature data
    tc_values = tc_temps_limited(:, i);   % TC temperature data
    pinn_values = pinn_temps(i, :) - 1;   % PINN temperature data (bias-corrected)
    
    % Calculate R^2 and RMSE at each time point
    for t = 1:length(time_tc_limited)
        % Residuals
        residuals = tc_values(t) - pinn_values(t);
        
        % RMSE (Root Mean Square Error)
        rmse_time(t, i) = abs(residuals); % Single-point absolute error (temperature difference)
        
        % R^2 Calculation
        ss_res = residuals.^2;            % Residual sum of squares
        ss_tot = var(tc_values) * (length(tc_values) - 1); % Total sum of squares
        r2_time(t, i) = 1 - (ss_res / ss_tot);            % R^2 formula
    end

    % Plot temperature comparison
    figure;
    hold on;
    plot(time_tc_limited, tc_values, 'LineWidth', 1, 'DisplayName', ['TC', num2str(i)]);
    plot(times_pinn_limited, pinn_values, '--', 'LineWidth', 1.5, 'DisplayName', ['PINN x = ', num2str(x_values(i))]);
    
    % Add legend and labels
    legend('show', 'Location', 'best');
    xlabel('Time (s)', 'FontSize', 20);
    ylabel('Temperature (°C)', 'FontSize', 20);
    title(['Comparison of TC', num2str(i), ' and PINN at x = ', num2str(x_values(i))]);
    set(gca, 'FontSize', 18);
    grid on;
    hold off;
    % Save plot
    saveas(gcf, ['Temperature_Comparison_TC', num2str(i), '_x_', num2str(x_values(i)), '.png']);
    
    % % Plot R^2 vs Time
    % figure;
    % plot(time_tc_limited, r2_time(:, i), 'LineWidth', 1.5, 'DisplayName', ['R^2 x = ', num2str(x_values(i))]);
    % xlabel('Time (s)', 'FontSize', 20);
    % ylabel('R^2', 'FontSize', 20);
    % title(['R^2 vs Time at x = ', num2str(x_values(i))]);
    % set(gca, 'FontSize', 18);
    % grid on;
    % % Save plot
    % saveas(gcf, ['R2_vs_Time_x_', num2str(x_values(i)), '.png']);
    
    % Plot Absolute Error (RMSE as absolute error per point)
    figure;
    plot(time_tc_limited, rmse_time(:, i), 'LineWidth', 1.5, 'DisplayName', ['Absolute Error x = ', num2str(x_values(i))]);
    xlabel('Time (s)', 'FontSize', 20);
    ylabel('Absolute Error (°C)', 'FontSize', 20);
    title(['Absolute Error vs Time at x = ', num2str(x_values(i))]);
    set(gca, 'FontSize', 18);
    grid on;
    % Save plot
    saveas(gcf, ['Absolute_Error_vs_Time_x_', num2str(x_values(i)), '.png']);
    
    % % Zoomed-in Absolute Error plot
    % figure;
    % plot(time_tc_limited, rmse_time(:, i), 'LineWidth', 1.5, 'DisplayName', ['Absolute Error x = ', num2str(x_values(i))]);
    % xlabel('Time (s)', 'FontSize', 20);
    % ylabel('Absolute Error (°C)', 'FontSize', 20);
    % title(['Absolute Error (Zoomed) vs Time at x = ', num2str(x_values(i))]);
    % xlim([800, 1200]); % Focus on a specific time range
    % set(gca, 'FontSize', 18);
    % grid on;
    % % Save plot
    % saveas(gcf, ['Absolute_Error_vs_Time_x_Zoomed_', num2str(x_values(i)), '.png']);
end

%% Compute overall R^2 and RMSE
% Initialize arrays to store R^2 and RMSE for each x position
r2_overall = zeros(length(x_values), 1); % Store overall R^2 for each x position
rmse_overall = zeros(length(x_values), 1); % Store overall RMSE for each x position

for i = 1:length(x_values)
    % Retrieve TC and PINN temperature data
    tc_values = tc_temps_limited(:, i);   % TC temperature data
    pinn_values = pinn_temps(i, :)';      % PINN temperature data, transposed to a column vector

    % Ensure data lengths match
    if length(tc_values) ~= length(pinn_values)
        error('Mismatch in the lengths of TC and PINN data');
    end

    % Compute residuals and associated metrics
    residuals = tc_values - pinn_values;               % Residuals at all time points
    ss_res = sum(residuals.^2);                        % Sum of squared residuals
    ss_tot = sum((tc_values - mean(tc_values)).^2);    % Total variance of TC data

    % Calculate overall R^2 and RMSE
    r2_overall(i) = 1 - (ss_res / ss_tot);             % Overall R^2
    rmse_overall(i) = sqrt(mean(residuals.^2));        % Overall RMSE

    % Print overall R^2 and RMSE for each x position
    fprintf('x = %.2f, Overall R^2: %.4f, Overall RMSE: %.4f\n', x_values(i), r2_overall(i), rmse_overall(i));
end


