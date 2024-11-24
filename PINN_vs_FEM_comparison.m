% This script performs the following tasks:
% 1. Extracts temperature predictions from the PINN model at specific x-positions and a fixed y-position
% 2. Calculates statistical metrics (R² and RMSE) to quantify the agreement between the PINN predictions and FEM model
% 3. Outputs overall R² and RMSE metrics for each x,y-position, summarizing model performance.
%% Video for PINN Prediction
% Load the data
data = load('PINN_exp_data_Trasient_input_1007_v.mat');
x = data.x;
y = data.y;
u = data.u;
times = data.times;

% Convert x and y to row or column vectors
x_vec = linspace(0, 0.16, size(x, 2));  % Extend x-coordinates to cover 0-0.16m
y_vec = y(:, 1);  % Extract y-coordinates as a vector

% % Fix the temperature in the 20mm x 20mm region to 60°C
% for t = 1:length(times)
%     u(t, y_vec <= 0.02, x_vec <= 0.02) = 60;  % Set fixed region temperature
% end

% Create figure and axes
figure;
ax = gca;
% Set the axis to adjust for data range
axis tight manual 
ax.NextPlot = 'replaceChildren';
grid on; % Enable grid

% Preallocate video frames
frames(length(times)) = struct('cdata', [], 'colormap', []);

% Loop through time steps to plot the evolution of u
for t = 1:length(times)
    % Plot the 2D temperature field using imagesc
    set(gcf, 'papertype', 'a4', 'paperorientation', 'portrait', 'paperunits', 'centimeters', ...
        'paperposition', [0.63 0.63 28.41 19.72]);
    imagesc(x_vec, y_vec, squeeze(u(t, :, :)) - 3, 'Parent', ax);
    set(ax, 'yDir', 'normal'); % Set the y-axis direction upwards
    axis equal
    xlim(ax, [0, 0.16]); % Ensure x-axis range is 0-0.16m
    ylim(ax, [min(y_vec), max(y_vec)]); % Ensure appropriate y-axis range
    title(ax, sprintf('PINN Temperature at time = %.2f', times(t)));
    xlabel(ax, 'X (m)', 'Fontsize', 16);
    ylabel(ax, 'Y (m)', 'Fontsize', 16);
    set(gca, 'Fontsize', 10)
    colormap(ax, 'jet'); % Use 'jet' colormap; other options can also be used
    colorbar; % Display the colorbar
    grid(ax, 'on'); % Enable grid
    drawnow; 
    % Capture the current frame
    frames(t) = getframe(gcf);
end

% Save the animation using VideoWriter
v = VideoWriter('PINN_baseline_200s_3_points_Transient.mp4', 'MPEG-4');
v.FrameRate = 10; % Set the frame rate
open(v);
writeVideo(v, frames);
close(v);



%% Video for Difference 
% Load FEM results
data = load("C:\Users\theal\OneDrive\Desktop\DARPA_project\Thermal_pinn\baseline_0_200s_h_50_4_points\FEM_baseline_h_50_downsampled.mat");
nodes = data.nodes;  % Node coordinates (nx2 matrix)
temperature = data.sampled_temperature;  % Temperature at each node for all time steps (nxTimeSteps matrix)
time = data.sampled_time;  % Time vector

% Load predicted results (PINN)
pred_data = load("C:\Users\theal\OneDrive\Desktop\DARPA_project\Thermal_pinn\baseline_0_200s_h_50_4_points\PINN_baseline_h_50_downsampled.mat");
x = pred_data.x;  % x-coordinates (200x200 matrix)
y = pred_data.y;  % y-coordinates (200x200 matrix)
u = pred_data.u;  % Predicted temperature results (200x200x200 matrix)
times = pred_data.times;  % Time vector

% Select 200 uniformly spaced time points
num_samples = 200;
sampled_time_indices = round(linspace(1, length(time), num_samples));  % Uniformly sample 200 points
sampled_time = time(sampled_time_indices);  % Extract sampled time points
sampled_temperature = temperature(:, sampled_time_indices);  % Extract sampled temperature data

%% Downsampling Spatial Data (x, y) to a 200x200 Grid
% Create a new uniform 200x200 grid for x and y
x_grid = linspace(min(nodes(:,1)), max(nodes(:,1)), 200);  % 200 evenly spaced x-coordinates
y_grid = linspace(min(nodes(:,2)), max(nodes(:,2)), 200);  % 200 evenly spaced y-coordinates

% Generate the mesh grid
[X_grid, Y_grid] = meshgrid(x_grid, y_grid);

% Initialize storage for the downsampled temperature data (200x200x200)
downsampled_temperature = zeros(length(sampled_time),200, 200);

% Iterate over each time step to interpolate temperature data
for t_idx = 1:length(sampled_time)
    % Extract the temperature data for the current time step
    temp_data = sampled_temperature(:, t_idx);
    
    % Perform 2D interpolation using griddata to map node temperature data
    % to the new 200x200 grid
    downsampled_temperature(t_idx,:, :) = griddata(nodes(:,1), nodes(:,2), temp_data, X_grid, Y_grid, 'linear');
end


%% Calculate the Difference Matrix
difference = zeros(size(u));  % Initialize the difference matrix

% Iterate over each time step and grid point to calculate the difference between FEM and PINN
for t_idx = 1:length(sampled_time)
    for i = 1:size(X_grid, 1)
        for j = 1:size(X_grid, 2)
            % Calculate the difference between FEM and PINN predictions
            % The difference is expressed as the relative error
            difference(t_idx, j, i) = (downsampled_temperature(t_idx,j, i) - u(t_idx, j, i)) / downsampled_temperature(t_idx, j, i);
        end
    end
end

%% Create Video
v = VideoWriter('Difference_baseline_2000s_h_50.mp4', 'MPEG-4');  % Create an MPEG-4 video file
v.FrameRate = 10;  % Set the frame rate to 10 frames per second
open(v);

% Create figure
figure;

% Iterate over each time step to generate the difference distribution at each time point
for t_idx = 1:length(sampled_time)
    % Create 2D difference distribution plot
    set(gcf,'papertype','a4','paperorientation','portrait','paperunits','centimeters', ...
        'paperposition',[0.63 0.63 28.41 19.72]);
    
    % Plot the difference map
    imagesc(x_grid, y_grid, squeeze(difference(t_idx, :, :)));
    set(gca, 'YDir', 'normal');  % Set the y-axis direction to be upwards
    colormap(jet);  % Select the 'jet' colormap
    colorbar;  % Display the color bar
    title(sprintf('Difference at time = %.2f', sampled_time(t_idx)));
    axis equal;
    xlim([min(x_grid), max(x_grid)]);  % Ensure the x-axis range
    ylim([min(y_grid), max(y_grid)]);  % Ensure the y-axis range
    xlabel('X (m)', 'FontSize', 16);
    ylabel('Y (m)', 'FontSize', 16);
    set(gca, 'FontSize', 10);

    drawnow;
    
    % Capture the current frame and write to the video
    frame = getframe(gcf);
    writeVideo(v, frame);
end

% Close the video file
close(v);

disp('Difference evolution video has been saved as Difference_baseline_2000s_h_50.mp4');

%% 3D Error Distribution (Ensure the first dimension is time)
% Create grids for time, x, and y, ensuring dimensions are 200x200x200
[T_grid, X_grid, Y_grid] = ndgrid(sampled_time, x_grid, y_grid);

% Note: T_grid, X_grid, and Y_grid dimensions are [length(sampled_time), 200, 200]

% Flatten the grids into column vectors, ensuring dimensions are 200x200x200
T_flat = T_grid(:);  % Flatten to a 200*200*200 column vector
X_flat = X_grid(:);  % Flatten to a 200*200*200 column vector
Y_flat = Y_grid(:);  % Flatten to a 200*200*200 column vector

% Flatten the error matrix into a column vector, ensuring correct dimensions (time as the first dimension)
difference_flat = difference(:);  % Flatten to a 200*200*200 column vector in [time, x, y] order

% Ensure that the lengths of difference_flat, X_flat, Y_flat, and T_flat are all 200*200*200
assert(numel(T_flat) == numel(X_flat) && numel(X_flat) == numel(Y_flat) && numel(Y_flat) == numel(difference_flat), 'Flat vectors dimension mismatch');

% Create the figure
figure;

% Plot the 3D scatter plot with colors representing the error magnitude
scatter3(X_flat, Y_flat, T_flat, 20, difference_flat, 'filled');
clim([0, 0.05]);  % Set the color range from 0 to 0.05
% Set the colormap
colormap(jet);
colorbar;
xlabel('X (m)', 'FontSize', 16);
ylabel('Y (m)', 'FontSize', 16);
zlabel('Time (s)', 'FontSize', 16);
set(gca, 'FontSize', 12);

% Set graphical properties
grid on;
view(3);  % Set to 3D view

%% Test certain points 
% Obtain the absolute error, R², and RMSE for the four specific points (A, D, M, P)

% Select certain x and y coordinates for points A, D, M, P
selected_x = [0.08, 0.14];  % Example selected x coordinates (A, D, M, P)
selected_y = [0.01, 0.02];         % Example selected y coordinates (A, D, M, P)

% Store error evolution over time for these points
% error_over_time = zeros(length(sampled_time), length(selected_x), length(selected_y));

% Obtain the x and y coordinate vectors
x_vec = x(1, :);  % x-direction coordinate vector
y_vec = y(:, 1);  % y-direction coordinate vector

% Iterate through the selected x and y coordinates
for x_idx = 1:length(selected_x)
    for y_idx = 1:length(selected_y)
        % Find the indices of the grid points closest to the selected x and y coordinates
        [~, x_grid_idx(x_idx)] = min(abs(pred_data.x(1,:) - selected_x(x_idx)));
        [~, y_grid_idx(y_idx)] = min(abs(pred_data.y(:,1) - selected_y(y_idx)));
        
        % % Extract the error at this point for all time steps
        % for t_idx = 1:length(sampled_time)
        %     error_over_time(t_idx, x_idx, y_idx) = difference(t_idx, y_grid_idx, x_grid_idx);
        % end
    end
end

%% Calculate R² and RMSE for selected points

% Initialize R² and RMSE for the selected points
R2_selected_points = zeros(length(selected_x), length(selected_y));
RMSE_selected_points = zeros(length(selected_x), length(selected_y));

% Calculate R² and RMSE for each selected point
for x_idx = 1:length(selected_x)
    for y_idx = 1:length(selected_y)
        % Extract the errors at this point over time
        % error_at_point = error_over_time(:, x_idx, y_idx);
        
        % Extract FEM and PINN temperatures for this point
        FEM_temp_point = downsampled_temperature(:,x_grid_idx(x_idx), y_grid_idx(y_idx));
        PINN_temp_point = u(:, x_grid_idx(x_idx), y_grid_idx(y_idx));
        
        % Calculate total sum of squares (SST) for R²
        SST = sum((FEM_temp_point - mean(FEM_temp_point)).^2);
        
        % Calculate residual sum of squares (SSR) for R²
        SSR = sum((FEM_temp_point - PINN_temp_point).^2);
        
        % Calculate R²
        R2_selected_points(x_idx, y_idx) = 1 - (SSR / SST);
        
        % Calculate RMSE (Root Mean Square Error)
        MSE = mean((FEM_temp_point - PINN_temp_point).^2);
        RMSE_selected_points(x_idx, y_idx) = sqrt(MSE);
    end
end

disp('R² for selected points:');
disp(R2_selected_points);

disp('RMSE for selected points:');
disp(RMSE_selected_points);




