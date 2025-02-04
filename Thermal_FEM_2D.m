% This script performs the following tasks:
% 1. Defines a 2D transient thermal model for a rectangular beam using the MATLAB PDE Toolbox.
% 2. Sets up thermal properties, initial conditions, and boundary conditions for the model.
% 3. Solves the transient thermal PDE for the specified time range and generates results.
% 4. Visualizes the results (temperature distribution and evolution over time) and exports data.
% 5. Extracts temperature data from specific sensor points, formats it, and expands it for PINN training.

%% create model
thermalmodel = createpde('thermal', 'transient');

% 0.16 * 0.2 rectangular (beam)
R1 = [3, 4, 0, 0.16, 0.16, 0, 0, 0, 0.02, 0.02]';
gd = R1;
sf = 'R1';
ns = char('R1');
ns = ns';
g = decsg(gd,sf,ns);
geometryFromEdges(thermalmodel,g);

figure;
pdegplot(thermalmodel, 'EdgeLabels', 'on');

%% Thermal propertiers and BC 
% properties
k = 116;      % thermal conductivity (W/(m*K))
rho = 2640;  % mass density (kg/m^3)
cp = 100;    % heat capacity (J/(kg*K))
thermalProperties(thermalmodel, 'ThermalConductivity', k, ...
                               'MassDensity', rho, ...
                               'SpecificHeat', cp);


% Boundary condition
h = 50;         % convection coefficient
T_inf = 22;     % ambient temperature

% Define boundary condition (trasient input here)
boundaryTemp = @(region, state) 57.74 - (57.74 - 33.14) * exp(-0.0051 * state.time);

thermalBC(thermalmodel, 'Edge', 4, 'Temperature', boundaryTemp); 
% thermalBC(thermalmodel, 'Edge', 3, 'Temperature', T_inf);
% thermalBC(thermalmodel, 'Edge', 2, 'Temperature', T_inf);
% thermalBC(thermalmodel, 'Edge', 1, 'Temperature', T_inf);


% thermalBC(thermalmodel, 'Edge', 4, 'ConvectionCoefficient', h, 'AmbientTemperature', T_inf);
thermalBC(thermalmodel, 'Edge', 3, 'ConvectionCoefficient', h, 'AmbientTemperature', T_inf);
thermalBC(thermalmodel, 'Edge', 1, 'ConvectionCoefficient', h, 'AmbientTemperature', T_inf);
thermalBC(thermalmodel, 'Edge', 2, 'ConvectionCoefficient', h, 'AmbientTemperature', T_inf);

figure;
pdegplot(thermalmodel, 'EdgeLabels', 'on');
%% check the temp input
% from 0-2000s
t = linspace(0, 2000, 20000);  

temp_t = 57.74 - (57.74 - 33.14) * exp(-0.0051 * t);

figure;
plot(t, temp_t, 'LineWidth', 2);
xlabel('Time (seconds)', 'FontSize', 14);
ylabel('Temperature (°C)', 'FontSize', 14);
title('Temperature vs Time', 'FontSize', 16);
grid on;

%% IC
% Initial condition (22 constant temp)
thermalIC(thermalmodel, 22);
% mesh size = 1e-3
generateMesh(thermalmodel, 'Hmax', 1e-3);  

figure;
pdeplot(thermalmodel);
%% solve the PDE

tlist = linspace(0,2000,2000);  
result = solve(thermalmodel, tlist);

%% plot final temp distribution 
figure;
pdeplot(thermalmodel, 'XYData', result.Temperature(:,end), 'Contour', 'off');
colormap(jet); 
colorbar; 
xlabel('X (m)', 'Fontsize', 16);
ylabel('Y (m)', 'Fontsize', 16);  
set(gca, 'Fontsize', 16)

%% video 

v = VideoWriter('FEM_baseline_2000s_h_50_Transient_input.mp4', 'MPEG-4'); 
v.FrameRate = 10; 
open(v);

figure;

for i = 1:length(tlist)/50
    set(gcf,'papertype','a4','paperorientation','portrait','paperunits','centimeters',...
    'paperposition',[0.63 0.63 28.41 19.72]);
    axis tight manual 
    pdeplot(thermalmodel, 'XYData', result.Temperature(:,i*50), 'Contour', 'off');
    title(['FEM Temperature at time = ', num2str(tlist(i*50))]);
    axis equal
    xlim([0,0.16]);
    ylim([0,0.02]);
    colormap(jet); 
    colorbar;  
    xlabel('X (m)', 'Fontsize', 16);
    ylabel('Y (m)', 'Fontsize', 16);
    set(gca, 'Fontsize', 10)
    
    frame = getframe(gcf);
    writeVideo(v, frame);
end

close(v);

disp('Video has been saved as FEM_baseline_2000s_h_50_Transient_input.mp4');

%% data export (down sampled)
% down sampled data time steps from 2k -> 200
% This is used for the comparision between FEM results and PINN prediction,
% since the resolution of the PINN prediction will be 200*200*200 (t*x*y)

num_time_steps = length(tlist);

% select 200 time steps
num_samples = 200;
sampled_time_indices = round(linspace(1, num_time_steps, num_samples));  
sampled_time = tlist(sampled_time_indices)';  

sampled_temperature = result.Temperature(:, sampled_time_indices);  

nodes = thermalmodel.Mesh.Nodes'; 

save('FEM_baseline_2000s_h_50_Transient_input_down_sampled.mat', 'nodes', 'sampled_temperature', 'sampled_time');

disp('Down sampled data has been saved as FEM_baseline_2000s_h_50_Transient_input_down_sampled.mat');

%% Extract temp information from FEM results 
% This is used for PINN data input
% Assume 3 points' FEM results as 'exp' data input 
% Extract 200*3 t,x,y data, construct them into 600*3 matrix
% Extract 200*3 T data, construct them into 600*1 matrix

% Reshape temperature data: result.Temperature is of size [n_nodes, n_timepoints]
% Extract mesh information
nodes = thermalmodel.Mesh.Nodes;
n_nodes = size(nodes, 2);  % Number of nodes
n_timepoints = length(tlist);
temperature_data = reshape(result.Temperature, n_nodes, n_timepoints);

% Define x and y coordinates for sampling
x_coords = [0.0165,0.08266, 0.146];
y_coord = 0.01;  % Constant y coordinate
dx = 0.001;  % Mesh resolution in x direction
dy = 0.001;  % Mesh resolution in y direction

% Find the nearest node indices corresponding to the x and y coordinates
x_indices = round(x_coords / dx) + 1;  % Indices for x coordinates
y_index = round(y_coord / dy) + 1;     % Index for y coordinate

% Now find the node indices for the specific x, y coordinates in the mesh
sampled_node_indices = [];
for x_i = x_coords
    dist = sqrt((nodes(1,:) - x_i).^2 + (nodes(2,:) - y_coord).^2);
    [~, min_index] = min(dist);
    sampled_node_indices = [sampled_node_indices, min_index];
end

% Initialize txy and T matrices for 200 time points
num_points = 200;
txy_data = zeros(num_points * 3, 3);  % [600, 3]: Each row contains [time, x, y] -> [t1 x1 y1; t2 x1 y1; ...]
T_data = zeros(num_points * 3, 1);    % [600, 1]: Each row contains the corresponding temperature value -> [T1; T2; ...]

% Sample 200 time points
sampled_time_indices = round(linspace(1, n_timepoints, num_points));  % 200 sampled time points

% Global index for storing data
global_index = 1;

% Extract temperature data for each x coordinate and sampled time points
for point = 1:3  % Iterate over 4 x coordinates
    for idx = 1:length(sampled_time_indices)
        time_step = sampled_time_indices(idx);
        txy_data(global_index, 1) = tlist(time_step);  % Time
        txy_data(global_index, 2) = x_coords(point);  % X coordinate
        txy_data(global_index, 3) = y_coord;  % Y coordinate (fixed)

        % Get corresponding temperature at this node and time
        T_data(global_index) = temperature_data(sampled_node_indices(point), time_step);

        % Update global index
        global_index = global_index + 1;
    end
end

% Save txy_data and T_data to .mat file
save('formatted_FEM_temperature_data_2000s_h_50_3_sensor_Transient_input.mat', 'txy_data', 'T_data');
fprintf('txy_data shape: [%d, %d]\n', size(txy_data, 1), size(txy_data, 2));
fprintf('T_data shape: [%d, 1]\n', length(T_data));

% Use the repmat function to repeat the matrix 10 times (PINN input requires 6000*3 input)
expanded_txy_data = repmat(txy_data, 10, 1);  % Expand txy_data to [6000, 3]
expanded_T_data = repmat(T_data, 10, 1);      % Expand T_data to [6000, 1]

% Verify the size of the expanded matrices
fprintf('expanded_txy_data shape: [%d, %d]\n', size(expanded_txy_data, 1), size(expanded_txy_data, 2));
fprintf('expanded_T_data shape: [%d, 1]\n', size(expanded_T_data, 1));

% Save the expanded data to a .mat file
save('expanded_FEM_temperature_data_2000s_h_50_3_sensor_Transient_input.mat', 'expanded_txy_data', 'expanded_T_data');




