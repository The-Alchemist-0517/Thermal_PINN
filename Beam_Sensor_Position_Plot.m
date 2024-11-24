% This script visualizes the geometry of a 2D beam and the positions of sensors located on it. 

% Define the dimensions of the beam
Lx = 0.16;  % Length of the beam (m)
Ly = 0.02;  % Width of the beam (m)

% Sensor positions
x_coords = [0, 0.00641, 0.07002, 0.133];  % x-coordinates of the sensors
y_coord = 0.01;  % y-coordinate of the sensors (fixed at 0.01 m)

% Plot the beam shape
figure;
hold on;
set(gcf, 'papertype', 'a4', 'paperorientation', 'portrait', 'paperunits', 'centimeters', ...
    'paperposition', [0.63, 0.63, 28.41, 19.72]);
rectangle('Position', [0, 0, Lx, Ly], 'EdgeColor', 'b', 'LineWidth', 2);  % Draw the beam boundary
axis equal;  % Maintain equal scaling for x and y axes
xlabel('X (m)', 'FontSize', 16);
ylabel('Y (m)', 'FontSize', 16);
set(gca, 'FontSize', 16);

% Mark sensor positions on the beam
scatter(x_coords, repmat(y_coord, size(x_coords)), 100, 'r', 'filled');  % Plot sensor locations

% Add a horizontal dashed line to indicate y = 0.01
plot([0 Lx], [y_coord, y_coord], 'k--');  % Draw horizontal dashed line at y = 0.01
text(Lx + 0.005, y_coord, [num2str(y_coord)], 'HorizontalAlignment', 'left', 'FontSize', 14);  % Label y = 0.01

% Draw vertical dashed lines and annotate x-coordinates for each sensor
for i = 1:length(x_coords)
    % Draw a vertical dashed line from y = 0 to y = 0.01
    plot([x_coords(i), x_coords(i)], [0, y_coord], 'k--');  
    
    % Annotate x-coordinate at the bottom of the dashed line
    text(x_coords(i), -0.002, [num2str(x_coords(i))], 'HorizontalAlignment', 'center', 'FontSize', 12);
end

% Adjust the display limits
xlim([0 Lx]);
ylim([-0.005 Ly]);  % Leave space below y = 0 for x-coordinate labels

hold off;
