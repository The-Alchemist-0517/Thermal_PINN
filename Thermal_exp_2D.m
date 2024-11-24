% This script performs the following tasks:
% 1. Extracts and downsamples the temperature data from multiple RTDs, performing interpolation to create uniformly spaced data points for later use in PINN.
% 3. Fits an exponential decay model to the temperature data of RTD1(located in temp input), providing a smooth input for transient boundary conditions in the thermal model.
% 4. Generates several plots to visualize the temperature variation over time across different sensors and locations, including RTD and TC data comparisons.
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

%% Sample from RTD Data
% This is used for PINN data input
% Down sample RTD2 RTD3 RTD4 data (RTD1 is regarded as input) -> 2000
% Since this is the requirement of PINN input

data = readtable(datafiles{8});

% Define the number of sampling points for RTD2 and RTD3
num_samples = 2000;

% Uniformly sample RTD2 and RTD3
t_sampled_RTD2 = linspace(data.time(5), data.time(end), num_samples);  % Time sampling
T_sampled_RTD2 = interp1(data.time(5:end), data.RTD2(5:end), t_sampled_RTD2);  % Interpolate and sample RTD2

t_sampled_RTD3 = linspace(data.time(5), data.time(end), num_samples);  % Time sampling
T_sampled_RTD3 = interp1(data.time(5:end), data.RTD3(5:end), t_sampled_RTD3);  % Interpolate and sample RTD3

t_sampled_RTD4 = linspace(data.time(5), data.time(end), num_samples);  % Time sampling
T_sampled_RTD4 = interp1(data.time(5:end), data.RTD4(5:end), t_sampled_RTD4);  % Interpolate and sample RTD4

% Define the position coordinates (x and y), with y-coordinate fixed at 0.01
x_RTD2 = 0.00641;  % x-coordinate of RTD2
x_RTD3 = 0.07002;  % x-coordinate of RTD3
x_RTD4 = 0.133;    % x-coordinate of RTD4
y_coord = 0.01;    % y-coordinate fixed

% Create txy data for RTD2
txy_RTD2 = [t_sampled_RTD2', repmat(x_RTD2, num_samples, 1), repmat(y_coord, num_samples, 1)];

% Create txy data for RTD3
txy_RTD3 = [t_sampled_RTD3', repmat(x_RTD3, num_samples, 1), repmat(y_coord, num_samples, 1)];

% Create txy data for RTD4
txy_RTD4 = [t_sampled_RTD4', repmat(x_RTD4, num_samples, 1), repmat(y_coord, num_samples, 1)];

% Combine the txy data of RTD2, RTD3, and RTD4
txy_data = [txy_RTD2; txy_RTD3; txy_RTD4];

% Combine the temperature data of RTD2, RTD3, and RTD4
T_data = [T_sampled_RTD2'; T_sampled_RTD3'; T_sampled_RTD4'];

% Save the txy and T data into a MAT file
save('RTD_Temperature_Data_separate.mat', 'txy_data', 'T_data');

disp('RTD2, RTD3, and RTD4 sampled data have been saved in txy_data and T_data.');

%% Trasient input function fit
% Fit an exponential decay model to the temperature data from RTD1 over time

data = readtable(datafiles{8});

% Get the time and temperature data for RTD1
time_RTD1 = data.time(5:end);  % Starting from the 5th data point
temp_RTD1 = data.RTD1(5:end);  % Temperature data for RTD1

% Define the exponential fitting model (T(t) = a * exp(b * t) + c)
exp_fit_model = fittype('a*exp(b*t) + c', 'independent', 't', 'coefficients', {'a', 'b', 'c'});

% Set fitting options, define initial values and bounds for the parameters
fit_options = fitoptions('Method', 'NonlinearLeastSquares', ...
                         'StartPoint', [1, -0.001, 20], ...  % Initial values
                         'Lower', [-Inf, -Inf, -Inf], ...     % Parameter lower bounds
                         'Upper', [Inf, 0, Inf]);             % Parameter upper bounds

% Perform the nonlinear least squares fitting
[fit_result, gof] = fit(time_RTD1, temp_RTD1, exp_fit_model, fit_options);

% Generate fitted temperature data
fitted_temp_RTD1 = feval(fit_result, time_RTD1);

% Plot the actual data and the fitted curve
figure;
hold on;
plot(time_RTD1, temp_RTD1, 'r', 'LineWidth', 1.5, 'DisplayName', 'Real BC input');
plot(time_RTD1, fitted_temp_RTD1, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Fitted curve');
legend('show', 'Location', 'best');
xlabel('Time (s)', 'fontsize', 20);
ylabel('Temperature (°C)', 'fontsize', 20);
set(gca, 'fontsize', 18);
grid on;
hold off;

% Display the fitting result
disp('Exponential fitting result:');
disp(fit_result);

% Display the goodness of fit
disp('Goodness of fit:');
disp(gof);

%% RT-80_RTD temp Vs time

data = readtable(datafiles{8});

figure
hold on;

plot(data.time(5:end), data.RTD1(5:end), Color='r', LineWidth=1);
plot(data.time(5:end), data.RTD2(5:end), Color='g', LineWidth=1);
plot(data.time(5:end), data.RTD3(5:end), Color='b', LineWidth=1);
plot(data.time(5:end), data.RTD4(5:end), Color='k', LineWidth=1);
plot(data.time(5:end), data.RTD6(5:end), '--r', LineWidth=1);
plot(data.time(5:end), data.RTD5(5:end),'--g', LineWidth=1);
plot(data.time(5:end), data.RTD7(5:end), '--b', LineWidth=1);
plot(data.time(5:end), data.RTD8(5:end), '--k', LineWidth=1);

hold off;
legend({'RTD1 L1 ','RTD2 L1','RTD3 L1','RTD4 L1','RTD1 L2','RTD2 L2','RTD3 L2','RTD4 L2'}, location='best',NumColumns=2);
xlabel('Time (s)','fontsize', 20)
ylabel('Temperature (°C)','fontsize', 20)
set(gca,'fontsize', 16);

grid on;


%% RT-80_TC temp Vs time

data = readtable(datafiles{6});

figure
hold on;

% plot(data.time(2:end)-data.time(1), data.TC1(2:end), LineWidth=1);
% plot(data.time(2:end)-data.time(1), data.TC2(2:end), LineWidth=1);
% plot(data.time(2:end)-data.time(1), data.TC3(2:end), LineWidth=1);
% plot(data.time(2:end)-data.time(1), data.TC4(2:end), LineWidth=1);
plot(data.time(3:end)-data.time(1), data.TC5(3:end), LineWidth=1);
% plot(data.time(2:end)-data.time(1), data.TC6(2:end), LineWidth=1);

hold off;
legend({'TC1','TC2','TC3','TC4','TC5','TC6'}, location='best');
xlabel('Time (s)','fontsize', 20)
ylabel('Temperature (°C)','fontsize', 20)
set(gca,'fontsize', 18);

grid on;

%% compare 0 mm temp along depth
data_RTD = readtable(datafiles{5});  
data_TC = readtable(datafiles{6});  

figure
hold on;

plot(data_RTD.time(5:end), data_RTD.RTD1(5:end), 'Color', 'r', 'LineWidth', 1, 'DisplayName', 'RTD1');
plot(data_TC.time(2:end)-data_TC.time(1), data_TC.TC1(2:end), 'Color', 'g', 'LineWidth', 1, 'DisplayName', 'TC1');
plot(data_TC.time(2:end)-data_TC.time(1), data_TC.TC5(2:end), 'Color', 'b', 'LineWidth', 1, 'DisplayName', 'TC5');

legend('show', 'Location', 'best');
xlabel('Time (s)', 'fontsize', 20);
ylabel('Temperature (°C)', 'fontsize', 20);
set(gca, 'fontsize', 18);

grid on;
hold off;

%% plot 160 mm 2 faces' temp

data = readtable(datafiles{1});

figure
hold on;

plot(data.time(2:end), data.heater(2:end), LineWidth=1);
plot(data.time(2:end), data.back(2:end), LineWidth=1);

hold off;
legend({'TC-heater','TC-back'}, location='best');
xlabel('Time (s)','fontsize', 20)
ylabel('Temperature (°C)','fontsize', 20)
set(gca,'fontsize', 18);

grid on;

%% plot RTD 2 temp along depth

data_RTD = readtable(datafiles{7});  
data_TC = readtable(datafiles{2});   

figure
hold on;

plot(data_TC.time(2:end), data_TC.heater(2:end), 'Color', 'g', 'LineWidth', 1, 'DisplayName', 'TC1');
plot(data_RTD.time(5:end), data_RTD.RTD1(5:end), 'Color', 'r', 'LineWidth', 1, 'DisplayName', 'RTD1');
plot(data_TC.time(2:end), data_TC.back(2:end), 'Color', 'b', 'LineWidth', 1, 'DisplayName', 'TC5');

legend({'TC-heater','RTD2','TC-back'}, location='best');
xlabel('Time (s)', 'fontsize', 20);
ylabel('Temperature (°C)', 'fontsize', 20);
set(gca, 'fontsize', 18);

grid on;
hold off;





