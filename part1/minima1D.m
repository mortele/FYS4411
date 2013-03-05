%% Clean up
close all;
clear all;
clc;


%% Load in data
% Determine the name of the data file.
folderName = '~/FYS4411/Part1/part1-build-Desktop-Release/';
inFileName = 'data_only_varying_a__N_1500_points_from_0_5_to_3_0.dat';
inFile     = [folderName inFileName];

data = load(inFile);


%% Reorganize data
alpha = data(:,1);
E     = data(:,3);


%% Plot data
figure(1);
hold on;

% Plot E as a function of alpha.
plot(alpha,E);
xlabel('$\alpha$', 'interpreter', 'latex','FontSize', 15);
ylabel('$E$', 'interpreter', 'latex','FontSize', 15);
legendStr = sprintf('MinE = %f', min(E));
legend(legendStr);

% Find and mark the global minia of E.
[minE minE_i] = min(E);
plot(alpha(minE_i), minE, 'ro', 'MarkerSize', 20);
plot(alpha(minE_i), minE, 'r.', 'MarkerSize', 10);


