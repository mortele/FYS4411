%% Clean up
close all;
clear all;
clc;


%% Specify which file of data to load
inFileName = 'data_varying_a_and_b__N_100_points_from_0_5_to_3_0.dat';
%inFileName = 'data_varying_a_and_b__N_100_points_from_0_9_to_2_0.dat';


%% Load in data
% Determine the name of the data file.
folderName = '~/FYS4411/Part1/part1-build-Desktop-Release/';
inFile     = [folderName inFileName];

% Load data.
data = load(inFile);

% Reorganize data into separate arrays.
alpha           = data(:,1);
beta            = data(:,2);
E               = data(:,3);
acceptanceRatio = data(:,4);

% Find the global minima of E.
[minE, minE_index] = min(E);


%% Restructure E array as a matrix
% Find the dimensions N of the NxN grid of alpha and beta values.
i = 2;
while (alpha(i) > alpha(i-1)) 
    i = i + 1;
end
N = i-1;

% Shave the excess off alpha, beta, and E arrays.
alpha = alpha (1:end-(mod(length(alpha),N)));
beta  = beta  (1:end-(mod(length(beta ),N)));
E     = E     (1:end-(mod(length(E    ),N)));

% Reorganize E as a matrix.
Ematrix = zeros(N,length(alpha)/N);
for i=1:(length(alpha)/N)-1
    Ematrix(:,i) = E(i*N+1:(i+1)*N);
end


%% Plot data
% Set up the figure.
hFig = figure(1);
cameratoolbar('Show');
p1 = axes('OuterPosition', [0.171 0.25 0.677 0.614], 'units','normalized');
hold on;

% Plot the energy as a function of alpha and beta.
[a, b] = meshgrid(beta(1:floor(length(data(:,1))/N)),alpha(1:N));
scatter3(alpha,beta,E); %mesh(a, b, Ematrix);
xlabel('$\alpha$', 'interpreter', 'latex', 'FontSize', 15);
ylabel('$\beta$',  'interpreter', 'latex', 'FontSize', 15);
zlabel('$E$',      'interpreter', 'latex', 'FontSize', 15);

% Mark the minimum of the energy.  
plot3(  alpha(minE_index), beta(minE_index), E(minE_index),...
        'ro', 'MarkerSize', 15);
plot3(  alpha(minE_index), beta(minE_index), E(minE_index),...
        'r.', 'MarkerSize', 10);

% Set legends.
str  =  sprintf('$E_{min}=(%.3f,%.3f)=%f$',...
        alpha(minE_index),beta(minE_index),minE);
str2 =  sprintf('$|-2.9037 - E_{min}| = %.5f$', abs(-2.9037 - minE));
h    =  legend('$E(\alpha,\beta)$', str, str2);
set     (h, 'interpreter', 'latex', 'FontSize', 14);








