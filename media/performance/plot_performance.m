%% Tunnel performance with and without openmp
mesh = [100, 300, 500];
default = [201 2095 6667];
omp = [103 805 2342];

close all

plot(mesh, default, '-o', 'Color', [0.3 0.4 0.5], 'LineWidth', 2, ...
    'MarkerSize', 8, 'MarkerFaceColor', [0.3 0.4 0.5]); hold on;
plot(mesh, omp, '-s', 'Color', [0.9 0.4 0], 'LineWidth', 2, ...
    'MarkerSize', 8, 'MarkerFaceColor', [0.9 0.4 0]);
hold off;
legend( '-O3', '-O3 with OpenMP', location='northwest')
grid on;
grid minor
xlabel('Side size [Nodes]')
ylabel('Simulation time [s]')
title ('Simulation time on a variable size squared grid')

figure
plot(mesh, default./omp, '-s', 'Color', [0.9 0.4 0], 'LineWidth', 2, ...
    'MarkerSize', 8, 'MarkerFaceColor', [0.9 0.4 0])
grid on
grid minor 
xlabel('Side size [Nodes]')
title ('SpeedUp Graph')
legend('OpenMP vs. Default')

%% Lid driven performance with openmp and cuda
mesh = [100, 300, 700, 1000];
default = [20 167 931 1641];
omp = [12 88 420 746];
cuda = [1 7 41 84];


figure

plot(mesh, default, '-o', 'Color', [0.3 0.4 0.5], 'LineWidth', 2, ...
    'MarkerSize', 8, 'MarkerFaceColor', [0.3 0.4 0.5]); hold on;
plot(mesh, omp, '-s', 'Color', [0.9 0.4 0], 'LineWidth', 2, ...
    'MarkerSize', 8, 'MarkerFaceColor', [0.9 0.4 0]);
plot(mesh, cuda, '-s', 'Color', [0.2 0.5 0.3], 'LineWidth', 2, ...
    'MarkerSize', 8, 'MarkerFaceColor', [0.2 0.5 0.3]);
hold off;
legend( '-O3', '-O3 with OpenMP', 'CUDA' , location='northwest')
grid on;
grid minor
xlabel('Side size [Nodes]')
ylabel('Simulation time [s]')
title ('Simulation time on a variable size squared grid')

figure
plot(mesh, default./omp, '-s', 'Color', [0.9 0.4 0], 'LineWidth', 2, ...
    'MarkerSize', 8, 'MarkerFaceColor', [0.9 0.4 0])
hold on
plot(mesh, default./cuda, '-s', 'Color', [0.2 0.5 0.3], 'LineWidth', 2, ...
    'MarkerSize', 8, 'MarkerFaceColor', [0.2 0.5 0.3]);
hold off
grid on
grid minor
xlabel('Side size [Nodes]')
title ('SpeedUp Graph')
legend('OpenMP vs. Default', 'CUDA vs. Default', location='best')

%% Lid driven performance with openmp and cuda, varying time steps
time_steps = [2000, 5000, 10000, 15000];
default2 = [122 256 583 794];
omp2 = [43 116 287 482];
cuda2 = [5 12 25 38];

figure

plot(mesh, default2, '-o', 'Color', [0.3 0.4 0.5], 'LineWidth', 2, ...
    'MarkerSize', 8, 'MarkerFaceColor', [0.3 0.4 0.5]); hold on;
plot(mesh, omp2, '-s', 'Color', [0.9 0.4 0], 'LineWidth', 2, ...
    'MarkerSize', 8, 'MarkerFaceColor', [0.9 0.4 0]);
plot(mesh, cuda2, '-s', 'Color', [0.2 0.5 0.3], 'LineWidth', 2, ...
    'MarkerSize', 8, 'MarkerFaceColor', [0.2 0.5 0.3]);
hold off;
legend( '-O3', '-O3 with OpenMP', 'CUDA' , location='northwest')
grid on;
grid minor
xlabel('Number of iterations')
ylabel('Simulation time [s]')
title ('Simulation time on a 500x500 grid')

figure
plot(mesh, default2./omp2, '-s', 'Color', [0.9 0.4 0], 'LineWidth', 2, ...
    'MarkerSize', 8, 'MarkerFaceColor', [0.9 0.4 0])
hold on
plot(mesh, default2./cuda2, '-s', 'Color', [0.2 0.5 0.3], 'LineWidth', 2, ...
    'MarkerSize', 8, 'MarkerFaceColor', [0.2 0.5 0.3]);
hold off
grid on
grid minor
xlabel('Number of Iterations')
title ('SpeedUp Graph')
legend('OpenMP vs. Default', 'CUDA vs. Default', location='best')
