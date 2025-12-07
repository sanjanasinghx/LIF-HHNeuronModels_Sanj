%% Run_HodgkinHuxley.m
% Driver script for Hodgkin–Huxley neuron with multiple step currents.
% Produces a publication-style figure with:
%   (1) Input currents
%   (2) Membrane potential for 3 current levels

clear; close all; clc;

%% Simulation settings
dt    = 0.01;       % time step (ms)
T_end = 120;        % total duration (ms)
t     = 0:dt:T_end;
nSteps = numel(t);

% Step current window
t_on  = 10;         % ms
t_off = 90;         % ms

% Three current amplitudes (µA/cm^2)
I_amps = [4, 7, 10];   % subthreshold, regular spiking, stronger drive

% Preallocate current matrix: 3 x nSteps
I = zeros(3, nSteps);

for k = 1:3
    mask = (t >= t_on) & (t <= t_off);
    I(k, mask) = I_amps(k);
end

%% Run Hodgkin–Huxley model for each current
[t_sim, V1, m1, h1, n1] = HodgkinHuxley(I(1, :), dt); %#ok<ASGLU>
[~,     V2, m2, h2, n2] = HodgkinHuxley(I(2, :), dt); %#ok<ASGLU>
[~,     V3, m3, h3, n3] = HodgkinHuxley(I(3, :), dt); %#ok<ASGLU>

%% Plot results
figure(1); clf;
set(gcf, 'Color', 'w');   % white figure background

axisFont = 12;
titleFont = 13;

% --- (A) Input currents ---
subplot(2,1,1);
hold on;

plot(t_sim, I(1,:), 'LineWidth', 1.5);
plot(t_sim, I(2,:), 'LineWidth', 1.5);
plot(t_sim, I(3,:), 'LineWidth', 1.5);

xlabel('Time (ms)');
ylabel('Input current (µA/cm^2)');
title('Hodgkin–Huxley neuron: step current protocols');

legend({ ...
    sprintf('Subthreshold (%g µA/cm^2)', I_amps(1)), ...
    sprintf('Moderate (%g µA/cm^2)',    I_amps(2)), ...
    sprintf('Strong (%g µA/cm^2)',      I_amps(3))}, ...
    'Location', 'northoutside', 'Orientation', 'horizontal');

set(gca, 'FontSize', axisFont, 'Box', 'off', 'Layer', 'top');
xlim([0 T_end]);

% --- (B) Membrane potential responses ---
subplot(2,1,2);
hold on;

plot(t_sim, V1, 'LineWidth', 1.5);          % subthreshold
plot(t_sim, V2, 'LineWidth', 1.5);          % regular spiking
plot(t_sim, V3, 'LineWidth', 1.5);          % stronger current

% Resting potential reference line
yline(-65, '--', 'Resting', ...
    'Color', [0.5 0.5 0.5], ...
    'LabelHorizontalAlignment', 'left');

xlabel('Time (ms)');
ylabel('Membrane potential (mV)');
title('Hodgkin–Huxley membrane responses to different step currents');

legend({ ...
    sprintf('Subthreshold (%g µA/cm^2)', I_amps(1)), ...
    sprintf('Regular spiking (%g µA/cm^2)', I_amps(2)), ...
    sprintf('Stronger drive (%g µA/cm^2)', I_amps(3)), ...
    'Resting potential'}, ...
    'Location', 'northeast');

set(gca, 'FontSize', axisFont, 'Box', 'off', 'Layer', 'top');
xlim([0 T_end]);
grid on;
