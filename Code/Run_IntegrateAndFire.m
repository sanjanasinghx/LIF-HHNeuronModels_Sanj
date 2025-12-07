%% Run_IntegrateAndFire.m
% Driver script for leaky integrate-and-fire neuron
clear; close all; clc;

%% Simulation settings
dt    = 0.1;          % time step (ms)
T_end = 300;          % total duration (ms)
t     = 0:dt:T_end;   % time vector
nSteps = numel(t);

%% Input current: step pulse (50–200 ms)
I      = zeros(size(t));
t_on   = 50;          % ms, step onset
t_off  = 200;         % ms, step offset
I_amp  = 2;           % nA, step amplitude
I(t >= t_on & t <= t_off) = I_amp;

%% LIF parameters (biologically reasonable-ish)
params.E_L     = -70;      % mV, resting potential
params.V_reset = -65;      % mV, reset potential
params.V_th    = -55;      % mV, threshold
params.V_spike =  30;      % mV, spike peak (for plotting)
params.R_m     =  10;      % MΩ, membrane resistance
params.tau_m   =  10;      % ms, membrane time constant
params.n_ref   = round(2/dt);  % refractory steps (~2 ms)

%% Run the model
[t_sim, V_raw, spikes] = IntegrateAndFire(I, dt, params);

%% Plot results (publication-style)
figure(1); clf;
set(gcf,'Color','w');      % white figure background

%% Top: input current
subplot(2,1,1);
plot(t_sim, I, 'LineWidth', 1.5, 'Color', [0 0.2 0.8]);
ylabel('Input current (nA)');
title('Step Current Protocol', ...
      'FontWeight','bold', 'Color','k');
xlim([0 T_end]);
ylim([0 max(I_amp*1.2, 0.1)]);
grid on;
set(gca, ...
    'Box', 'on', ...
    'XColor', 'k', 'YColor', 'k', ...   % black axes
    'Color', 'w', ...                   % white axes background
    'FontSize', 11);

%% Bottom: membrane potential
subplot(2,1,2); hold on;

% Membrane potential trace
hV = plot(t_sim, V_raw, 'LineWidth', 1.5, 'Color', [0 0.6 0.2]);

% Resting & threshold lines
hRest = yline(params.E_L, '--', 'Resting', ...
    'Color', [0.3 0.3 0.3], ...
    'LabelHorizontalAlignment', 'left', ...
    'FontSize', 10);

hTh = yline(params.V_th, '--', 'Threshold', ...
    'Color', [0.8 0 0], ...
    'LabelHorizontalAlignment', 'left', ...
    'FontSize', 10);

% Mark spike times (red dots)
if any(spikes)
    spike_times = t_sim(spikes);
    spike_vals  = V_raw(spikes);
    hSpk = plot(spike_times, spike_vals, 'o', ...
        'MarkerFaceColor', [1 0.4 0.4], ...
        'MarkerEdgeColor', [0.8 0 0], ...
        'MarkerSize', 5);
else
    hSpk = [];  % in case no spikes
end

xlabel('Time (ms)');
ylabel('Membrane potential (mV)');
title('Leaky Integrate-and-Fire Neuron Response', ...
      'FontWeight','bold', 'Color','k');
xlim([0 T_end]);
ylim([-80 40]);
grid on;
set(gca, ...
    'Box', 'on', ...
    'XColor', 'k', 'YColor', 'k', ...   % black axes
    'Color', 'w', ...                   % white axes background
    'FontSize', 11);

% Legend (top-right) with white background
if isempty(hSpk)
    lgd = legend([hV hRest hTh], ...
        {'Membrane potential', 'Resting', 'Threshold'}, ...
        'Location', 'northeast');
else
    lgd = legend([hV hSpk hRest hTh], ...
        {'Membrane potential', 'Spikes', 'Resting', 'Threshold'}, ...
        'Location', 'northeast');
end
set(lgd, 'Box', 'on', 'Color', 'w');    % white legend box

hold off;
