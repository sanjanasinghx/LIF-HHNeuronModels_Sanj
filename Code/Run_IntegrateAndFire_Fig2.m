% Run_IntegrateAndFire_Fig2.m
% Subthreshold vs suprathreshold responses for the leaky integrate–and–fire neuron

clear; close all; clc;

%% Simulation settings
dt    = 0.1;          % time step (ms)
T_end = 300;          % total duration (ms)
t     = 0:dt:T_end;   % time vector
nSteps = numel(t);

%% Input currents: same step window, different amplitudes
t_on  = 50;           % ms, step onset
t_off = 200;          % ms, step offset

% Biophysically motivated: threshold current ~ 1.5 nA
I_sub_amp = 1.2;      % nA, subthreshold (no spikes)
I_sup_amp = 2.0;      % nA, suprathreshold (repetitive spiking)

I_sub = zeros(size(t));
I_sup = zeros(size(t));

I_sub(t >= t_on & t <= t_off) = I_sub_amp;
I_sup(t >= t_on & t <= t_off) = I_sup_amp;

% LIF parameters
params.E_L     = -70;        % mV, resting potential
params.V_reset = -65;        % mV
params.V_th    = -55;        % mV
params.V_spike = 30;         % mV
params.R_m     = 10;         % MΩ
params.tau_m   = 10;         % ms

% Refractory period (the missing field!)
params.n_ref   = round(2/dt);   % ~2 ms refractory steps

%% Run model for both currents
[~, V_sub,  spikes_sub]  = IntegrateAndFire(I_sub, dt, params);
[~, V_sup,  spikes_sup]  = IntegrateAndFire(I_sup, dt, params);

%% Figure 2: Subthreshold vs suprathreshold responses
figure(2); clf;
set(gcf, 'Color', 'w');    % white figure background

% Plot membrane potentials
plot(t, V_sub, 'LineWidth', 1.5, 'Color', [0 0.45 0.74]); hold on;  % blue-ish
plot(t, V_sup, 'LineWidth', 1.5, 'Color', [0.85 0.33 0.1]);        % orange-ish

% Mark spike times
if any(spikes_sub)
    plot(t(spikes_sub), V_sub(spikes_sub), 'o', ...
        'MarkerEdgeColor', [0 0.45 0.74], ...
        'MarkerFaceColor', 'w', ...
        'MarkerSize', 5);
end

if any(spikes_sup)
    plot(t(spikes_sup), V_sup(spikes_sup), 'o', ...
        'MarkerEdgeColor', [0.85 0.33 0.1], ...
        'MarkerFaceColor', 'w', ...
        'MarkerSize', 5);
end

% Reference lines: resting & threshold
yline(params.E_L,  '--', 'Resting',   ...
      'Color', [0.4 0.4 0.4], 'LabelHorizontalAlignment', 'left');
yline(params.V_th, '--', 'Threshold', ...
      'Color', [0.8 0   0  ], 'LabelHorizontalAlignment', 'left');

% Axes & labels
xlabel('Time (ms)');
ylabel('Membrane potential (mV)');
title('Subthreshold vs suprathreshold responses in LIF neuron');

xlim([0 T_end]);
ylim([-80 40]);
grid on;
set(gca, 'FontSize', 11, 'Box', 'on', 'XColor', 'k', 'YColor', 'k');

% Legend (top-right, white box, black text)
legend({'Subthreshold (1.2 nA)', ...
        'Suprathreshold (2.0 nA)', ...
        'Resting', ...
        'Threshold'}, ...
        'Location', 'northeast', ...
        'Box', 'on');

hold off;
