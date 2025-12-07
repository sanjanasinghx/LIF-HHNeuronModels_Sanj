function [t, V, spikes] = IntegrateAndFire(I, dt, params)
% Clean and minimal leaky integrate-and-fire model

% Extract parameters
E_L   = params.E_L;      % resting potential (mV)
V_reset = params.V_reset;
V_th  = params.V_th;     % threshold
V_spike = params.V_spike;
tau_m = params.tau_m;    % membrane time constant (ms)
R_m   = params.R_m;      % membrane resistance (MΩ)
n_ref = params.n_ref;    % refractory steps

% Time vector comes from input I
N = length(I);
t = (0:N-1) * dt;

% Allocate
V = zeros(1, N);
V(1) = E_L;
spikes = false(1, N);

% Loop
ref_count = 0;

for i = 1:N-1
    
    if ref_count > 0
        V(i+1) = V_reset;
        ref_count = ref_count - 1;
        continue
    end

    % membrane ODE
    dVdt = (-(V(i) - E_L) + R_m * I(i)) / tau_m;
    V(i+1) = V(i) + dt * dVdt;

    % Spike condition
    if V(i+1) >= V_th
        V(i) = V_spike;
        V(i+1) = V_reset;
        spikes(i) = true;
        ref_count = n_ref;
    end
end
end
