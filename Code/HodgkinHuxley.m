function [t, V, m, h, n] = HodgkinHuxley(I_ext, dt)
% HodgkinHuxley.m
% Classic Hodgkin–Huxley single-compartment neuron with step current input.
%
% INPUTS
%   I_ext : row or column vector of external current (µA/cm^2)
%   dt    : time step (ms)
%
% OUTPUTS
%   t     : time vector (ms)
%   V     : membrane potential (mV)
%   m,h,n : gating variables for Na+ (m,h) and K+ (n)

    % Ensure row vector
    I_ext = I_ext(:)';      % 1 x N
    nSteps = numel(I_ext);
    t = (0:nSteps-1) * dt;

    % --- Biophysical parameters (classic HH squid axon) ---
    C_m  = 1.0;     % µF/cm^2
    g_Na = 120.0;   % mS/cm^2
    g_K  = 36.0;    % mS/cm^2
    g_L  = 0.3;     % mS/cm^2

    E_Na = 50.0;    % mV
    E_K  = -77.0;   % mV
    E_L  = -54.4;   % mV
    V_rest = -65.0; % mV

    % --- Allocate vectors ---
    V = zeros(1, nSteps);
    m = zeros(1, nSteps);
    h = zeros(1, nSteps);
    n = zeros(1, nSteps);

    % --- Initial conditions (standard HH starting point) ---
    V(1) = V_rest;
    m(1) = 0.05;
    h(1) = 0.60;
    n(1) = 0.32;

    % --- Time stepping (explicit Euler) ---
    for k = 1:nSteps-1
        Vk = V(k);
        mk = m(k);
        hk = h(k);
        nk = n(k);

        % ---- Voltage-dependent rate constants (ms^-1) ----
        % small "fixes" avoid 0/0 numerical issues
        if abs(Vk + 40) < 1e-6
            alpha_m = 1.0;
        else
            alpha_m = 0.1 * (Vk + 40) / (1 - exp(-(Vk + 40) / 10));
        end
        beta_m  = 4.0 * exp(-(Vk + 65) / 18);

        alpha_h = 0.07 * exp(-(Vk + 65) / 20);
        beta_h  = 1.0 ./ (1 + exp(-(Vk + 35) / 10));

        if abs(Vk + 55) < 1e-6
            alpha_n = 0.1;
        else
            alpha_n = 0.01 * (Vk + 55) / (1 - exp(-(Vk + 55) / 10));
        end
        beta_n  = 0.125 * exp(-(Vk + 65) / 80);

        % ---- Gating variable dynamics ----
        dm = alpha_m * (1 - mk) - beta_m * mk;
        dh = alpha_h * (1 - hk) - beta_h * hk;
        dn = alpha_n * (1 - nk) - beta_n * nk;

        m(k+1) = mk + dt * dm;
        h(k+1) = hk + dt * dh;
        n(k+1) = nk + dt * dn;

        % ---- Ionic currents ----
        I_Na = g_Na * (mk^3) * hk * (Vk - E_Na);
        I_K  = g_K  * (nk^4) * (Vk - E_K);
        I_L  = g_L  * (Vk - E_L);

        % ---- Membrane equation ----
        dV = ( I_ext(k) - I_Na - I_K - I_L ) / C_m;
        V(k+1) = Vk + dt * dV;
    end
end
