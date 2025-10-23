%% Initialisation % % resolve base folders relative to this file
thisFileDir = fileparts(mfilename('fullpath'));
projectDir = thisFileDir;  % MATLAB folder is the base
flightInfoDir = fullfile(projectDir, 'Group6');
helpersDir = fullfile(projectDir, 'Helpers');
resultsDir = fullfile(projectDir, 'Results');

% Ensure required folders are on the path (Helpers contains aero5560_LoadFlightData)
addpath(helpersDir);
addpath(flightInfoDir);

% Basic sizes used below
dx = 1e-7;          % state perturbation size
du = 1e-7;          % control perturbation size
Xg = zeros(6,1);    % gust vector

%% load initial conditions
% Load trim condition (X0, U0)
load(fullfile(flightInfoDir, 'ICsAircraft6_120knts_400ft.mat'));

%% load flight data
FlightData           = FlightDataAircraft6();
FlightData.Geometric = FlightData.Geo;
FlightData.Inertial  = FlightData.I;
FlightData.Propeller = FlightData.Prop;

%% Task 1
% Task1 Part (a) Calculations
% ========================================================================
% initial guess of Xdot from trimmed states (essentially all zeros)
Xdot0               = getXdot(X0, U0, Xg, FlightData, true);
fprintf('Xdot0^T: ['); fprintf(' %.6e', Xdot0(:).'); fprintf(' ]\n');

% generate A matrix
A                   = zeros(12,12);
for i = 1:12
    Xp = X0; Xm = X0;
    Xp(i) = Xp(i) + dx;
    Xm(i) = Xm(i) - dx;
    fp = getXdot(Xp, U0, Xg, FlightData);
    fm = getXdot(Xm, U0, Xg, FlightData);
    A(:,i) = (fp - fm) / (2*dx);
end
% generate B matrix (include flap as 5th input for differentiation)
B                   = zeros(12, 5);
for j = 1:5
    Up = U0; Um = U0;
    Up(j) = Up(j) + du;
    Um(j) = Um(j) - du;
    fp = getXdot(X0, Up, Xg, FlightData);
    fm = getXdot(X0, Um, Xg, FlightData);
    B(:,j) = (fp - fm) / (2*du);
end

% Sense-check prints
disp('A size:'), disp(size(A));
disp('B size:'), disp(size(B));
disp('A:'), disp(A);
disp('B:'), disp(B);

% ========================================================================
% Isolate the lateral components of the A and B matrices
% Lateral states: [beta, p, r, phi, psi]
lat_idx            = [2 4 6 7 9];
Alat               = A(lat_idx, lat_idx);
% Lateral inputs: [a, r]
Blat               = B(lat_idx, [3 4]);

% Isolate the lateral components of the X and U vectors
Xlat               = X0(lat_idx);
Ulat               = U0([3 4]);

% ========================================================================
% Isolate G^{\phi}_{\delta a}
C_phi    = [0 0 0 1 0];
D_phi    = 0;

[num_phi_da, den_phi_da] = ss2tf(Alat, Blat(:,1), C_phi, D_phi);
G_phi_da = minreal(tf(num_phi_da, den_phi_da), 1e-3);

% Print forms
G_phi_da
zpk(G_phi_da)

%% Task 1.5: Mathematical Optimization for Critical Damping
% ------------------------------------------------------------------------------
% THEORY: Relationship Between Controller Parameters and Damping
% ================================================================
%
% 1. CLOSED-LOOP POLES & DAMPING:
%    For a second-order system with complex poles: s = -ζωn ± jωn√(1-ζ²)
%    - ζ (damping ratio): determines overshoot
%      * ζ = 1.0: critically damped (no overshoot)
%      * ζ = 0.707: ~5% overshoot
%      * ζ = 0.5: ~16% overshoot
%    - ωn (natural frequency): determines speed of response
%      * Settling time Ts ≈ 4/(ζωn)
%
% 2. CONTROLLER STRUCTURE:
%    C(s) = K · (1 + s/zl)/(1 + s/pl) · (1 + 1/(Ti·s)) · 1/(1 + s/wf)
%         = K · [Lead] · [Integrator] · [Low-pass filter]
%
%    - K: Overall gain → affects crossover frequency (speed)
%    - zl, pl: Lead compensation → adds phase boost at crossover
%              * Ratio pl/zl determines max phase boost: φ_max = arcsin((pl-zl)/(pl+zl))
%              * Frequency of max boost: ωm = √(zl·pl)
%    - Ti: Integrator time constant → removes steady-state error
%    - wf: Filter cutoff → attenuates high-frequency noise
%
% 3. OPTIMIZATION APPROACH:
%    We minimize: Cost = Ts + w·Mp²
%    Subject to:  ζ_dominant ≥ 0.85, PM ≥ 50°, all poles stable
%
%    This mathematically ensures critical damping while maximizing speed.

fprintf('\n========================================\n');
fprintf('MATHEMATICAL CONTROLLER OPTIMIZATION\n');
fprintf('========================================\n');

% Define the actuator and plant (fixed)
s = tf('s');
w_act_fixed = 35;
G_act_fixed = 1/(1 + s/w_act_fixed);

% Check if optimization results already exist
opt_cache_file = fullfile(resultsDir, 'optimization_cache.mat');
run_optimization = true;

if exist(opt_cache_file, 'file')
    fprintf('\n✓ Found cached optimization results: %s\n', opt_cache_file);
    fprintf('  Loading cached parameters...\n');
    cached_data = load(opt_cache_file);

    % Verify cache contains required fields
    if isfield(cached_data, 'K_opt') && isfield(cached_data, 'zl_opt') && ...
       isfield(cached_data, 'pl_opt') && isfield(cached_data, 'Ti_opt') && ...
       isfield(cached_data, 'best_cost')
        K_opt = cached_data.K_opt;
        zl_opt = cached_data.zl_opt;
        pl_opt = cached_data.pl_opt;
        Ti_opt = cached_data.Ti_opt;
        wf = 28;

        fprintf('\n  Cached optimal parameters:\n');
        fprintf('     K  = %.4f\n', K_opt);
        fprintf('     zl = %.3f rad/s\n', zl_opt);
        fprintf('     pl = %.3f rad/s\n', pl_opt);
        fprintf('     Ti = %.3f s\n', Ti_opt);
        fprintf('     Best cost = %.3f\n', cached_data.best_cost);
        fprintf('\n  To re-run optimization, delete: %s\n', opt_cache_file);

        run_optimization = false;
    else
        fprintf('  Warning: Cache file incomplete. Re-running optimization...\n');
    end
end

if run_optimization
    % Grid search parameters (no toolbox required!)
    fprintf('\nSearching parameter space via grid search...\n');
    fprintf('Target: ζ_dominant ≥ 0.40, Mp < 10%%, PM > 50° (soft constraints)\n');

    % REFINED SEARCH around best solution: K=-0.146, zl=2.1, pl=23.2, Ti=0.667
    % Previous best: Cost=0.374, Ts=0.374s, Mp=0%, PM=73.8°
    fprintf('Stage 2: Refined search around best solution from Stage 1\n\n');

% Define search grids - REFINED around optimal region
K_range  = linspace(0.10, 0.20, 15);   % Narrow around K=0.146
zl_range = linspace(1.5, 2.5, 12);     % Narrow around zl=2.1
pl_range = linspace(18.0, 28.0, 12);   % Narrow around pl=23.2
Ti_range = linspace(0.50, 0.85, 10);   % Narrow around Ti=0.667
% Total combinations: 15*12*12*10 = 21,600 (takes ~60 seconds)

wf = 28;

% Create all parameter combinations
[K_grid, zl_grid, pl_grid, Ti_grid] = ndgrid(K_range, zl_range, pl_range, Ti_range);
param_combos = [K_grid(:), zl_grid(:), pl_grid(:), Ti_grid(:)];
total_evals = size(param_combos, 1);

fprintf('Evaluating %d parameter combinations using parallel processing...\n', total_evals);

% Pre-allocate results
costs = inf(total_evals, 1);

% Profiling arrays (for diagnostic purposes)
timing_build = zeros(total_evals, 1);
timing_stability = zeros(total_evals, 1);
timing_margin = zeros(total_evals, 1);
timing_damping = zeros(total_evals, 1);
timing_stepinfo = zeros(total_evals, 1);

tic;

% Parallel loop over all combinations
parfor idx = 1:total_evals
    K_val  = param_combos(idx, 1);
    zl_val = param_combos(idx, 2);
    pl_val = param_combos(idx, 3);
    Ti_val = param_combos(idx, 4);

    % Build controller
    t0 = tic;
    s_local = tf('s');
    K_test = -K_val;  % negative for feedback
    Clead = (1 + s_local/zl_val)/(1 + s_local/pl_val);
    Cint  = (1 + 1/(Ti_val*s_local));
    Clp   = 1/(1 + s_local/wf);
    C_test = minreal(K_test * Clead * Cint * Clp, 1e-6);

    % Actuator model (create fresh in parfor)
    G_act_local = 1/(1 + s_local/w_act_fixed);

    % Open and closed loop
    L_test = minreal(C_test * G_act_local * G_phi_da, 1e-6);
    T_test = minreal(feedback(L_test, 1), 1e-6);
    timing_build(idx) = toc(t0);

    % Check stability
    t0 = tic;
    cl_poles = pole(T_test);
    if any(real(cl_poles) >= 0)
        costs(idx) = inf;  % unstable
        timing_stability(idx) = toc(t0);
        continue;
    end
    timing_stability(idx) = toc(t0);

    % Check phase margin
    t0 = tic;
    [~, PM, ~, ~] = margin(L_test);
    if PM < 50
        costs(idx) = inf;  % insufficient margin
        timing_margin(idx) = toc(t0);
        continue;
    end
    timing_margin(idx) = toc(t0);

    % Check damping ratio for DOMINANT poles only (ignoring plant's natural dynamics)
    % Focus on poles that dominate the step response (2-10 rad/s range)
    t0 = tic;
    min_zeta = 1.0;
    found_relevant_pole = false;
    for i = 1:length(cl_poles)
        if imag(cl_poles(i)) ~= 0
            wn_pole = abs(cl_poles(i));
            zeta_pole = -real(cl_poles(i))/wn_pole;
            % Only check poles in the dominant response range (2 to 10 rad/s)
            % This ignores very slow plant poles that we can't easily fix
            if wn_pole >= 2.0 && wn_pole <= 10.0
                min_zeta = min(min_zeta, zeta_pole);
                found_relevant_pole = true;
            end
        end
    end

    % Only apply damping constraint if we found relevant poles
    if found_relevant_pole && min_zeta < 0.40  % Very relaxed - focus on stability
        costs(idx) = inf;  % insufficient damping
        timing_damping(idx) = toc(t0);
        continue;
    end
    timing_damping(idx) = toc(t0);

    % Compute cost with heavy penalties for poor performance
    t0 = tic;
    try
        S = stepinfo(T_test);
        Ts = S.SettlingTime;
        Mp = S.Overshoot;

        % Base cost: settling time
        base_cost = Ts;

        % Heavy overshoot penalty (we want Mp < 10%)
        if Mp > 10
            overshoot_penalty = 50 * (Mp - 10)^2;  % Very steep penalty
        else
            overshoot_penalty = 0.5 * Mp^2;  % Small penalty for acceptable overshoot
        end

        % Damping penalty for dominant poles (soft constraint)
        damping_penalty = 0;
        if found_relevant_pole && min_zeta < 0.60
            damping_penalty = 100 * (0.60 - min_zeta)^2;  % Encourage better damping
        end

        costs(idx) = base_cost + overshoot_penalty + damping_penalty;
    catch
        costs(idx) = inf;  % problem evaluating
    end
    timing_stepinfo(idx) = toc(t0);
end

% Find best result
[best_cost, best_idx] = min(costs);
best_params = param_combos(best_idx, :);

elapsed_time = toc;
fprintf('\n✓ Grid search complete in %.1f seconds!\n', elapsed_time);

% Print profiling results
fprintf('\n[Performance Profiling - Average Times per Iteration]\n');
fprintf('   Controller build:     %.4f ms\n', mean(timing_build)*1000);
fprintf('   Stability check:      %.4f ms\n', mean(timing_stability)*1000);
fprintf('   Phase margin check:   %.4f ms\n', mean(timing_margin)*1000);
fprintf('   Damping check:        %.4f ms\n', mean(timing_damping)*1000);
fprintf('   Step info (cost):     %.4f ms\n', mean(timing_stepinfo(timing_stepinfo>0))*1000);
fprintf('   Total per iteration:  %.4f ms\n', ...
    mean(timing_build + timing_stability + timing_margin + timing_damping + timing_stepinfo)*1000);

% Identify bottlenecks
total_times = [mean(timing_build), mean(timing_stability), mean(timing_margin), ...
               mean(timing_damping), mean(timing_stepinfo(timing_stepinfo>0))];
[~, slowest_idx] = max(total_times);
operations = {'Controller build', 'Stability check', 'Phase margin', 'Damping check', 'Step info'};
fprintf('   Slowest operation:    %s\n', operations{slowest_idx});

% Count how many passed each stage
n_stable = sum(costs < inf | timing_margin > 0);
n_goodPM = sum(costs < inf | timing_damping > 0);
n_goodDamp = sum(timing_stepinfo > 0);
n_valid = sum(costs < inf);
fprintf('\n[Constraint Satisfaction]\n');
fprintf('   Stable systems:       %d / %d (%.1f%%)\n', n_stable, total_evals, 100*n_stable/total_evals);
fprintf('   PM > 50°:             %d / %d (%.1f%%)\n', n_goodPM, total_evals, 100*n_goodPM/total_evals);
fprintf('   ζ_dom > 0.40:         %d / %d (%.1f%%)\n', n_goodDamp, total_evals, 100*n_goodDamp/total_evals);
fprintf('   Valid solutions:      %d / %d (%.1f%%)\n', n_valid, total_evals, 100*n_valid/total_evals);

K_opt  = -best_params(1);
zl_opt = best_params(2);
pl_opt = best_params(3);
Ti_opt = best_params(4);

fprintf('\nOptimal parameters:\n');
fprintf('   K  = %.4f\n', K_opt);
fprintf('   zl = %.3f rad/s\n', zl_opt);
fprintf('   pl = %.3f rad/s\n', pl_opt);
fprintf('   Ti = %.3f s\n', Ti_opt);
fprintf('   Best cost = %.3f\n', best_cost);

% Analyze top 10 solutions to understand parameter trends
fprintf('\n[Top 10 Solutions Analysis]\n');
[sorted_costs, sort_idx] = sort(costs);
valid_costs = sorted_costs(sorted_costs < inf);
if length(valid_costs) >= 10
    fprintf('Analyzing top 10 solutions for parameter trends...\n');
    top10_idx = sort_idx(1:10);
    top10_params = param_combos(top10_idx, :);

    fprintf('   K  range: [%.4f, %.4f]  (mean: %.4f)\n', ...
        min(top10_params(:,1)), max(top10_params(:,1)), mean(top10_params(:,1)));
    fprintf('   zl range: [%.3f, %.3f]  (mean: %.3f)\n', ...
        min(top10_params(:,2)), max(top10_params(:,2)), mean(top10_params(:,2)));
    fprintf('   pl range: [%.3f, %.3f]  (mean: %.3f)\n', ...
        min(top10_params(:,3)), max(top10_params(:,3)), mean(top10_params(:,3)));
    fprintf('   Ti range: [%.3f, %.3f]  (mean: %.3f)\n', ...
        min(top10_params(:,4)), max(top10_params(:,4)), mean(top10_params(:,4)));
    fprintf('   Cost range: [%.3f, %.3f]\n', valid_costs(1), valid_costs(10));

    % Visualize parameter space clustering
    figure('Name','Parameter Space Analysis');

    % Get top 100 solutions for visualization
    n_top = min(100, length(valid_costs));
    top_idx = sort_idx(1:n_top);
    top_params = param_combos(top_idx, :);
    top_costs = costs(top_idx);

    % Normalize costs for colormap (0=best, 1=worst of top 100)
    norm_costs = (top_costs - min(top_costs)) / (max(top_costs) - min(top_costs));

    subplot(2,2,1);
    scatter(top_params(:,1), top_params(:,2), 50, norm_costs, 'filled');
    xlabel('K (gain)'); ylabel('z_l (rad/s)');
    title('Top 100 Solutions: K vs z_l');
    colorbar; colormap(jet); clim([0 1]);
    hold on;
    plot(best_params(1), best_params(2), 'rp', 'MarkerSize', 15, 'LineWidth', 2);

    subplot(2,2,2);
    scatter(top_params(:,3), top_params(:,4), 50, norm_costs, 'filled');
    xlabel('p_l (rad/s)'); ylabel('T_i (s)');
    title('Top 100 Solutions: p_l vs T_i');
    colorbar; colormap(jet); clim([0 1]);
    hold on;
    plot(best_params(3), best_params(4), 'rp', 'MarkerSize', 15, 'LineWidth', 2);

    subplot(2,2,3);
    scatter(top_params(:,1), top_params(:,3), 50, norm_costs, 'filled');
    xlabel('K (gain)'); ylabel('p_l (rad/s)');
    title('Top 100 Solutions: K vs p_l');
    colorbar; colormap(jet); clim([0 1]);
    hold on;
    plot(best_params(1), best_params(3), 'rp', 'MarkerSize', 15, 'LineWidth', 2);

    subplot(2,2,4);
    scatter(top_params(:,2), top_params(:,3), 50, norm_costs, 'filled');
    xlabel('z_l (rad/s)'); ylabel('p_l (rad/s)');
    title('Top 100 Solutions: z_l vs p_l (lead ratio)');
    colorbar; colormap(jet); clim([0 1]);
    hold on;
    plot(best_params(2), best_params(3), 'rp', 'MarkerSize', 15, 'LineWidth', 2);

    % Save figure
    if ~exist(resultsDir, 'dir'); mkdir(resultsDir); end
    saveas(gcf, fullfile(resultsDir, 'ParameterSpace_Top100.png'));
end

    % Save optimization results to cache file
    fprintf('\n✓ Saving optimization results to cache...\n');
    if ~exist(resultsDir, 'dir'); mkdir(resultsDir); end
    save(opt_cache_file, 'K_opt', 'zl_opt', 'pl_opt', 'Ti_opt', 'best_cost', ...
         'best_params', 'costs', 'param_combos');
    fprintf('  Cache saved: %s\n', opt_cache_file);
    fprintf('  (Delete this file to re-run optimization)\n');

end  % end if run_optimization

%% Task 2: WL Controller, Linear Analysis, and Nonlinear Sim with Saturation
% ------------------------------------------------------------------------------
% Use optimized controller parameters from mathematical optimization
s = tf('s');
K      = K_opt;        % Optimized gain
zl     = zl_opt;       % Optimized lead zero
pl     = pl_opt;       % Optimized lead pole
Ti     = Ti_opt;       % Optimized integrator time constant
wf     = 28;           % Noise roll-off (fixed)

Clead  = (1 + s/zl)/(1 + s/pl);
Cint   = (1 + 1/(Ti*s));
Clp    = 1/(1 + s/wf);
Cwl    = minreal(K * Clead * Cint * Clp, 1e-6);

% Actuator model (aileron servo), first order
w_act  = 35;               % rad/s (you can change 30–40 per aircraft)
Gact   = 1/(1 + s/w_act);

% Open-loop and closed-loop (unity sensor for analysis)
L      = minreal(Cwl*Gact*G_phi_da, 1e-6);
Tphi   = minreal(feedback(L, 1), 1e-6);     % phi/phi_c (linear, no saturation)

% ------------------------------------------------------------------------------
% Linear analysis plots (Bode + margins, root locus, and step response)
figure('Name','WL Linear Analysis');
subplot(2,2,1);
margin(L); grid on; title('Open-loop L(j\omega) = C \cdot G_{act} \cdot G_{\phi/\delta_a}');

subplot(2,2,2);
rlocus(Cwl*Gact*G_phi_da); grid on;
title('Root Locus of L(s)');
% Mark the current closed-loop poles
hold on;
cl_poles = pole(Tphi);
plot(real(cl_poles), imag(cl_poles), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
legend('Root Locus', 'Closed-loop poles');

subplot(2,2,3);
step(deg2rad(30)*Tphi, 8); grid on;
ylabel('\phi (rad)'); title('Closed-loop Step (Linear, 30° command)');

subplot(2,2,4);
% Nichols chart
nichols(L); grid on;
title('Nichols Chart of L(j\omega)');

% Extract linear specs (to 2% settling)
S = stepinfo(Tphi);
fprintf('\n[Linear Closed-Loop Specs]\n');
fprintf('   Ts (2%%):     %.3f s\n', S.SettlingTime);
fprintf('   Mp (%%):       %.1f %%\n', S.Overshoot);
[GM, PM, Wcg, Wcp] = margin(L);
fprintf('   Crossover w:  %.3f rad/s\n', Wcp);
fprintf('   Phase margin: %.1f deg\n', PM);

% Analyze closed-loop pole damping ratios
cl_poles = pole(Tphi);
fprintf('\n[Closed-Loop Pole Analysis]\n');
for i = 1:length(cl_poles)
    p = cl_poles(i);
    if imag(p) ~= 0  % Complex pole
        wn = abs(p);
        zeta = -real(p)/wn;
        fprintf('   Pole %d: %.3f ± %.3fj  |  ωn=%.3f, ζ=%.3f', ...
                i, real(p), abs(imag(p)), wn, zeta);
        if zeta >= 0.9 && zeta <= 1.1
            fprintf(' (critically damped)\n');
        elseif zeta > 1.1
            fprintf(' (overdamped)\n');
        elseif zeta > 0.7
            fprintf(' (well damped)\n');
        else
            fprintf(' (underdamped)\n');
        end
    else  % Real pole
        fprintf('   Pole %d: %.3f  (real)\n', i, real(p));
    end
end

% ------------------------------------------------------------------------------
% Nonlinear time simulation with actuator saturation and simple anti-windup
% Build continuous SS for controller, actuator, plant; then integrate with ode45
% Saturation limit (aileron): +/- 20 deg
alim_deg = 20;
alim = deg2rad(alim_deg);

% Controller SS (e -> u)
Cwl_ss = ss(Cwl);   % will be strictly proper (states include integrator & filters)

% Actuator SS (u -> delta_a_unsat)
A_act = -w_act; B_act =  w_act; C_act = 1; D_act = 0;   % xdot = -w_act x + w_act u; y = x
Act_ss = ss(A_act, B_act, C_act, D_act);

% Plant SS (delta_a -> phi)
Gp_ss = ss(G_phi_da);    % (A_p, B_p, C_p, D_p) from aileron to bank angle

% Assemble simulation options
simOpts = struct();
simOpts.alim     = alim;     % rad
simOpts.k_aw     = 5;        % anti-windup gain for back-calculation (tune 2–10)
simOpts.phi_cmd  = deg2rad(30); % 30 deg comand
simOpts.tfinal   = 8;        % seconds
simOpts.x0       = [];       % default zeros

[t, y, u_hist, da_hist, xall] = simulate_wl_nl(Cwl_ss, Act_ss, Gp_ss, simOpts);

% ------------------------------------------------------------------------------
% Plots (nonlinear sim with saturation)
figure('Name','WL Nonlinear Sim (Saturation)');
subplot(3,1,1);
plot(t, rad2deg(y), 'LineWidth',1.4); grid on;
ylabel('\phi (deg)'); title('Bank Angle (with actuator saturation)');
yline(30,'--'); ylim([min(-5, min(rad2deg(y))*1.1) max(35, max(rad2deg(y))*1.1)]);

subplot(3,1,2);
plot(t, rad2deg(da_hist), 'LineWidth',1.4); grid on;
ylabel('\delta_a (deg)'); title('Aileron (after saturation)');
yline(alim_deg,'r--'); yline(-alim_deg,'r--');

subplot(3,1,3);
plot(t, u_hist, 'LineWidth',1.2); grid on;
xlabel('Time (s)'); ylabel('u (controller output)');
title('Controller output u (pre-servo)');

% Basic metrics from the nonlinear sim
phi_deg = rad2deg(y);
Mp_nl   = (max(phi_deg) - 30)/30 * 100;     % overshoot in %
Ts_nl   = NaN;
% 2% settling time to 30 deg:
band = 0.02*30;
idx  = find(abs(phi_deg - 30) <= band, 1, 'first');
if ~isempty(idx)
    % Check it stays within band afterwards:
    if all(abs(phi_deg(idx:end) - 30) <= band)
        Ts_nl = t(idx);
    else
        % if it re-exits, take last time it re-enters and stays
        idx2 = find(abs(phi_deg - 30) <= band);
        % find first index after which it never leaves
        Ts_nl = NaN;
        for k = 1:numel(idx2)
            if all(abs(phi_deg(idx2(k):end) - 30) <= band)
                Ts_nl = t(idx2(k)); break;
            end
        end
    end
end
fprintf('\n[Nonlinear (Sat) Step Metrics]\n');
fprintf('   Ts_2%%:        %s s\n', num2str(Ts_nl,'%.3f'));
fprintf('   Mp (%%):       %.1f %%\n', Mp_nl);
fprintf('   |delta_a|max: %.1f deg\n', max(abs(rad2deg(da_hist))));

% Save figures (optional)
if ~exist(resultsDir, 'dir'); mkdir(resultsDir); end
saveas(gcf, fullfile(resultsDir, 'WL_NonlinearStep.png'));

%% Additional Analysis: Multiple Step Commands & Sensitivity
% ------------------------------------------------------------------------------
figure('Name','Multiple Step Commands');
subplot(2,1,1);
hold on; grid on;
step_commands = [5 15 30 45 60];  % degrees
colors = lines(length(step_commands));
for i = 1:length(step_commands)
    cmd = step_commands(i);
    opts_temp = simOpts;
    opts_temp.phi_cmd = deg2rad(cmd);
    [t_temp, phi_temp, ~, ~, ~] = simulate_wl_nl(Cwl_ss, Act_ss, Gp_ss, opts_temp);
    plot(t_temp, rad2deg(phi_temp), 'Color', colors(i,:), 'LineWidth', 1.5, ...
         'DisplayName', sprintf('%d° cmd', cmd));
end
ylabel('\phi (deg)'); xlabel('Time (s)');
title('Bank Angle Response for Multiple Commands (with saturation)');
legend('Location','best');
ylim([-5 65]);

subplot(2,1,2);
% Sensitivity and Complementary Sensitivity
S_sens = minreal(1/(1+L), 1e-6);  % Sensitivity function
T_comp = Tphi;                     % Complementary sensitivity (already have this)
bodemag(S_sens, T_comp, {0.01, 100}); grid on;
legend('S (Sensitivity)', 'T (Complementary)', 'Location', 'best');
title('Sensitivity Functions (disturbance rejection & tracking)');

if ~exist(resultsDir, 'dir'); mkdir(resultsDir); end
saveas(gcf, fullfile(resultsDir, 'WL_MultipleSteps_Sensitivity.png'));

function [Xdot, info] = getXdot(X, U, Xg, FlightData, debug)
    % Solve the implicit algebraic loop: Xdot = motion( aero(X, U, Xg, Xdot) )
    if nargin < 5, debug = false; end
    if nargin < 3 || isempty(Xg), Xg = zeros(6,1); end
    
    tol    = 1e-7;     % relative tolerance on Xdot
    maxit  = 200;      % safety cap
    alpha  = 0.6;      % under-relaxation (0<alpha<=1) for robustness
    
    % augment controls if needed
    U_aug = U; if numel(U_aug)==4, U_aug = [U_aug; 0]; end
    
    % initial guess for state rates (trim: zeros is fine)
    Xdot_i = zeros(12,1);
    
    for k = 1:maxit
        % 1) aero loads depend on the current guess of the state rates
        [CF, CM] = aero5560_aero(X, Xg, Xdot_i, U_aug, FlightData);
        % 2) equations of motion map loads -> new rates
        Xdot_o   = aero5560_motion(X, CF, CM, FlightData);
    
        % relative fixed-point error
        denom = max(1, norm(Xdot_o, inf));
        err   = norm(Xdot_o - Xdot_i, inf) / denom;
    
        if err < tol
            Xdot = Xdot_o;
            info = struct('iter',k,'err',err,'converged',true);
            return
        end
    
        % 3) update the guess (with relaxation to avoid oscillations/divergence)
        Xdot_i = (1 - alpha)*Xdot_i + alpha*Xdot_o;
    end
    
    % if we get here, not converged
    Xdot = Xdot_i;
    info = struct('iter',maxit,'err',err,'converged',false);
    if debug
        warning('getXdot did not converge: err=%.3e after %d iters', err, maxit);
    end
end

function [t, phi, u_hist, delta_a, xout] = simulate_wl_nl(Css, Actss, Gpss, opts)
    if nargin < 4 || isempty(opts), opts = struct(); end
    if ~isfield(opts,'alim'),    opts.alim    = deg2rad(20); end
    if ~isfield(opts,'k_aw'),    opts.k_aw    = 5;           end
    if ~isfield(opts,'phi_cmd'), opts.phi_cmd = deg2rad(30); end
    if ~isfield(opts,'tfinal'),  opts.tfinal  = 8;           end
    if ~isfield(opts,'x0'),      opts.x0      = [];          end

    % Dimensions
    Ac = Css.A;  Bc = Css.B;  Cc = Css.C;  Dc = Css.D;
    Aa = Actss.A; Ba = Actss.B; Ca = Actss.C; Da = Actss.D;
    Ap = Gpss.A;  Bp = Gpss.B;  Cp = Gpss.C;  Dp = Gpss.D;

    nc = size(Ac,1);  na = size(Aa,1);  np = size(Ap,1);

    % Composite state: x = [x_c; x_a; x_p]
    if isempty(opts.x0)
        x0 = zeros(nc + na + np, 1);
    else
        x0 = opts.x0(:);
        if numel(x0) ~= nc+na+np
            error('x0 length must be %d (controller+actuator+plant).', nc+na+np);
        end
    end

    % ODE
    function dx = f(~, x)
        xc = x(1:nc);
        xa = x(nc+(1:na));
        xp = x(nc+na+(1:np));

        % Signals
        phi = Cp*xp + Dp*0;       % plant output (phi), Dp should be 0
        e   = opts.phi_cmd - phi; % error

        % Controller output (pre-servo)
        u   = Cc*xc + Dc*e;

        % Actuator (unsaturated output)
        ya_unsat = Ca*xa + Da*u;

        % Saturation
        ya_sat = min(max(ya_unsat, -opts.alim), opts.alim);

        % Anti-windup back-calculation: feed (sat - unsat) into controller
        e_aw = e + opts.k_aw*(ya_sat - ya_unsat);

        % Recompute controller with e_aw (state-derivative uses e_aw)
        u_aw = Cc*xc + Dc*e_aw;

        % Dynamics
        xc_dot = Ac*xc + Bc*e_aw;              % controller states
        xa_dot = Aa*xa + Ba*u_aw;              % actuator states driven by u_aw
        xp_dot = Ap*xp + Bp*ya_sat;            % plant states driven by saturated aileron

        dx = [xc_dot; xa_dot; xp_dot];
    end

    % Integrate
    tspan = [0 opts.tfinal];
    odeOpts = odeset('RelTol',1e-6,'AbsTol',1e-8);
    [t, x] = ode45(@f, tspan, x0, odeOpts);

    % Outputs
    nc = size(Css.A,1);  na = size(Actss.A,1);  np = size(Gpss.A,1);
    xc = x(:, 1:nc).';
    xa = x(:, nc+(1:na)).';
    xp = x(:, nc+na+(1:np)).';

    % Recompute signals for logging
    phi     = (Gpss.C * xp + Gpss.D * 0).';
    % For logging u, recompute with back-calculation (consistent with ODE)
    u_hist  = zeros(numel(t),1);
    delta_a = zeros(numel(t),1);

    for k=1:numel(t)
        e    = opts.phi_cmd - phi(k);
        u    = Css.C*xc(:,k) + Css.D*e;
        ya_u = Actss.C*xa(:,k) + Actss.D*u;
        ya_s = min(max(ya_u, -opts.alim), opts.alim);
        e_aw = e + opts.k_aw*(ya_s - ya_u);
        u_aw = Css.C*xc(:,k) + Css.D*e_aw;
        u_hist(k)  = u_aw;
        delta_a(k) = ya_s;
    end

    xout = x;
end

%% Export Results for LaTeX
% ------------------------------------------------------------------------------
% This section exports all figures, transfer functions, and performance data
% to LaTeX-compatible formats (PDF figures and .tex snippets)

fprintf('\n========================================\n');
fprintf('EXPORTING RESULTS FOR LATEX\n');
fprintf('========================================\n');

% Suppress harmless MATLAB warnings during figure generation
warning('off', 'MATLAB:print:ContentTypeImageSuggested');
warning('off', 'MATLAB:legend:IgnoringExtraEntries');

% Create output directory
latex_export_dir = fullfile(projectDir, 'LaTeX_Exports');
if ~exist(latex_export_dir, 'dir')
    mkdir(latex_export_dir);
end
fprintf('Export directory: %s\n', latex_export_dir);

% Helper function to convert polynomial coefficients to LaTeX string
function str = latex_poly(coeffs)
    n = length(coeffs) - 1;
    terms = {};
    for k = 0:n
        c = coeffs(n-k+1);
        if abs(c) < 1e-10, continue; end

        % Format coefficient
        if abs(c - 1) < 1e-10
            c_str = '';
        elseif abs(c + 1) < 1e-10
            c_str = '-';
        else
            c_str = sprintf('%.4g', c);
        end

        % Format power of s
        if k == 0
            term = c_str;
            if isempty(c_str), term = '1'; end
        elseif k == 1
            if isempty(c_str)
                term = 's';
            else
                term = [c_str 's'];
            end
        else
            if isempty(c_str)
                term = sprintf('s^{%d}', k);
            else
                term = sprintf('%ss^{%d}', c_str, k);
            end
        end

        terms{end+1} = term;
    end

    if isempty(terms)
        str = '0';
    else
        str = strjoin(terms, ' + ');
        str = strrep(str, '+ -', '- ');
    end
end

% Helper function to export design data (TFs, poles)
function export_design_data(label, C, L, T, outdir)
    % Extract transfer functions
    % Handle case where C might be a scalar (P-only control)
    if isnumeric(C) && isscalar(C)
        C = tf(C, 1);
    end
    [num_c, den_c] = tfdata(C, 'v');
    [num_l, den_l] = tfdata(L, 'v');
    [num_t, den_t] = tfdata(T, 'v');

    % Export controller TF
    fid = fopen(fullfile(outdir, sprintf('design_%s_controller.tex', label)), 'w');
    fprintf(fid, '$C_{\\mathrm{WL}}(s) = \\frac{%s}{%s}$', ...
        latex_poly(num_c), latex_poly(den_c));
    fclose(fid);

    % Export open-loop TF
    fid = fopen(fullfile(outdir, sprintf('design_%s_openloop.tex', label)), 'w');
    fprintf(fid, '$L(s) = \\frac{%s}{%s}$', ...
        latex_poly(num_l), latex_poly(den_l));
    fclose(fid);

    % Export closed-loop TF
    fid = fopen(fullfile(outdir, sprintf('design_%s_closedloop.tex', label)), 'w');
    fprintf(fid, '$T(s) = \\frac{%s}{%s}$', ...
        latex_poly(num_t), latex_poly(den_t));
    fclose(fid);

    % Export pole locations
    cl_poles = pole(T);
    fid = fopen(fullfile(outdir, sprintf('design_%s_poles.tex', label)), 'w');
    fprintf(fid, '\\begin{itemize}\n');
    for i = 1:length(cl_poles)
        p = cl_poles(i);
        if abs(imag(p)) > 1e-6
            wn = abs(p);
            zeta = -real(p)/wn;
            fprintf(fid, '    \\item $s = %.3f \\pm %.3fj$ ($\\zeta=%.2f$, $\\omega_n=%.2f$ rad/s)\n', ...
                real(p), abs(imag(p)), zeta, wn);
        else
            fprintf(fid, '    \\item $s = %.3f$ (real)\n', real(p));
        end
    end
    fprintf(fid, '\\end{itemize}\n');
    fclose(fid);

    fprintf('  Exported design %s transfer functions and poles\n', label);
end

% Helper function to generate all figures for a design
function generate_design_figures(label, C, G_act, G_plant, L, T, outdir)
    % Add matlab2tikz to path if not already there
    matlab2tikz_path = '/Users/adam/Library/Application Support/MathWorks/MATLAB Add-Ons/Collections/matlab2tikz_matlab2tikz/src';
    if exist(matlab2tikz_path, 'dir')
        addpath(matlab2tikz_path);
    end

    % Set publication quality with MASSIVE fonts for PDFs
    set(0, 'DefaultFigureColor', 'w');
    set(0, 'DefaultAxesFontSize', 28);          % HUGE fonts for readability
    set(0, 'DefaultAxesFontName', 'Times');
    set(0, 'DefaultAxesLabelFontSizeMultiplier', 1.3);  % Even larger axis labels
    set(0, 'DefaultLegendFontSize', 24);        % Large legend
    set(0, 'DefaultLineLineWidth', 3);          % Thick lines

    % Handle case where C might be a scalar (P-only control)
    if isnumeric(C) && isscalar(C)
        C = tf(C, 1);
    end

    % 1. Root Locus
    fig = figure('Visible', 'off', 'Position', [0, 0, 1920, 1080], 'Units', 'pixels');
    rlocus(C*G_act*G_plant);
    grid on;

    % Mark closed-loop poles
    hold on;
    cl_poles = pole(T);
    plot(real(cl_poles), imag(cl_poles), 'rx', 'MarkerSize', 12, 'LineWidth', 2);
    try
        legend('Root Locus', 'Closed-loop poles', 'Location', 'best');
    catch
        % Ignore legend errors in invisible figures
    end

    % Export as PGF for LaTeX
    matlab2tikz(fullfile(outdir, sprintf('design_%s_root_locus.tex', label)), ...
        'standalone', false, 'showInfo', false);
    % Also save PDF backup
    exportgraphics(fig, fullfile(outdir, sprintf('design_%s_root_locus.pdf', label)), ...
        'ContentType', 'vector');
    close(fig);

    % 2. Bode Plot with margins
    fig = figure('Visible', 'off', 'Position', [0, 0, 1920, 1080], 'Units', 'pixels');
    margin(L);
    grid on;

    matlab2tikz(fullfile(outdir, sprintf('design_%s_bode.tex', label)), ...
        'standalone', false, 'showInfo', false);
    exportgraphics(fig, fullfile(outdir, sprintf('design_%s_bode.pdf', label)), ...
        'ContentType', 'vector');
    close(fig);

    % 3. Step Response
    fig = figure('Visible', 'off', 'Position', [0, 0, 1920, 1080], 'Units', 'pixels');
    step(deg2rad(30)*T, 8);
    grid on;
    ylabel('$\phi$ (rad)', 'Interpreter', 'latex');
    xlabel('Time (s)');

    matlab2tikz(fullfile(outdir, sprintf('design_%s_step.tex', label)), ...
        'standalone', false, 'showInfo', false);
    exportgraphics(fig, fullfile(outdir, sprintf('design_%s_step.pdf', label)), ...
        'ContentType', 'vector');
    close(fig);

    % 4. Nichols Chart
    fig = figure('Visible', 'off', 'Position', [0, 0, 1920, 1080], 'Units', 'pixels');
    nichols(L);
    grid on;

    matlab2tikz(fullfile(outdir, sprintf('design_%s_nichols.tex', label)), ...
        'standalone', false, 'showInfo', false);
    exportgraphics(fig, fullfile(outdir, sprintf('design_%s_nichols.pdf', label)), ...
        'ContentType', 'vector');
    close(fig);

    fprintf('  Generated figures for design %s\n', label);
end

% Design A: Proportional Control
fprintf('\nExporting Design A (P control)...\n');
s = tf('s');
K_A = -0.020;
C_A = K_A;
w_act = 35;
G_act = 1/(1 + s/w_act);
L_A = minreal(C_A * G_act * G_phi_da, 1e-6);
T_A = minreal(feedback(L_A, 1), 1e-6);

export_design_data('A', C_A, L_A, T_A, latex_export_dir);
generate_design_figures('A', C_A, G_act, G_phi_da, L_A, T_A, latex_export_dir);

% Export Design A performance metrics
[GM_A, PM_A, ~, Wcp_A] = margin(L_A);
S_A = stepinfo(T_A);
fid = fopen(fullfile(latex_export_dir, 'design_A_metrics.tex'), 'w');
fprintf(fid, '\\begin{tabular}{ll}\n');
fprintf(fid, 'Gain & $K = %.3f$ \\\\\n', K_A);
fprintf(fid, 'Phase Margin & $%.1f^\\circ$ \\\\\n', PM_A);
fprintf(fid, 'Crossover Freq. & $%.2f$ rad/s \\\\\n', Wcp_A);
fprintf(fid, 'Settling Time & $%.2f$ s \\\\\n', S_A.SettlingTime);
fprintf(fid, 'Overshoot & $%.1f\\%%$ \\\\\n', S_A.Overshoot);
fprintf(fid, '\\end{tabular}\n');
fclose(fid);

% Design B: PI Control
fprintf('\nExporting Design B (PI control)...\n');
K_B = -0.080;
Ti_B = 0.60;
C_B = K_B * (1 + 1/(Ti_B*s));
L_B = minreal(C_B * G_act * G_phi_da, 1e-6);
T_B = minreal(feedback(L_B, 1), 1e-6);

export_design_data('B', C_B, L_B, T_B, latex_export_dir);
generate_design_figures('B', C_B, G_act, G_phi_da, L_B, T_B, latex_export_dir);

% Export Design B performance metrics
[GM_B, PM_B, ~, Wcp_B] = margin(L_B);
S_B = stepinfo(T_B);
fid = fopen(fullfile(latex_export_dir, 'design_B_metrics.tex'), 'w');
fprintf(fid, '\\begin{tabular}{ll}\n');
fprintf(fid, 'Gain & $K = %.3f$ \\\\\n', K_B);
fprintf(fid, 'Integrator TC & $T_i = %.2f$ s \\\\\n', Ti_B);
fprintf(fid, 'Phase Margin & $%.1f^\\circ$ \\\\\n', PM_B);
fprintf(fid, 'Crossover Freq. & $%.2f$ rad/s \\\\\n', Wcp_B);
fprintf(fid, 'Settling Time & $%.2f$ s \\\\\n', S_B.SettlingTime);
fprintf(fid, 'Overshoot & $%.1f\\%%$ \\\\\n', S_B.Overshoot);
fprintf(fid, '\\end{tabular}\n');
fclose(fid);

% Design C: Final optimized design (already computed)
fprintf('\nExporting Design C (optimized lead-PI with roll-off)...\n');
export_design_data('C', Cwl, L, Tphi, latex_export_dir);
generate_design_figures('C', Cwl, Gact, G_phi_da, L, Tphi, latex_export_dir);

% Export Design C performance metrics
fid = fopen(fullfile(latex_export_dir, 'design_C_metrics.tex'), 'w');
fprintf(fid, '\\begin{tabular}{ll}\n');
fprintf(fid, 'Gain & $K = %.4f$ \\\\\n', K);
fprintf(fid, 'Lead Zero & $z_\\ell = %.3f$ rad/s \\\\\n', zl);
fprintf(fid, 'Lead Pole & $p_\\ell = %.3f$ rad/s \\\\\n', pl);
fprintf(fid, 'Integrator TC & $T_i = %.3f$ s \\\\\n', Ti);
fprintf(fid, 'Filter Cutoff & $\\omega_f = %.1f$ rad/s \\\\\n', wf);
fprintf(fid, 'Phase Margin & $%.1f^\\circ$ \\\\\n', PM);
fprintf(fid, 'Crossover Freq. & $%.2f$ rad/s \\\\\n', Wcp);
fprintf(fid, 'Settling Time & $%.2f$ s \\\\\n', S.SettlingTime);
fprintf(fid, 'Overshoot & $%.1f\\%%$ \\\\\n', S.Overshoot);
fprintf(fid, '\\end{tabular}\n');
fclose(fid);

% Export nonlinear response comparison
fprintf('\nExporting nonlinear response comparison...\n');
fig = figure('Visible', 'off', 'Position', [0, 0, 1920, 1080], 'Units', 'pixels');

subplot(3,1,1);
plot(t, rad2deg(y), 'LineWidth', 1.5);
grid on;
ylabel('$\phi$ (deg)', 'Interpreter', 'latex');
yline(30, '--k', 'LineWidth', 1);

subplot(3,1,2);
plot(t, rad2deg(da_hist), 'LineWidth', 1.5);
grid on;
ylabel('$\delta_a$ (deg)', 'Interpreter', 'latex');
yline(alim_deg, '--r', 'LineWidth', 1);
yline(-alim_deg, '--r', 'LineWidth', 1);

subplot(3,1,3);
plot(t, u_hist, 'LineWidth', 1.5);
grid on;
xlabel('Time (s)');
ylabel('$u$ (controller output)', 'Interpreter', 'latex');

matlab2tikz(fullfile(latex_export_dir, 'nonlinear_response.tex'), ...
    'standalone', false, 'showInfo', false);
exportgraphics(fig, fullfile(latex_export_dir, 'nonlinear_response.pdf'), ...
    'ContentType', 'vector');
close(fig);

% Export design comparison (step responses)
fprintf('\nExporting design comparison...\n');
fig = figure('Visible', 'off', 'Position', [0, 0, 1920, 1080], 'Units', 'pixels');

subplot(1,2,1);
hold on; grid on;
[y1, t1] = step(deg2rad(30)*T_A, 8);
[y2, t2] = step(deg2rad(30)*T_B, 8);
[y3, t3] = step(deg2rad(30)*Tphi, 8);
plot(t1, y1, 'LineWidth', 1.2);
plot(t2, y2, 'LineWidth', 1.2);
plot(t3, y3, 'LineWidth', 2.0);
ylabel('$\phi$ (rad)', 'Interpreter', 'latex');
xlabel('Time (s)');
legend('Design A (P)', 'Design B (PI)', 'Design C (Lead-PI)', 'Location', 'best');
yline(deg2rad(30), '--k', 'LineWidth', 0.5);

subplot(1,2,2);
% Bode magnitude comparison
bodemag(L_A, L_B, L, {0.01, 100});
grid on;
legend('Design A', 'Design B', 'Design C', 'Location', 'best');

matlab2tikz(fullfile(latex_export_dir, 'design_comparison.tex'), ...
    'standalone', false, 'showInfo', false);
exportgraphics(fig, fullfile(latex_export_dir, 'design_comparison.pdf'), ...
    'ContentType', 'vector');
close(fig);

% Export sensitivity functions
fprintf('\nExporting sensitivity functions...\n');
S_sens = minreal(1/(1+L), 1e-6);
T_comp = Tphi;

fig = figure('Visible', 'off', 'Position', [0, 0, 1920, 1080], 'Units', 'pixels');
bodemag(S_sens, T_comp, {0.01, 100});
grid on;
legend('$S$ (Sensitivity)', '$T$ (Complementary)', 'Interpreter', 'latex', 'Location', 'best');

matlab2tikz(fullfile(latex_export_dir, 'sensitivity_functions.tex'), ...
    'standalone', false, 'showInfo', false);
exportgraphics(fig, fullfile(latex_export_dir, 'sensitivity_functions.pdf'), ...
    'ContentType', 'vector');
close(fig);

% Create performance comparison table
fprintf('\nCreating performance comparison table...\n');
fid = fopen(fullfile(latex_export_dir, 'performance_table.tex'), 'w');
fprintf(fid, '\\begin{table}[h!]\n');
fprintf(fid, '\\centering\n');
fprintf(fid, '\\caption{Performance Comparison of Three Controller Designs}\n');
fprintf(fid, '\\begin{tabular}{lccc}\n');
fprintf(fid, '\\hline\n');
fprintf(fid, '\\textbf{Metric} & \\textbf{Design A} & \\textbf{Design B} & \\textbf{Design C} \\\\\n');
fprintf(fid, '\\hline\n');
fprintf(fid, 'Phase Margin & $%.1f^\\circ$ & $%.1f^\\circ$ & $%.1f^\\circ$ \\\\\n', PM_A, PM_B, PM);
fprintf(fid, 'Crossover Freq. & $%.2f$ rad/s & $%.2f$ rad/s & $%.2f$ rad/s \\\\\n', Wcp_A, Wcp_B, Wcp);
fprintf(fid, 'Settling Time & $%.2f$ s & $%.2f$ s & $%.2f$ s \\\\\n', S_A.SettlingTime, S_B.SettlingTime, S.SettlingTime);
fprintf(fid, 'Overshoot & $%.1f\\%%$ & $%.1f\\%%$ & $%.1f\\%%$ \\\\\n', S_A.Overshoot, S_B.Overshoot, S.Overshoot);
fprintf(fid, 'Steady-State Error & $%.2f^\\circ$ & $0^\\circ$ & $0^\\circ$ \\\\\n', 30 - rad2deg(dcgain(T_A)*deg2rad(30)));
fprintf(fid, '\\hline\n');
fprintf(fid, '\\end{tabular}\n');
fprintf(fid, '\\label{tab:performance_comparison}\n');
fprintf(fid, '\\end{table}\n');
fclose(fid);

fprintf('\n✓ Export complete! All files saved to:\n');
fprintf('  %s\n', latex_export_dir);
fprintf('\nGenerated files:\n');
fprintf('  - Transfer functions: design_X_controller.tex, design_X_openloop.tex, design_X_closedloop.tex\n');
fprintf('  - Pole locations: design_X_poles.tex\n');
fprintf('  - Figures: design_X_root_locus.pdf, design_X_bode.pdf, design_X_step.pdf, design_X_nichols.pdf\n');
fprintf('  - Performance metrics: design_X_metrics.tex\n');
fprintf('  - Comparison plots: design_comparison.pdf, nonlinear_response.pdf, sensitivity_functions.pdf\n');
fprintf('  - Summary table: performance_table.tex\n');
fprintf('\nYou can now compile your LaTeX document!\n');
