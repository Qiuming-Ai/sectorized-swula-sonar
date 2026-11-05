%% Cross-correlation peak to autocorrelation peak ratio (dB) vs. pulse width T
% CDM-Zadoff-Chu multibeam excitation
% Bandwidth: 100 kHz; Variable pulse width T; M = [2, 4, 6, 8] beams
clear; clc;

% -------- Parameters --------
B    = 100e3;                    % Signal bandwidth (Hz)
M_list = [2, 4, 6, 8];          % Number of simultaneously active beams
L_min = 16;                      % Minimum code length
L_max = 8192;                    % Maximum code length
M_points = 100;                  % Number of points for T axis
L_vec = round(logspace(log10(L_min), log10(L_max), M_points));  % Code length vector

% Remove duplicates (logspace might generate duplicate integers)
L_vec = unique(L_vec);

% Convert code length to pulse width: T = L / B
T_vec = L_vec / B;               % Pulse width vector (s)

% Different line styles and markers for different M
styles  = {'-', '--', '-.', ':'};
markers = {'o', 's', '^', 'd'};
mkEvery = round(M_points/20);   % Marker every N points

% -------- Compute normalized peak ratio (dB) for each M --------
% According to theory:
% Single interferer: R_{i->k}^{(1)} = 1/sqrt(L), dB = -10*log10(L)
% M beams worst-case total: R_k^(tot) <= (M-1)/sqrt(L)
% dB: 20*log10(M-1) - 10*log10(L)

R_dB = zeros(length(M_list), length(L_vec));

for iM = 1:length(M_list)
    M = M_list(iM);
    % Worst-case normalized cross-correlation peak ratio
    % R = (M-1) / sqrt(L)
    % dB = 20*log10(M-1) - 10*log10(L)
    R_dB(iM, :) = 20*log10(M-1) - 10*log10(L_vec);
end

% Note: Theoretical bound is used directly.
% According to ZC sequence properties, for prime length L and different roots:
% - Autocorrelation peak: |r_kk[0]| = L
% - Cross-correlation magnitude: |r_ik[ell]| = sqrt(L) (constant for all ell)
% - Worst-case M-beam total interference: (M-1)*sqrt(L)
% - Normalized ratio: (M-1)/sqrt(L), dB = 20*log10(M-1) - 10*log10(L)

% -------- Plot (dB) --------
figure('Color','w', 'Position', [100, 100, 900, 600]); hold on;

% Plot theoretical curves
for iM = 1:length(M_list)
    M = M_list(iM);
    plot(T_vec*1e3, R_dB(iM,:), styles{iM}, 'LineWidth', 1.8, ...
        'Marker', markers{iM}, 'MarkerIndices', 1:mkEvery:length(T_vec), ...
        'MarkerSize', 6, 'DisplayName', sprintf('M = %d beams', M));
end

grid on; box on;
xlabel('Pulse width T (ms)', 'FontSize', 12);
ylabel('Peak amplitude ratio (dB)', 'FontSize', 12);
title(sprintf('CDM-Zadoff-Chu: Cross-correlation peak ratio vs. pulse width (B = %.0f kHz)', B/1e3), ...
    'FontSize', 13, 'FontWeight', 'bold');

% Set x-axis to linear scale (default)
% set(gca, 'XScale', 'log');  % Use log scale if needed

% Set y-axis limits
ylim([-60, 5]);
xlim([T_vec(1)*1e3, T_vec(end)*1e3]);

legend('Location', 'southwest', 'FontSize', 10);

% Add grid
set(gca, 'GridAlpha', 0.3, 'MinorGridAlpha', 0.1);
set(gca, 'XMinorGrid', 'on', 'YMinorGrid', 'on');

% -------- Optional: Export results --------
Result.L_vec = L_vec;
Result.T_vec = T_vec;
Result.M_list = M_list;
Result.R_dB = R_dB;
Result.B = B;
assignin('base', 'Result', Result);

fprintf('\n========== Summary ==========\n');
fprintf('Bandwidth: %.0f kHz\n', B/1e3);
fprintf('Code length range: %d - %d\n', L_min, L_max);
fprintf('Pulse width range: %.3f - %.3f ms\n', T_vec(1)*1e3, T_vec(end)*1e3);
fprintf('Number of beams tested: ');
fprintf('%d ', M_list);
fprintf('\n===========================\n');


