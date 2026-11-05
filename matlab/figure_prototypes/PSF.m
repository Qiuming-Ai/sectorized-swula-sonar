% Publication package: read the archived normalized slices, or the optional
% original full workspace when it is present next to this script.
here = fileparts(mfilename('fullpath'));
matFile = fullfile(here,'PSF.mat');
rows = [851 791 911];
if isfile(matFile)
    D = load(matFile,'PSF_Rx','AF_Tx','G2','grid');
    angle = D.grid.tgt_deg(:);
    rx = (D.PSF_Rx(rows,:)/max(D.PSF_Rx(:))).';
    tx = (D.AF_Tx(1,:)/max(D.AF_Tx(:))).';
    product = (D.G2(rows,:)/max(D.G2(:))).';
else
    C = readtable(fullfile(here,'..','..','data','design','figure2_legacy_psf_slices.csv'));
    angle = C.angle_deg;
    rx = C{:,{'rx_focus_0','rx_focus_minus6','rx_focus_plus6'}};
    tx = C.tx_focus_0;
    product = C{:,{'product_focus_0','product_focus_minus6','product_focus_plus6'}};
end
% Keep the archived global normalization; do not renormalize individual traces.
Rx_PSF = 20*log10(rx);
Tx_PSF = 20*log10(tx);
G2_PSF = 20*log10(product);
styles = {'-','--',':'};

% Create figure with proper settings for publication
figure('Position', [100, 100, 800, 600], 'Color', 'white');

% Subplot 1: Tx&Rx PSF Slice
subplot(2,1,1);
hold on;
for k = 1:3
    plot(angle, Rx_PSF(:,k), 'LineWidth', 2, ...
         'DisplayName', sprintf('Rx Focus=%.1f°',angle(rows(k))), 'LineStyle', styles{k});
end
plot(angle, Tx_PSF, 'LineWidth', 2, 'DisplayName', 'Tx Focus=0.0°', 'LineStyle', '-.');
hold off;
title('Tx \& Rx PSF Slice', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Angle \theta (deg)', 'FontSize', 11);
ylabel('Normalized Amplitude (dB)', 'FontSize', 11);
legend('Location', 'best', 'FontSize', 10, 'Box', 'on');
grid on;
xlim([-40, 40]);
ylim([-60, 0]);
set(gca, 'FontSize', 10, 'TickDir', 'out', 'Box', 'on');

% Subplot 2: Two-way PSF Slice
subplot(2,1,2);
hold on;
for k = 1:3
    plot(angle, G2_PSF(:,k), 'LineWidth', 2, ...
         'DisplayName', sprintf('Rx Focus=%.1f°',angle(rows(k))), 'LineStyle', styles{k});
end
hold off;
title('Two-way PSF Slice', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Angle \theta (deg)', 'FontSize', 11);
ylabel('Normalized Amplitude (dB)', 'FontSize', 11);
legend('Location', 'best', 'FontSize', 10, 'Box', 'on');
grid on;
xlim([-40, 40]);
ylim([-60, 0]);
set(gca, 'FontSize', 10, 'TickDir', 'out', 'Box', 'on');

