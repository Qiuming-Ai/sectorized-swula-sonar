function metrics = compute_image_metrics(truth, estimate, guard_radius)
%COMPUTE_IMAGE_METRICS Toolbox-free image fidelity and leakage metrics.
if nargin < 3
    guard_radius = [2,1];
end
t = abs(truth);
p = abs(estimate);
t = t / (max(t(:)) + eps);

% Nonnegative least-squares amplitude alignment avoids rewarding gain alone.
alpha = max(0, sum(t(:).*p(:)) / (sum(p(:).^2) + eps));
p = alpha*p;
nmse = sum((p(:)-t(:)).^2) / (sum(t(:).^2) + eps);

support = t > max(1e-8, 0.01*max(t(:)));
kernel = ones(2*guard_radius(1)+1, 2*guard_radius(2)+1);
guard = conv2(double(support), kernel, 'same') > 0;
inside_energy = sum(p(guard).^2);
outside_energy = sum(p(~guard).^2);
support_leakage = outside_energy / (inside_energy + outside_energy + eps);

if any(guard(:)), main_peak = max(p(guard)); else, main_peak = max(p(:)); end
if any(~guard(:)), ghost_peak = max(p(~guard)); else, ghost_peak = 0; end
ghost_to_main_db = 20*log10((ghost_peak+eps)/(main_peak+eps));

background = p(~guard);
target = p(guard);
if isempty(background), background = 0; end
if isempty(target), target = 0; end
cnr = (mean(target)-mean(background)) / (std(background)+eps);

% Global SSIM form. It is intentionally toolbox-free and labeled G-SSIM.
c1 = 0.01^2;
c2 = 0.03^2;
mu_t = mean(t(:)); mu_p = mean(p(:));
var_t = var(t(:),1); var_p = var(p(:),1);
cov_tp = mean((t(:)-mu_t).*(p(:)-mu_p));
gssim = ((2*mu_t*mu_p+c1)*(2*cov_tp+c2)) / ...
    ((mu_t^2+mu_p^2+c1)*(var_t+var_p+c2));

metrics.NMSE = nmse;
metrics.NMSEdB = 10*log10(nmse+eps);
metrics.GlobalSSIM = gssim;
metrics.SupportLeakage = support_leakage;
metrics.GhostToMainDb = ghost_to_main_db;
metrics.CNR = cnr;
metrics.ScaleFactor = alpha;
end
