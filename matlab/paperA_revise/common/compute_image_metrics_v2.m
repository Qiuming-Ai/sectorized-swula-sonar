function metrics = compute_image_metrics_v2(truth,estimate,truth_support,clean_reference,spec)
%COMPUTE_IMAGE_METRICS_V2 Unified image metrics with explicit mask semantics.
% Truth fidelity is evaluated on magnitude after a nonnegative scalar fit.
% Noise NMSE compares the complex noisy output with the same method's
% noise-free output after a single complex scalar fit.
arguments
    truth
    estimate
    truth_support logical
    clean_reference = []
    spec.CoreRadius (1,2) double = [2,1]
    spec.GuardRadius (1,2) double = [5,2]
end
if ~isequal(size(truth),size(estimate),size(truth_support))
    error('truth, estimate, and truth_support must have identical sizes.');
end

t = abs(truth);
p = abs(estimate);
t = t/(max(t(:))+eps);
alpha_truth = max(0,sum(t(:).*p(:))/(sum(p(:).^2)+eps));
p_aligned = alpha_truth*p;
truth_nmse = sum((p_aligned(:)-t(:)).^2)/(sum(t(:).^2)+eps);

core = dilate_mask(truth_support,spec.CoreRadius);
guard = dilate_mask(truth_support,spec.GuardRadius);
background = ~guard;
if ~any(core(:)) || ~any(background(:))
    cnr = NaN;
else
    cnr = (mean(p(core))-mean(p(background)))/(std(p(background),1)+eps);
end
inside_energy = sum(abs(estimate(guard)).^2);
outside_energy = sum(abs(estimate(background)).^2);
leakage = outside_energy/(inside_energy+outside_energy+eps);
if any(background(:))
    max_outside_db = 20*log10((max(p(background))+eps)/(max(p(:))+eps));
else
    max_outside_db = NaN;
end

truth_corr = normalized_correlation(t,p_aligned);
gssim = global_ssim(t,p_aligned);

noise_nmse = NaN;
noise_scale = NaN;
if ~isempty(clean_reference)
    if ~isequal(size(clean_reference),size(estimate))
        error('clean_reference and estimate must have identical sizes.');
    end
    denominator = sum(abs(estimate(:)).^2)+eps;
    noise_scale = sum(conj(estimate(:)).*clean_reference(:))/denominator;
    residual = noise_scale*estimate-clean_reference;
    noise_nmse = sum(abs(residual(:)).^2)/(sum(abs(clean_reference(:)).^2)+eps);
end

metrics.TruthNMSEdB = 10*log10(truth_nmse+eps);
metrics.NoiseNMSEdB = 10*log10(noise_nmse+eps);
metrics.TruthCorrelation2D = truth_corr;
metrics.GlobalSSIM = gssim;
metrics.CNR = cnr;
metrics.SupportLeakage = leakage;
metrics.MaxOutsideDb = max_outside_db;
metrics.TruthScaleFactor = alpha_truth;
metrics.NoiseScaleFactor = noise_scale;
metrics.CorePixelCount = nnz(core);
metrics.GuardPixelCount = nnz(guard);
metrics.BackgroundPixelCount = nnz(background);
end

function expanded = dilate_mask(mask,radius)
kernel = ones(2*radius(1)+1,2*radius(2)+1);
expanded = conv2(double(mask),kernel,'same')>0;
end

function value = normalized_correlation(a,b)
a=a(:)-mean(a(:)); b=b(:)-mean(b(:));
value=sum(a.*b)/sqrt((sum(a.^2)+eps)*(sum(b.^2)+eps));
end

function value = global_ssim(a,b)
c1=0.01^2; c2=0.03^2;
mu_a=mean(a(:)); mu_b=mean(b(:));
var_a=var(a(:),1); var_b=var(b(:),1);
cov_ab=mean((a(:)-mu_a).*(b(:)-mu_b));
value=((2*mu_a*mu_b+c1)*(2*cov_ab+c2))/ ...
    ((mu_a^2+mu_b^2+c1)*(var_a+var_b+c2));
end
