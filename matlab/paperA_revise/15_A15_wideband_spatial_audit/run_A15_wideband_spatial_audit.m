function metrics=run_A15_wideband_spatial_audit(mode)
%RUN_A15_WIDEBAND_SPATIAL_AUDIT Quantify the centre-frequency approximation.
if nargin<1,mode="quick";end
this_dir=fileparts(mfilename('fullpath'));root=fileparts(this_dir);addpath(fullfile(root,'common'));
cfg=config_A15(mode);out_dir=fullfile(this_dir,'results');if ~exist(out_dir,'dir'),mkdir(out_dir);end
angles=cfg.angle_grid_deg(:).';items=repmat(struct(),0,1);profiles=repmat(struct(),0,1);
for im=1:numel(cfg.methods)
    method=cfg.methods(im);[waveforms,meta]=make_waveform_bank(cfg,method);
    weights=waveform_power_weights(waveforms,cfg.waveform_fs,cfg.frequency_grid_hz-cfg.fc);
    for it=1:numel(cfg.target_angles_deg)
        target=cfg.target_angles_deg(it);sector=it;
        frequency0=meta.center_frequency_hz(sector);
        approx=spatial_cut(cfg,frequency0,angles,target,sector);
        reference=complex(zeros(size(angles)));
        for jf=1:numel(cfg.frequency_grid_hz)
            frequency=cfg.frequency_grid_hz(jf);
            reference=reference+weights(jf,sector).*spatial_cut(cfg,frequency,angles,target,sector);
        end
        masks=sector_masks(cfg,angles);
        approx(~masks(sector,:))=0;
        reference(~masks(sector,:))=0;
        approx=abs(approx);reference=abs(reference);approx=approx/(max(approx)+eps);reference=reference/(max(reference)+eps);
        qa=cut_metrics(approx,angles,target);qr=cut_metrics(reference,angles,target);
        aligned=max(0,sum(approx.*reference)/(sum(approx.^2)+eps))*approx;
        nmse=sum((aligned-reference).^2)/(sum(reference.^2)+eps);
        corr=normalized_corr(approx,reference);
        item.Method=method;item.Sector=sector;item.TargetAngleDeg=target;item.CentreFrequencyKHz=frequency0/1e3;
        item.ApproxPeakAngleDeg=qa.PeakAngleDeg;item.ReferencePeakAngleDeg=qr.PeakAngleDeg;item.PeakAngleDifferenceDeg=qa.PeakAngleDeg-qr.PeakAngleDeg;
        item.ApproxHPBWDeg=qa.HPBWDeg;item.ReferenceHPBWDeg=qr.HPBWDeg;item.HPBWDifferenceDeg=qa.HPBWDeg-qr.HPBWDeg;
        item.ApproxPSLRDb=qa.PSLRDb;item.ReferencePSLRDb=qr.PSLRDb;item.PSLRDifferenceDb=qa.PSLRDb-qr.PSLRDb;
        item.Correlation=corr;item.NMSEdB=10*log10(nmse+eps);items=append(items,item);
        profiles=append_profile(profiles,method,sector,target,angles,approx,reference);
    end
end

function cut=spatial_cut(cfg,frequency,scan_angles,target_angle,sector)
% Compute only one target-angle cut; avoid forming an Nangle-by-Nangle matrix.
lambda_c=cfg.c/cfg.fc;k=2*pi*frequency/cfg.c;
ptx=((0:cfg.n_tx-1)-(cfg.n_tx-1)/2).'*cfg.tx_spacing_lambda*lambda_c;
prx=((0:cfg.n_rx_sparse-1)-(cfg.n_rx_sparse-1)/2).'*cfg.rx_spacing_lambda*lambda_c;
wtx=local_window("hamming",cfg.n_tx);wtx=wtx/sum(wtx);
wrx=local_window("hamming",cfg.n_rx_sparse);wrx=wrx/sum(wrx);
u0=sind(cfg.focus_angles_deg(sector));ut=sind(target_angle);us=sind(scan_angles);
tx_scalar=sum(wtx.*exp(1j*k*ptx*(ut-u0)));
rx_cut=(wrx.'*exp(1j*k*prx*(ut-us)));
tx_width=cfg.tx_element_width_lambda*lambda_c;
rx_width=cfg.rx_element_width_lambda*lambda_c;
tx_scalar=tx_scalar*normalized_sinc_scalar(tx_width*frequency/cfg.c*ut);
rx_cut=rx_cut*normalized_sinc_scalar(rx_width*frequency/cfg.c*ut);
cut=tx_scalar.*rx_cut;
end

function y=normalized_sinc_scalar(x)
if abs(x)<=eps,y=1;else,y=sin(pi*x)/(pi*x);end
end
metrics=struct2table(items);profile_table=struct2table(profiles);
writetable(metrics,fullfile(out_dir,'A15_wideband_metrics.csv'));writetable(profile_table,fullfile(out_dir,'A15_wideband_profiles.csv'));
parameter_snapshot(cfg,fullfile(out_dir,'A15_parameters.json'));save(fullfile(out_dir,'A15_workspace.mat'),'cfg','metrics','profile_table');disp(metrics);
end

function weights=waveform_power_weights(waveforms,fs,offset_grid)
nfft=2^nextpow2(16*size(waveforms,1));f=(-nfft/2:nfft/2-1).'*fs/nfft;weights=zeros(numel(offset_grid),size(waveforms,2));
for k=1:size(waveforms,2)
 spectrum=abs(fftshift(fft(waveforms(:,k),nfft))).^2;weights(:,k)=interp1(f,spectrum,offset_grid(:),'linear',0);weights(:,k)=weights(:,k)/(sum(weights(:,k))+eps);
end
end
function q=cut_metrics(curve,angles,target)
candidate=find(abs(angles-target)<=3);[~,p]=max(curve(candidate));peak=candidate(p);[lo,hi]=first_minima(curve,peak);main=false(size(curve));main(lo:hi)=true;
q.PeakAngleDeg=angles(peak);q.HPBWDeg=threshold_width(curve,angles,peak,1/sqrt(2));q.PSLRDb=20*log10(max(curve(~main))+eps);
end
function [lo,hi]=first_minima(x,peak),lo=1;hi=numel(x);for k=peak-1:-1:2,if x(k)<=x(k-1)&&x(k)<=x(k+1),lo=k;break;end,end;for k=peak+1:numel(x)-1,if x(k)<=x(k-1)&&x(k)<=x(k+1),hi=k;break;end,end,end
function width=threshold_width(x,a,p,t),l=find(x(1:p)<=t,1,'last');r0=find(x(p:end)<=t,1,'first');if isempty(l),l=1;end;if isempty(r0),r=numel(x);else,r=p+r0-1;end;width=a(r)-a(l);end
function value=normalized_corr(a,b),a=a(:)-mean(a);b=b(:)-mean(b);value=sum(a.*b)/sqrt((sum(a.^2)+eps)*(sum(b.^2)+eps));end
function rows=append(rows,item),if isempty(rows),rows=item;else,rows(end+1,1)=item;end,end
function rows=append_profile(rows,method,sector,target,angles,approx,reference)
for k=1:numel(angles)
 item.Method=method;item.Sector=sector;item.TargetAngleDeg=target;item.ScanAngleDeg=angles(k);item.CentreFrequencyMagnitudeDb=max(-100,20*log10(approx(k)+eps));item.WaveformWeightedMagnitudeDb=max(-100,20*log10(reference(k)+eps));rows=append(rows,item);
end
end
