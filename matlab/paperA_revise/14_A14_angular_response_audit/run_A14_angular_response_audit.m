function metrics=run_A14_angular_response_audit(mode)
%RUN_A14_ANGULAR_RESPONSE_AUDIT Separate receive geometry from effective PSF.
if nargin<1,mode="quick";end
this_dir=fileparts(mfilename('fullpath'));root=fileparts(this_dir);addpath(fullfile(root,'common'));
cfg=config_A14(mode);out_dir=fullfile(this_dir,'results');if ~exist(out_dir,'dir'),mkdir(out_dir);end
angles=cfg.angle_grid_deg(:).';lambda_c=cfg.c/cfg.fc;

% Pure receive responses: CDM/FDM share exactly the same SW-ULA geometry.
receive_names=["dense","sw_ula"];
receive_labels=["Dense ULA","SW-ULA (common to CDM and FDM)"];
receive_curves=zeros(numel(angles),numel(receive_names),numel(cfg.target_angles_deg));
effective_curves=zeros(numel(angles),numel(cfg.methods),numel(cfg.target_angles_deg));
items=repmat(struct(),0,1);profiles=repmat(struct(),0,1);

for it=1:numel(cfg.target_angles_deg)
    target=cfg.target_angles_deg(it);
    for ig=1:numel(receive_names)
        if receive_names(ig)=="dense",nrx=cfg.n_rx_dense;spacing=cfg.dense_rx_spacing_lambda;else,nrx=cfg.n_rx_sparse;spacing=cfg.rx_spacing_lambda;end
        positions=((0:nrx-1)-(nrx-1)/2).'*spacing*lambda_c;
        taper=local_window("hamming",nrx);taper=taper/sum(taper);
        curve=abs((taper.*exp(-1j*2*pi*cfg.fc/cfg.c*positions*sind(target))).'* ...
            exp(1j*2*pi*cfg.fc/cfg.c*positions*sind(angles)));
        curve=curve/(max(curve)+eps);receive_curves(:,ig,it)=curve(:);
        q=cut_metrics(curve,angles,target);
        item.ResponseType="pure_receive";item.Method=receive_names(ig);item.MethodLabel=receive_labels(ig);item.TargetAngleDeg=target;
        item.HPBWDeg=q.HPBWDeg;item.FNBWDeg=q.FNBWDeg;item.PSLRDb=q.PSLRDb;item.ISLRDb=q.ISLRDb;items=append(items,item);
        profiles=append_profile(profiles,"pure_receive",receive_names(ig),target,angles,curve);
    end

    for im=1:numel(cfg.methods)
        method=cfg.methods(im);
        if method=="dense"
            curve=receive_curves(:,1,it).';
        else
            [~,meta]=make_waveform_bank(cfg,method);sector=find_sector(cfg,target);
            frequency=meta.center_frequency_hz(sector);
            [tx,rx]=array_responses(cfg,frequency,angles);
            [~,target_idx]=min(abs(angles-target));
            curve=abs(tx(sector,target_idx).*rx(:,target_idx).');
            mask=sector_masks(cfg,angles);curve(~mask(sector,:))=0;
            curve=curve/(max(curve)+eps);
        end
        effective_curves(:,im,it)=curve(:);q=cut_metrics(curve,angles,target);
        item.ResponseType="effective_angular_psf";item.Method=method;item.MethodLabel=cfg.method_labels(im);item.TargetAngleDeg=target;
        item.HPBWDeg=q.HPBWDeg;item.FNBWDeg=q.FNBWDeg;item.PSLRDb=q.PSLRDb;item.ISLRDb=q.ISLRDb;items=append(items,item);
        profiles=append_profile(profiles,"effective_angular_psf",method,target,angles,curve);
    end
end
metrics=struct2table(items);profile_table=struct2table(profiles);
writetable(metrics,fullfile(out_dir,'A14_angular_metrics.csv'));writetable(profile_table,fullfile(out_dir,'A14_angular_profiles.csv'));
parameter_snapshot(cfg,fullfile(out_dir,'A14_parameters.json'));
save(fullfile(out_dir,'A14_workspace.mat'),'cfg','metrics','profile_table','receive_curves','effective_curves');disp(metrics);
end

function sector=find_sector(cfg,angle)
masks=sector_masks(cfg,angle);
sector=find(masks(:,1),1,'first');
end
function q=cut_metrics(curve,angles,target)
curve=abs(curve(:).');curve=curve/(max(curve)+eps);
% A super-wavelength receive ULA has equal-height grating replicas. Select
% the local lobe associated with the declared target rather than whichever
% numerically equal replica happens to be returned by a global max.
candidate=find(abs(angles-target)<=3);
[~,local_peak]=max(curve(candidate));peak=candidate(local_peak);
[lo,hi]=first_minima(curve,peak);main=false(size(curve));main(lo:hi)=true;
q.HPBWDeg=threshold_width(curve,angles,peak,1/sqrt(2));q.FNBWDeg=angles(hi)-angles(lo);
q.PSLRDb=20*log10(max(curve(~main))+eps);q.ISLRDb=10*log10((sum(curve(~main).^2)+eps)/(sum(curve(main).^2)+eps));
if abs(angles(peak)-target)>0.25,warning('Local peak shifted by %.3f deg.',angles(peak)-target);end
end
function [lo,hi]=first_minima(x,peak)
lo=1;hi=numel(x);for k=peak-1:-1:2,if x(k)<=x(k-1)&&x(k)<=x(k+1),lo=k;break;end,end
for k=peak+1:numel(x)-1,if x(k)<=x(k-1)&&x(k)<=x(k+1),hi=k;break;end,end
end
function width=threshold_width(x,axis_values,peak,threshold)
l=find(x(1:peak)<=threshold,1,'last');r0=find(x(peak:end)<=threshold,1,'first');if isempty(l),l=1;end;if isempty(r0),r=numel(x);else,r=peak+r0-1;end;width=axis_values(r)-axis_values(l);
end
function rows=append(rows,item),if isempty(rows),rows=item;else,rows(end+1,1)=item;end,end
function rows=append_profile(rows,type,method,target,angles,curve)
for k=1:numel(angles)
 item.ResponseType=type;item.Method=method;item.TargetAngleDeg=target;item.ScanAngleDeg=angles(k);item.Magnitude=curve(k);item.MagnitudeDb=max(-100,20*log10(curve(k)+eps));
 rows=append(rows,item);
end
end
