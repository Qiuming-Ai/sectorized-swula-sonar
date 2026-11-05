function summary=run_A19_sector_count_optimization(mode)
%RUN_A19_SECTOR_COUNT_OPTIMIZATION Minimum feasible broadband sector count.
if nargin<1,mode="quick";end
this_dir=fileparts(mfilename('fullpath'));
root=fileparts(this_dir);
addpath(fullfile(root,'common'));
cfg=config_A19(mode);
out_dir=fullfile(this_dir,'results');
if ~exist(out_dir,'dir'),mkdir(out_dir);end

items=[];
for inb=1:numel(cfg.candidate_sector_counts)
    nb=cfg.candidate_sector_counts(inb);
    edges=linspace(cfg.fov_deg(1),cfg.fov_deg(2),nb+1);
    centers=zeros(1,nb);
    for k=1:nb
        candidates=edges(k):cfg.design_angle_step_deg:edges(k+1);
        score=-inf(size(candidates));
        for ic=1:numel(candidates)
            score(ic)=sector_worst_coverage(cfg,candidates(ic),edges(k:k+1));
        end
        [~,best]=max(score);
        centers(k)=candidates(best);
    end
    item=audit_design(cfg,"optimized_N"+nb,centers,edges);
    if isempty(items),items=item;else,items(end+1,1)=item;end %#ok<AGROW>
end

current_edges=linspace(cfg.fov_deg(1),cfg.fov_deg(2),numel(cfg.focus_angles_deg)+1);
items(end+1,1)=audit_design(cfg,"implemented_N6",cfg.focus_angles_deg,current_edges);
summary=struct2table(items);
writetable(summary,fullfile(out_dir,'A19_sector_count_feasibility.csv'));
parameter_snapshot(cfg,fullfile(out_dir,'A19_parameters.json'));
save(fullfile(out_dir,'A19_workspace.mat'),'cfg','summary','-v7.3');
disp(summary);
end

function worst_db=sector_worst_coverage(cfg,center_deg,edges_deg)
angles=cfg.response_angle_grid_deg;
inside=angles>=edges_deg(1)&angles<=edges_deg(2);
worst=inf;
local_cfg=cfg;
local_cfg.focus_angles_deg=center_deg;
for f=cfg.frequency_grid
    tx=tx_responses(local_cfg,f,angles);
    magnitude=abs(tx(1,:));
    magnitude=magnitude/(max(magnitude)+eps);
    worst=min(worst,min(magnitude(inside)));
end
worst_db=20*log10(worst+eps);
end

function item=audit_design(cfg,label,centers,edges)
% Final audit: include every prescribed sector boundary explicitly, then add
% a dense uniform FOV grid.  Receive aliases are evaluated at their exact
% continuous angles rather than rounded to the audit grid.
angles=unique([cfg.fov_deg(1):cfg.final_audit_angle_step_deg:cfg.fov_deg(2), ...
    edges(:).']);
local_cfg=cfg;
local_cfg.focus_angles_deg=centers;
worst_coverage=inf;
worst_grating=-inf;
worst_grating_frequency=NaN;
worst_grating_target=NaN;
worst_grating_angle=NaN;
worst_coverage_frequency=NaN;
worst_coverage_angle=NaN;
for f=cfg.final_audit_frequency_grid
    tx=tx_responses(local_cfg,f,angles);
    for k=1:numel(centers)
        assigned=angles>=edges(k)&angles<=edges(k+1);
        target_angles=angles(assigned);
        full_beam=abs(tx(k,:));
        full_beam=full_beam/(max(full_beam)+eps);
        [coverage_at_f,local_index]=min(full_beam(assigned));
        if coverage_at_f<worst_coverage
            worst_coverage=coverage_at_f;
            worst_coverage_frequency=f/1e3;
            worst_coverage_angle=target_angles(local_index);
        end
        desired=abs(single_tx_response(cfg,f,centers(k),target_angles));
        desired=max(desired,eps);
        for q=cfg.receive_alias_orders
            u_alias=sind(target_angles)+q*cfg.fc/(cfg.rx_spacing_lambda*f);
            visible=abs(u_alias)<=1;
            alias_angles=nan(size(target_angles));
            alias_angles(visible)=asind(u_alias(visible));
            visible=visible & alias_angles>=cfg.fov_deg(1) & ...
                alias_angles<=cfg.fov_deg(2);
            if ~any(visible),continue;end
            alias_response=abs(single_tx_response(cfg,f,centers(k), ...
                alias_angles(visible)));
            ratio_db=20*log10((alias_response+eps)./desired(visible));
            [candidate,index]=max(ratio_db);
            if candidate>worst_grating
                visible_targets=target_angles(visible);
                visible_aliases=alias_angles(visible);
                worst_grating=candidate;
                worst_grating_frequency=f/1e3;
                worst_grating_target=visible_targets(index);
                worst_grating_angle=visible_aliases(index);
            end
        end
    end
end
coverage_db=20*log10(worst_coverage+eps);
item.Design=label;
item.SectorCount=numel(centers);
item.SteeringAnglesDeg=strjoin(compose('%.2f',centers),';');
item.WorstAssignedCoverageDb=coverage_db;
item.WorstCoverageFrequencyKHz=worst_coverage_frequency;
item.WorstCoverageAngleDeg=worst_coverage_angle;
item.WorstVisibleGratingRatioDb=worst_grating;
item.WorstGratingFrequencyKHz=worst_grating_frequency;
item.WorstGratingTargetDeg=worst_grating_target;
item.WorstGratingAngleDeg=worst_grating_angle;
item.CoveragePass=coverage_db>=cfg.coverage_threshold_db;
item.GratingPass=worst_grating<=cfg.grating_threshold_db;
item.Feasible=item.CoveragePass&&item.GratingPass;
item.DesignCandidateSpacingDeg=cfg.design_angle_step_deg;
item.SearchAuditAngleSpacingDeg=cfg.search_audit_angle_step_deg;
item.SearchFrequencySpacingKHz=median(diff(cfg.search_frequency_grid))/1e3;
item.FinalAuditAngleSpacingDeg=cfg.final_audit_angle_step_deg;
item.FinalAuditFrequencySpacingKHz=median(diff(cfg.final_audit_frequency_grid))/1e3;
end

function response=single_tx_response(cfg,frequency_hz,focus_deg,angles_deg)
lambda_c=cfg.c/cfg.fc;
k=2*pi*frequency_hz/cfg.c;
if lower(string(cfg.steering_mode))=="ttd"
    k_steer=k;
else
    k_steer=2*pi*cfg.fc/cfg.c;
end
positions=((0:cfg.n_tx-1)-(cfg.n_tx-1)/2).' * ...
    cfg.tx_spacing_lambda*lambda_c;
weights=local_window("hamming",cfg.n_tx);
weights=weights/sum(weights);
steering=exp(-1j*k_steer*positions*sind(focus_deg));
response=(weights.*steering).'*exp(1j*k*positions*sind(angles_deg));
if isfield(cfg,'include_element_pattern')&&cfg.include_element_pattern
    x=(cfg.tx_element_width_lambda*lambda_c*frequency_hz/cfg.c).*sind(angles_deg);
    element=ones(size(x));
    nonzero=abs(x)>eps;
    element(nonzero)=sin(pi*x(nonzero))./(pi*x(nonzero));
    response=response.*element;
end
end

function tx=tx_responses(cfg,frequency_hz,angles_deg)
% Compute only the transmit response; the common helper also forms a large
% receive scan matrix that is unnecessary in this steering search.
lambda_c=cfg.c/cfg.fc;
k=2*pi*frequency_hz/cfg.c;
if lower(string(cfg.steering_mode))=="ttd"
    k_steer=k;
else
    k_steer=2*pi*cfg.fc/cfg.c;
end
positions=((0:cfg.n_tx-1)-(cfg.n_tx-1)/2).' * ...
    cfg.tx_spacing_lambda*lambda_c;
weights=local_window("hamming",cfg.n_tx);
weights=weights/sum(weights);
u=sind(angles_deg);
manifold=exp(1j*k*positions*u);
tx=complex(zeros(numel(cfg.focus_angles_deg),numel(angles_deg)));
for kbeam=1:numel(cfg.focus_angles_deg)
    steering=exp(-1j*k_steer*positions*sind(cfg.focus_angles_deg(kbeam)));
    tx(kbeam,:)=(weights.*steering).'*manifold;
end
if isfield(cfg,'include_element_pattern')&&cfg.include_element_pattern
    x=(cfg.tx_element_width_lambda*lambda_c*frequency_hz/cfg.c).*u;
    element=ones(size(x));
    nonzero=abs(x)>eps;
    element(nonzero)=sin(pi*x(nonzero))./(pi*x(nonzero));
    tx=tx.*element;
end
end
