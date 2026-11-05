function rows = compute_detection_assignment_metrics(image_data, range_axis_m, angle_axis_deg, truth_catalog, thresholds_db, varargin)
%COMPUTE_DETECTION_ASSIGNMENT_METRICS One-to-one target/peak evaluation.
% Peaks are extracted by regional-maximum detection followed by physical
% range--angle non-maximum suppression. A global linear assignment is then
% solved, and only pairs inside the declared elliptical gate are accepted.

p = inputParser;
addParameter(p,'NMSRangeM',0.015);
addParameter(p,'NMSAngleDeg',1.0);
addParameter(p,'GateRangeM',0.04);
addParameter(p,'GateAngleDeg',2.0);
parse(p,varargin{:});
o=p.Results;

mag=abs(image_data);
mag=mag/(max(mag(:))+eps);
db=20*log10(mag+eps);
truth_range=truth_catalog.RangeM(:);
truth_angle=truth_catalog.AngleDeg(:);
rows=repmat(struct(),numel(thresholds_db),1);

for it=1:numel(thresholds_db)
    threshold_db=thresholds_db(it);
    regional=imregionalmax(mag,8) & db>=threshold_db;
    [rr,aa]=find(regional);
    peak_amp=mag(regional);
    [peak_amp,ord]=sort(peak_amp,'descend');
    rr=rr(ord); aa=aa(ord);

    keep=false(size(rr)); kept=[];
    for k=1:numel(rr)
        if isempty(kept)
            accept=true;
        else
            dr=abs(range_axis_m(rr(k))-range_axis_m(rr(kept)));
            da=abs(angle_axis_deg(aa(k))-angle_axis_deg(aa(kept)));
            accept=all((dr/o.NMSRangeM).^2+(da/o.NMSAngleDeg).^2>1);
        end
        if accept
            keep(k)=true;
            kept(end+1)=k; %#ok<AGROW>
        end
    end
    rr=rr(keep); aa=aa(keep); peak_amp=peak_amp(keep);
    cand_range=range_axis_m(rr(:)); cand_range=cand_range(:);
    cand_angle=angle_axis_deg(aa(:)); cand_angle=cand_angle(:);

    if isempty(cand_range)
        pairs=zeros(0,2); pair_cost=[];
    else
        cost=sqrt(((truth_range-cand_range.')./o.GateRangeM).^2 + ...
                  ((truth_angle-cand_angle.')./o.GateAngleDeg).^2);
        [pairs,~]=matchpairs(cost,2,'min');
        pair_cost=cost(sub2ind(size(cost),pairs(:,1),pairs(:,2)));
        valid=pair_cost<=1;
        pairs=pairs(valid,:); pair_cost=pair_cost(valid);
    end

    if isempty(pairs)
        angle_abs=[]; range_abs=[];
    else
        angle_abs=abs(truth_angle(pairs(:,1))-cand_angle(pairs(:,2)));
        range_abs=abs(truth_range(pairs(:,1))-cand_range(pairs(:,2)));
    end

    rows(it).ThresholdDb=threshold_db;
    rows(it).TargetCount=numel(truth_range);
    rows(it).CandidateCount=numel(cand_range);
    rows(it).MatchedCount=size(pairs,1);
    rows(it).MissedCount=numel(truth_range)-size(pairs,1);
    rows(it).FalsePeakCount=numel(cand_range)-size(pairs,1);
    rows(it).DetectionFraction=size(pairs,1)/numel(truth_range);
    rows(it).MeanAngleAbsErrorDeg=mean(angle_abs,'omitnan');
    rows(it).RMSEAngleDeg=sqrt(mean(angle_abs.^2,'omitnan'));
    rows(it).MeanRangeAbsErrorMm=1000*mean(range_abs,'omitnan');
    rows(it).RMSERangeMm=1000*sqrt(mean(range_abs.^2,'omitnan'));
    rows(it).MeanAcceptedNormalizedCost=mean(pair_cost,'omitnan');
    rows(it).StrongestUnassignedPeakDb=strongest_unassigned(peak_amp,pairs);
end
end

function value=strongest_unassigned(peak_amp,pairs)
if isempty(peak_amp), value=-Inf; return; end
assigned=false(numel(peak_amp),1);
if ~isempty(pairs),assigned(pairs(:,2))=true;end
if all(assigned),value=-Inf;else,value=20*log10(max(peak_amp(~assigned))+eps);end
end
