"""Author-approved detector post-processing, without changing primary results.

All nine physical NMS combinations are frozen before this script computes its
outcomes. Original null thresholds, source images, seeds and association gate
are retained. No new physical echo or noise is generated.
"""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import csv,json,hashlib,sys,time
import numpy as np
from scipy import ndimage,optimize,stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import argparse
ROOT=Path(__file__).resolve().parents[2]
HERE=ROOT/'data'
OUT=ROOT/'outputs/nms_sensitivity'
METHODS=('dense','cdm','fdm')
WIDTHS_R=(.005,.01,.03)
WIDTHS_A=(.5,1.,2.)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,o):Path(p).write_text(json.dumps(o,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
def table(p,rows):
    with Path(p).open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def read(p):
    with Path(p).open(encoding='utf-8-sig') as f:return list(csv.DictReader(f))

def detect(image,r,a,threshold,tr,ta,half_r,half_a):
    p=abs(image);dr=float(np.median(np.diff(r)));da=float(np.median(np.diff(a)))
    sr=int(np.ceil(half_r/dr-1e-9));sa=int(np.ceil(half_a/da-1e-9))
    local=ndimage.maximum_filter(p,size=(2*sr+1,2*sa+1),mode='constant',cval=-np.inf)
    candidates=np.argwhere((p>=threshold)&(p==local));selected=[]
    if not len(candidates):return 0.,0,0,None,None
    order=np.argsort(p[candidates[:,0],candidates[:,1]])[::-1]
    for index in order:
        ir,ia=candidates[index]
        if any(abs(r[ir]-r[jr])<=half_r+1e-12 and abs(a[ia]-a[ja])<=half_a+1e-12 for jr,ja in selected):continue
        selected.append((ir,ia))
    found=np.asarray(selected);rr=r[found[:,0]];aa=a[found[:,1]]
    delta_r=abs(tr[:,None]-rr);delta_a=abs(ta[:,None]-aa)
    admissible=(delta_r<=.02+1e-12)&(delta_a<=1.+1e-12)
    cost=np.sqrt((delta_r/.02)**2+(delta_a/1.)**2);cost[~admissible]=1e9
    it,ic=optimize.linear_sum_assignment(cost);matched=int(np.sum(admissible[it,ic]))
    accepted=admissible[it,ic]
    mean_angle=float(np.mean(delta_a[it[accepted],ic[accepted]])) if matched else None
    mean_range_mm=float(np.mean(delta_r[it[accepted],ic[accepted]])*1000) if matched else None
    return matched/len(tr),len(selected)-matched,len(selected),mean_angle,mean_range_mm

def run(noise_cache):
    OUT.mkdir(parents=True,exist_ok=True)
    inputs=[HERE/'point30m/clean_images.npz',noise_cache,
        HERE/'point30m/all_trials.csv',HERE/'point30m/false_alarm_summary.csv',Path(__file__)]
    protocol={'purpose':'postprocess all nine NMS settings; primary frozen detector unchanged',
        'range_halfwidths_m':list(WIDTHS_R),'angle_halfwidths_deg':list(WIDTHS_A),
        'primary_halfwidths':{'range_m':.03,'angle_deg':2},
        'association_admissibility':'rectangle abs(range-error)<=.02m AND abs(angle-error)<=1deg',
        'association_cost':'sqrt((range-error/.02)^2+(angle-error/1)^2); one-to-one Hungarian',
        'localization':'per-frame mean absolute angle/range errors over accepted one-to-one matched targets only; not error of missed targets',
        'threshold':'same independent60train method-specific threshold, multiplied by original sigma; no retuning',
        'physical_echo':'same cached complex clean and30 unit-noise images, no new generation',
        'cached_precision':'noise images complex64, clean complex128; primary Pd/FP exact-match audit against original complex128 outcomes',
        'source_sha256':{str(p):sha(p) for p in inputs},'created_before_outcomes':True,
        'no_best_setting_selection':True,'author_approved_additional_postprocess':True}
    locked=OUT/'protocol_before_postprocess.json'
    if locked.exists():
        if json.loads(locked.read_text(encoding='utf-8'))!=protocol:raise RuntimeError('Postprocess source changed after freeze')
    else:dump(locked,protocol)
    source=np.load(inputs[0]);bank=np.load(inputs[1]);original=read(inputs[2]);r=source['ranges_m'];a=source['angles_deg'];tr=source['true_ranges_m'];ta=source['true_angles_deg']
    clean_cache={m:source[m] for m in METHODS}
    noise_cache={m:bank[m] for m in METHODS}
    unit_threshold={q['method']:float(q['unit_noise_threshold']) for q in read(inputs[3])}
    rows=[];started=time.time();same=[]
    for index,row in enumerate(original):
        m=row['method'];rep=int(row['noise_repetition']);sigma=float(row['pre_frontend_noise_sigma'])
        clean=clean_cache[m];noise=noise_cache[m][rep].astype(complex);image=clean+sigma*noise
        threshold=unit_threshold[m]*sigma
        for half_r in WIDTHS_R:
            for half_a in WIDTHS_A:
                pd,fp,total,angle_error,range_error=detect(image,r,a,threshold,tr,ta,half_r,half_a)
                rows.append({'method':m,'noise_normalization':row['noise_normalization'],'nominal_snr_db':float(row['nominal_snr_db']),
                    'noise_repetition':rep,'noise_seed':row['noise_seed'],'NMS_range_halfwidth_m':half_r,
                    'NMS_angle_halfwidth_deg':half_a,'detection_probability':pd,'false_detections':fp,'total_detections':total,
                    'matched_mean_absolute_angle_error_deg':angle_error,'matched_mean_absolute_range_error_mm':range_error})
                if half_r==.03 and half_a==2.:
                    same.append({'method':m,'noise_normalization':row['noise_normalization'],'nominal_snr_db':float(row['nominal_snr_db']),
                        'noise_repetition':rep,'Pd_equal':abs(pd-float(row['detection_probability']))<=1e-12,
                        'FP_equal':fp==int(row['false_detections']),'total_equal':total==int(row['total_detections'])})
        if index%90==0:print('NMS samples',index+1,'/',len(original),flush=True)
    table(OUT/'all_combinations_all_trials.csv',rows);table(OUT/'primary_cached_precision_audit.csv',same)
    group={}
    for row in rows:
        key=tuple(row[k] for k in ['method','noise_normalization','nominal_snr_db','NMS_range_halfwidth_m','NMS_angle_halfwidth_deg'])
        group.setdefault(key,[]).append(row)
    agg=[]
    for key,g in group.items():
        q=dict(zip(['method','noise_normalization','nominal_snr_db','NMS_range_halfwidth_m','NMS_angle_halfwidth_deg'],key));q['replications']=len(g)
        for metric in ['detection_probability','false_detections','total_detections','matched_mean_absolute_angle_error_deg','matched_mean_absolute_range_error_mm']:
            x=np.array([row[metric] for row in g if row[metric] is not None]);q[metric]=float(x.mean()) if len(x) else None
            q[metric+'_ci95_halfwidth']=float(stats.t.ppf(.975,len(x)-1)*x.std(ddof=1)/np.sqrt(len(x))) if len(x)>1 else None
        agg.append(q)
    table(OUT/'summary_mean_ci.csv',agg)
    plt.rcParams.update({'font.family':'serif','font.size':10,'pdf.fonttype':42})
    fig,axes=plt.subplots(2,3,figsize=(7.15,4.35))
    for j,m in enumerate(METHODS):
        q=[row for row in agg if row['method']==m and row['noise_normalization']=='common_per_element_noise' and row['nominal_snr_db']==0]
        for i,metric in enumerate(('detection_probability','false_detections')):
            z=np.asarray([[next(row[metric] for row in q if row['NMS_range_halfwidth_m']==rr and row['NMS_angle_halfwidth_deg']==aa) for aa in WIDTHS_A] for rr in WIDTHS_R])
            im=axes[i,j].imshow(z,origin='lower',aspect='auto',cmap='cividis',vmin=0,vmax=1 if i==0 else None)
            for ir in range(3):
                for ia in range(3):axes[i,j].text(ia,ir,f'{z[ir,ia]:.3f}' if i==0 else f'{z[ir,ia]:.1f}',ha='center',va='center',fontsize=9.6,color='white' if z[ir,ia]<np.mean(z) else 'black')
            axes[i,j].set_xticks(range(3),['0.5','1','2']);axes[i,j].set_yticks(range(3),['5','10','30'])
            axes[i,j].set_xlabel('NMS angle halfwidth (deg)')
            if j==0:axes[i,j].set_ylabel('NMS range halfwidth (mm)')
            if i==0:axes[i,j].set_title(m.upper())
            fig.colorbar(im,ax=axes[i,j],label='Pd' if i==0 else 'Unassigned peaks')
    fig.tight_layout()
    for ext in ('pdf','png','svg'):fig.savefig(OUT/('all_nms_0dB.'+ext),dpi=600,bbox_inches='tight')
    plt.close(fig)
    audit={'status':'computed_pending_author_validation','rows':len(rows),'settings':9,
        'primary_precision_audit_samples':len(same),'primary_precision_audit_failures':sum(not(q['Pd_equal'] and q['FP_equal'] and q['total_equal']) for q in same),
        'elapsed_seconds':time.time()-started,'unchanged_primary_results':True,'no_best_setting_selected':True,
        'interpretation':'NMS-induced suppression is distinct from waveform/beam resolution; unassigned signal-scene peaks are distinct from noise-only Pfa'}
    dump(OUT/'receipt.json',audit);print(audit,flush=True)
if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--noise-cache',type=Path,required=True,
        help='Full independent_test_noise_images.npz, obtained from a point rerun or the author.')
    options=parser.parse_args()
    if not options.noise_cache.is_file(): parser.error('The full noise image cache does not exist.')
    run(options.noise_cache.resolve())
