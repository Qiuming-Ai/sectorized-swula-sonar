"""Frozen final main-scene checks. Originals are read-only.

All array phases use physical frequency exactly once. Batched Bluestein/CZT
evaluation is an algebraic implementation of the explicit uniform-array sum,
not an approximate beam mask. The consistency gates precede formal outcomes.
"""
from __future__ import annotations
import os
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[name]='1'
import argparse, csv, hashlib, json, shutil, sys, time, traceback
from pathlib import Path
from dataclasses import asdict
from types import SimpleNamespace
from collections import defaultdict
import numpy as np
from scipy import fft, signal, stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

PROJECT=Path(__file__).resolve().parents[2]
CORE=PROJECT/'code'
SNAP=CORE/'common'
HERE=PROJECT/'outputs/scenes'
INPUTS=PROJECT/'data/inputs'
SHIP_ARRIVALS=PROJECT/'data/shipwreck/geometry/A11_raytrace_arrivals.csv'
METHODS=('dense','cdm','fdm')
LABELS={'dense':'Dense ULA','cdm':'CDM SW-ULA','fdm':'FDM SW-ULA'}
COLORS={'dense':'#0072B2','cdm':'#D55E00','fdm':'#009E73'}
STYLES={'dense':('-', 'o'),'cdm':('--','s'),'fdm':(':','^')}
plt.rcParams.update({'font.family':'serif','font.size':10,'axes.labelsize':10,
    'legend.fontsize':8,'axes.spines.top':False,'axes.spines.right':False,
    'pdf.fonttype':42,'ps.fonttype':42})

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,o):
    Path(p).parent.mkdir(parents=True,exist_ok=True)
    Path(p).write_text(json.dumps(o,ensure_ascii=False,indent=2,allow_nan=False,
        default=lambda x:x.item() if isinstance(x,np.generic) else str(x))+'\n',encoding='utf-8')
def table(p,rows):
    if not rows:return
    keys=list(dict.fromkeys(k for row in rows for k in row))
    with Path(p).open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
def savefig(fig,p):
    fig.tight_layout()
    for ext in ('pdf','png','svg'):fig.savefig(str(p)+'.'+ext,dpi=600,bbox_inches='tight')
    plt.close(fig)

PROTOCOL={
 'id':'RINENG_FINAL_MAIN_SCENES_v01_20261010',
 'created_before_outcomes':True,'no_favorable_result_gate':True,
 'originals_read_only':True,'primary_energy':'actual simultaneous Tx-element sum energy=1',
 'frontend':'common final 129-tap Kaiser8.6 Tx centre crop and common guarded Rx FIR',
 'cdm':'frequency ZC -> IFFT -> Kaiser5 polyphase; no temporal taper',
 'point':{'source':'08_A8_review_imaging_figures/data/target_coordinates.csv',
   'coordinate_shift_m':25.1,'target_amplitude':'all original amplitudes=1, phase=0; no reflectivity-energy rescale',
   'physical_geometry':'exact spherical Tx-target-Rx generation, plane-wave Rx focus',
   'range_ROI_m':[30,30.8],'range_step_m':.0025,'angle_step_deg':.25,
   'noise_SNR_db':[-15,-10,-5,0],'noise_repetitions':30,
   'pure_noise_train_blocks':60,'pure_noise_independent_test_blocks':30,
   'threshold':'method-specific95th percentile higher of60 pure-noise image maxima',
   'threshold_train_seed_base':20261015100,'noise_test_seed_base':20261015200,
   'full_ROI_detector':'same frozen physical NMS +/-0.03m,+/-2deg, one-to-one matching +/-0.02m,+/-1deg',
   'close_pair_warning':'20mm axial pair lies inside the frozen30mm NMS halfwidth; Pd is not an isolated resolution limit'},
 'dense':{'source':'A10_workspace.mat scenes{2}; exact saved complex coherent cells',
   'generative_recipe':'1000 reflectors,50 requested distance layers,40dB amplitude dynamic range, phase U(0,2pi)',
   'coordinate_shift_m':25.1,'original_range_ROI_m':[4.92,5.52],
   'range_ROI_m':[30.02,30.62],'source_angle_grid_deg':1,
   'model_grid':'201 equal sin(theta) nodes; source-angle quantization recorded',
   'noise_SNR_db':[-15,-10,-5,0],'noise_repetitions':30,'noise_seed_base':20261016100,
   'repeat_identity':'one original fixed geometry and phase realization,30 independent noise blocks'},
 'ship':{'source':'A11_raytrace_arrivals.csv RangeM,AngleDeg,MetricAmplitude',
   'range_ROI_m':[42,122],'no_geometry_or_view_change':True,
   'model_grid':'201 equal sin(theta) nodes; exact physical polar view (x=r sin theta,z=r cos theta)',
   'phase_repetitions':30,'phase_seed_base':20250821,'noise_seed_base':20261017100,
   'SNR_db':10,'amplitude':'unaltered MetricAmplitude; no metric TVG/gamma/smoothing',
   'repeat_identity':'one original ship/seabed geometry,30 independently generated arrival phases and paired receiver-noise blocks'},
 'noise_rules':['common_per_element_noise','equal_raw_channel_snr'],
 'noise_primary':'same per-channel CN(0,sigma^2) before common Rx frontend, same standardized128channel draw; sparse first32',
 'noise_SNR':'common-noise nominal SNR anchored to noiseless Dense mean raw-channel power; equalrawSNR supplemental',
 'metrics':'linear magnitude/power on fixed physical ROI; noiseless Dense reconstruction primary reference',
 'guard':{'range_halfwidth_m':.06,'angle_halfwidth_deg':1,'minimum_background_fraction':.05},
 'consistency_gates':{'CZT_vs_explicit_relative_L2':1e-10,'CZT_grid_vs_core_image_relative_L2':1e-10},
 'scope':'numerical waveform/array simulation; geometry-derived ship stress test, not acoustic acquisition',
}

def freeze(suite):
    sources=[CORE/'common'/x for x in ('numerical_model.py','scene_noise.py','wave_point.py')]
    sources += [CORE/'core_ready.json',Path(__file__)]
    if suite in ('point','all'): sources.append(INPUTS/'target_coordinates.csv')
    if suite in ('dense','all'): sources.append(INPUTS/'dense1000x50_source.npz')
    if suite == 'ship': sources.append(SHIP_ARRIVALS)
    gate=json.loads((CORE/'core_ready.json').read_text(encoding='utf-8'))
    if gate['model_sha256']!=sha(CORE/'common/numerical_model.py'):raise RuntimeError('Core gate mismatch')
    current={str(p):sha(p) for p in sources}
    existing=HERE/f'protocol_before_{suite}.json'
    if existing.exists():
        locked=json.loads(existing.read_text(encoding='utf-8'))
        if locked['source_sha256']!=current:raise RuntimeError('Source changed after freeze: preserve and create a new version')
    else:
        dump(existing,{'protocol':PROTOCOL,'source_sha256':current,'frozen_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())})
    sys.path.insert(0,str(SNAP))
    import numerical_model as nm
    import scene_noise as sn
    sn.ROOT=CORE;sn._core.cache_clear()
    return nm,sn,current

def czt_rows(x,q,a,m):
    """Each row sum x[n] exp(i*a*n) exp(i*q*n*j), j=0..m-1."""
    n=x.shape[1];pad=fft.next_fast_len(n+m-1);out=np.empty((len(x),m),complex)
    ni=np.arange(n);ji=np.arange(m);lag=np.arange(-n+1,m)
    for lo in range(0,len(x),128):
        sl=slice(lo,min(lo+128,len(x)));qq=q[sl,None];aa=a[sl,None]
        xx=x[sl]*np.exp(1j*(aa*ni+qq*ni**2/2))
        kernel=np.exp(-1j*qq*lag**2/2)
        conv=fft.ifft(fft.fft(xx,pad,axis=1)*fft.fft(kernel,pad,axis=1),axis=1)
        out[sl]=conv[:,n-1:n-1+m]*np.exp(1j*qq*ji**2/2)
    return out

class GridChain:
    def __init__(self,nm,cfg,method):
        self.nm,self.cfg,self.method=nm,cfg,method
        self.nr,self.na=len(cfg.ranges),len(cfg.u)
        self.nf=fft.next_fast_len(self.nr+2*len(nm.waveforms(cfg,method)[0])-2)
        self.send,self.matched,self.nu,self.meta=nm.spectra(cfg,method,self.nf)
        self.f=cfg.fc+self.nu;self.n=cfg.n_dense if method=='dense' else cfg.n_sparse
        self.d=(.5 if method=='dense' else 2)*cfg.c/cfg.fc
        self.scale=2*np.pi*self.f*self.d/cfg.c
        self.du=cfg.u[1]-cfg.u[0];self.u0=cfg.u[0]
        self.sec=np.clip(np.digitize(cfg.angles,cfg.edges[1:-1]),0,5)
        self.front=nm.frontend_response(cfg,self.nf)
        self.carrier=np.exp(-4j*np.pi*cfg.fc*(cfg.ranges-cfg.gate_range)/cfg.c)
        self.reverse=np.conj(self.carrier)
    def raw(self,scene):
        sc=fft.fft(scene*self.carrier[:,None],self.nf,axis=0)
        if self.method=='dense':illum=self.send[:,0,None]
        else:illum=np.sum(self.send[:,:,None]*self.nm.tx_pattern(self.cfg,self.f,self.cfg.angles),axis=1)
        source=sc*illum*np.sinc(self.cfg.element_width_lambda/self.cfg.fc*self.f[:,None]*self.cfg.u[None,:])
        # Sum source nodes for centred physical receiver positions.
        q=self.scale*self.du;a=-q*(self.n-1)/2
        raw=czt_rows(source,q,a,self.n)
        raw*=np.exp(1j*self.scale[:,None]*self.u0*(np.arange(self.n)[None,:]-(self.n-1)/2))
        return raw
    def reconstruct(self,raw):
        # Exact wideband uniform Rx sum; envelope delay and carrier phase included once.
        q=-self.scale*self.du;a=-self.scale*self.u0
        scan=czt_rows(raw*self.nm.weights(self.n)[None,:],q,a,self.na)
        scan*=np.exp(1j*self.scale[:,None]*(self.n-1)/2*self.cfg.u[None,:])
        matched=self.matched[:,0,None] if self.method=='dense' else self.matched[:,self.sec]
        image=fft.ifft(scan*self.front[:,None]*np.conj(matched),axis=0)[:self.nr]
        return image*self.reverse[:,None]

def gates(nm):
    out=HERE/'implementation_gates';out.mkdir(exist_ok=True);rows=[]
    gen=np.random.default_rng(20261018001)
    x=gen.normal(size=(7,13))+1j*gen.normal(size=(7,13));q=gen.uniform(-2,2,7);a=gen.uniform(-2,2,7)
    fast=czt_rows(x,q,a,19)
    direct=np.asarray([[np.sum(xx*np.exp(1j*(aa+qq*j)*np.arange(13))) for j in range(19)] for xx,aa,qq in zip(x,a,q)])
    error=float(np.linalg.norm(fast-direct)/np.linalg.norm(direct))
    rows.append({'gate':'CZT_vs_explicit','relative_L2':error,'tolerance':1e-10,'passed':error<=1e-10})
    cfg=nm.Config(range_span=.15,angular_nodes=41)
    scene,_=nm.scene_grid(cfg,20,5,.1,20,20261018002)
    for method in METHODS:
        chain=GridChain(nm,cfg,method);raw=chain.raw(scene);fast=chain.reconstruct(raw);direct=nm.image_grid(cfg,scene,method)[0]
        error=float(np.linalg.norm(fast-direct)/np.linalg.norm(direct))
        rows.append({'gate':'CZT_grid_vs_core_'+method,'relative_L2':error,'tolerance':1e-10,'passed':error<=1e-10})
    table(out/'consistency.csv',rows);dump(out/'receipt.json',{'all_passed':all(r['passed'] for r in rows),'rows':rows})
    if not all(r['passed'] for r in rows):raise RuntimeError('Frozen implementation consistency gate failed; all diagnostics retained')
    print('GATES PASS',rows,flush=True)

KEYS=('truth_magnitude_nmse_db','truth_magnitude_correlation','dense_magnitude_correlation','support_leakage',
      'cnr_linear_magnitude','noise_induced_magnitude_NMSE_db','range_profile_correlation_to_clean_Dense',
      'angle_profile_correlation_to_clean_Dense','detection_probability','false_detections')
def axis_config(r,a):return SimpleNamespace(ranges=r,angles=a,dr=float(np.median(np.diff(r))),range_origin=float(r[0]),range_span=float(r[-1]-r[0]))
def evaluate(nm,truth,image,cfg,dense,clean=None):
    q=nm.metrics(truth,image,cfg,dense=dense)
    valid=np.broadcast_to(((cfg.ranges>=cfg.range_origin)&(cfg.ranges<=cfg.range_origin+cfg.range_span))[:,None],image.shape)
    def rho(x,y):
        x=x.ravel()-x.mean();y=y.ravel()-y.mean();den=np.linalg.norm(x)*np.linalg.norm(y)
        return float(x@y/den) if den>0 else None
    mag=abs(image);ref=abs(dense)
    q['range_profile_correlation_to_clean_Dense']=rho(np.sum(mag*valid,axis=1),np.sum(ref*valid,axis=1))
    q['angle_profile_correlation_to_clean_Dense']=rho(np.sum(mag*valid,axis=0),np.sum(ref*valid,axis=0))
    if clean is not None:
        x=mag[valid];y=abs(clean)[valid];fit=float(x@y/max(y@y,1e-30));err=np.sum((x-fit*y)**2)/max(np.sum((fit*y)**2),1e-30)
        q['noise_induced_magnitude_NMSE_db']=float(10*np.log10(max(err,1e-30)));q['fitted_clean_scale']=fit
    return q
def aggregate(rows,by):
    grouped=defaultdict(list)
    for row in rows:grouped[tuple(row[k] for k in by)].append(row)
    result=[]
    for key,g in grouped.items():
        q=dict(zip(by,key));q['replications']=len(g)
        for metric in KEYS:
            v=np.array([row[metric] for row in g if row.get(metric) is not None],float)
            if not len(v):continue
            q[metric]=float(v.mean());q[metric+'_sd']=float(v.std(ddof=1)) if len(v)>1 else None
            q[metric+'_ci95_halfwidth']=float(stats.t.ppf(.975,len(v)-1)*v.std(ddof=1)/np.sqrt(len(v))) if len(v)>1 else None
        result.append(q)
    return result

def panel(out,name,r,a,images,common=False,sector=False,truth_r=None,truth_a=None):
    fig,axes=plt.subplots(1,3,figsize=(10.5,3.7))
    ref=max(abs(x).max() for x in images.values()) if common else None
    xx=r[:,None]*np.sin(np.deg2rad(a))[None,:];zz=r[:,None]*np.cos(np.deg2rad(a))[None,:]
    for ax,method in zip(axes,METHODS):
        z=20*np.log10(np.maximum(abs(images[method]),1e-30)/(float(abs(images[method]).max()) if ref is None else ref))
        if sector:
            mesh=ax.pcolormesh(xx,zz,np.maximum(z,-50),cmap='viridis',vmin=-50,vmax=0,shading='auto',rasterized=True)
            ax.set(xlabel='Cross-range x (m)',ylabel='Forward range z (m)',aspect='equal')
        else:
            mesh=ax.pcolormesh(a,r,np.maximum(z,-50),cmap='viridis',vmin=-50,vmax=0,shading='auto',rasterized=True)
            ax.set(xlabel='Angle (deg)',ylabel='Range (m)',xlim=(-50,50))
        if truth_r is not None:ax.scatter(truth_a,truth_r,facecolors='none',edgecolors='white',s=15,linewidth=.6)
        ax.set_title(LABELS[method],fontsize=10)
        fig.colorbar(mesh,ax=ax,label='Amplitude (dB)',shrink=.85)
    savefig(fig,out/name)

def curves(out,name,rows,rule):
    selected=[r for r in rows if r['noise_normalization']==rule]
    keys=[('noise_induced_magnitude_NMSE_db','Magnitude perturbation NMSE (dB)'),
          ('dense_magnitude_correlation','Noise-free Dense correlation'),('support_leakage','Support leakage')]
    fig,axes=plt.subplots(1,3,figsize=(10.5,3.1))
    for ax,(metric,label) in zip(axes,keys):
        for method in METHODS:
            q=sorted([r for r in selected if r['method']==method],key=lambda r:r['nominal_snr_db'])
            if not q or metric not in q[0]:continue
            ls,mark=STYLES[method]
            ax.errorbar([r['nominal_snr_db'] for r in q],[r[metric] for r in q],
                yerr=[r[metric+'_ci95_halfwidth'] or 0 for r in q],color=COLORS[method],ls=ls,marker=mark,capsize=2,label=LABELS[method])
        ax.set(xlabel='Nominal input SNR (dB)',ylabel=label);ax.grid(alpha=.15);ax.legend()
    savefig(fig,out/name)

def point(nm,sn):
    out=HERE/'point30m';out.mkdir(exist_ok=True);cfg=nm.Config();nf=8192
    source=np.loadtxt(INPUTS/'target_coordinates.csv',delimiter=',',skiprows=1)
    tr=source[:,1]+25.1;ta=source[:,0];coef=source[:,2].astype(complex)
    r=np.arange(30,30.800001,.0025);a=np.arange(-50,50.000001,.25);acfg=axis_config(r,a)
    truth=np.zeros((len(r),len(a)),complex)
    for rr,aa,c in zip(tr,ta,coef):truth[np.argmin(abs(r-rr)),np.argmin(abs(a-aa))]+=c
    clean={};raws={};matched={};powers={};meta={}
    for m in METHODS:
        raws[m],matched[m],meta[m]=nm.point_raw(cfg,m,tr,ta,coef,geometry='spherical',nfft=nf)
        clean[m]=nm.reconstruct_points(cfg,m,raws[m],matched[m],r,a)
        powers[m]=float(np.mean(abs(fft.ifft(raws[m],axis=0))**2))
    table(out/'target_coordinates_30m.csv',[{'target':i+1,'old_range_m':float(source[i,1]),'range_m':float(rr),'angle_deg':float(aa),'coefficient_real':float(coef[i].real),'coefficient_imag':0.} for i,(rr,aa) in enumerate(zip(tr,ta))])
    cleanrows=[{'method':m,**evaluate(nm,truth,clean[m],acfg,clean['dense'])} for m in METHODS]
    table(out/'clean_metrics.csv',cleanrows)
    np.savez_compressed(out/'clean_images.npz',ranges_m=r,angles_deg=a,truth_complex=truth,true_ranges_m=tr,true_angles_deg=ta,**{m:clean[m] for m in METHODS})
    panel(out,'point30m_per_image',r,a,clean,truth_r=tr,truth_a=ta);panel(out,'point30m_common_scale',r,a,clean,True,truth_r=tr,truth_a=ta)
    cal=[];null=[];bank={m:[] for m in METHODS};beg=time.time()
    for rep in range(90):
        seed=(20261015100+rep) if rep<60 else (20261015200+rep-60)
        gen=np.random.default_rng(seed);noise=(gen.normal(size=(nf,128))+1j*gen.normal(size=(nf,128)))/np.sqrt(2)
        nfnoise=fft.fft(noise,axis=0)
        for m in METHODS:
            nr=128 if m=='dense' else 32;im=nm.reconstruct_points(cfg,m,nfnoise[:,:nr],matched[m],r,a)
            maximum=float(abs(im).max())
            if rep<60:cal.append({'method':m,'training_block':rep,'noise_seed':seed,'max_amplitude':maximum})
            else:bank[m].append(im);null.append({'method':m,'test_block':rep-60,'noise_seed':seed,'max_amplitude':maximum})
        if rep%5==0 or rep==89:print(f'POINT pure noise {rep+1}/90 {time.time()-beg:.1f}s',flush=True)
    threshold={m:float(np.quantile([q['max_amplitude'] for q in cal if q['method']==m],.95,method='higher')) for m in METHODS}
    for q in null:q['threshold']=threshold[q['method']];q['exceeded']=q['max_amplitude']>q['threshold']
    table(out/'pure_noise_training.csv',cal);table(out/'independent_false_alarm_blocks.csv',null)
    falses=[]
    for m in METHODS:
        count=sum(q['exceeded'] for q in null if q['method']==m)
        falses.append({'method':m,'test_blocks':30,'exceeded_blocks':count,'empirical_Pfa':count/30,'CI95_lower':0. if count==0 else float(stats.beta.ppf(.025,count,31-count)),
            'CI95_upper':1. if count==30 else float(stats.beta.ppf(.975,count+1,30-count)),'unit_noise_threshold':threshold[m]})
    table(out/'false_alarm_summary.csv',falses)
    rows=[]
    for rule in PROTOCOL['noise_rules']:
        for snr in PROTOCOL['point']['noise_SNR_db']:
            for m in METHODS:
                sigma=np.sqrt((powers['dense'] if rule=='common_per_element_noise' else powers[m])/10**(snr/10))
                for rep,n in enumerate(bank[m]):
                    image=clean[m]+sigma*n;pd,fp,peaks=sn._blind_detection(image,r,a,threshold[m]*sigma,tr,ta)
                    rows.append({'method':m,'noise_normalization':rule,'nominal_snr_db':snr,'noise_repetition':rep,'noise_seed':20261015200+rep,
                        'pre_frontend_noise_sigma':float(sigma),'actual_raw_channel_snr_db':float(10*np.log10(powers[m]/sigma**2)),
                        'detection_probability':pd,'false_detections':fp,'total_detections':peaks,**evaluate(nm,truth,image,acfg,clean['dense'],clean[m])})
            print('POINT metrics',rule,snr,flush=True)
    agg=aggregate(rows,['noise_normalization','nominal_snr_db','method']);table(out/'all_trials.csv',rows);table(out/'summary_mean_ci.csv',agg)
    for rule in PROTOCOL['noise_rules']:curves(out,'noise_'+rule,agg,rule)
    np.savez_compressed(out/'independent_test_noise_images.npz',ranges_m=r,angles_deg=a,**{m:np.asarray(bank[m],np.complex64) for m in METHODS})
    dump(out/'metadata.json',{'raw_signal_power_by_method':powers,'energy_metadata':meta,'noise_free_Dense_reference':True,'noise_repetitions_are_not_geometry_replications':True,
        'mask':'truth-nearest pixel for coefficient envelope; physical guard .06m/1deg','precision':'complex128 computations; complex64 stored noise bank only'})
    return {'rows':len(rows),'noise_blocks':90,'noise_computations':270,'false_alarm_summary':falses,'shape':[len(r),len(a)]}

def grid_suite(nm,kind):
    out=HERE/kind;out.mkdir(exist_ok=True);rows=[];cleanrows=[];scene_meta=[]
    if kind=='dense1000x50':
        cfg=nm.Config(range_origin=30.02,range_span=.6,angular_nodes=201)
        with np.load(INPUTS/'dense1000x50_source.npz',allow_pickle=False) as data:
            source=data['complex_reflectivity_grid']
        ri,ai=np.nonzero(source);tr=4.92+ri*cfg.dr+25.1;ta=-50+ai;coefs=source[ri,ai]
        scenes=1;snrs=PROTOCOL['dense']['noise_SNR_db'];reps=30
        quant=float(np.max(abs(cfg.angles[np.rint((np.sin(np.deg2rad(ta))-cfg.u[0])/(cfg.u[1]-cfg.u[0])).astype(int)]-ta)))
        dump(out/'source_scene_identity.json',{'saved_source_shape':list(source.shape),'generated_reflectors':1000,'requested_layers':50,
             'original_coherent_occupied_pixels':len(ri),'coherent_power_before_mapping':float(np.sum(abs(coefs)**2)),
             'source_angle_quantization_max_deg':quant,'old_to_new_actual_coordinate_shift_m':25.1,'source_recipe_unchanged':True})
    else:
        cfg=nm.Config(range_origin=42,range_span=80,angular_nodes=201)
        source=np.loadtxt(SHIP_ARRIVALS,delimiter=',',skiprows=1,usecols=(2,3,6))
        keep=(source[:,0]>=42)&(source[:,0]<=122)&(abs(source[:,1])<=50)&(source[:,2]>0)
        tr,ta,amp=source[keep].T;scenes=30;snrs=[10];reps=1
        dump(out/'source_scene_identity.json',{'unchanged_arrivals':len(tr),'range_min_m':float(tr.min()),'range_max_m':float(tr.max()),
            'no_coordinate_translation':True,'amplitude':'source MetricAmplitude unchanged','origin_view':'original x-z receiver plane, displayed using physical polar coordinates'})
    chains={m:GridChain(nm,cfg,m) for m in METHODS};nf=chains['dense'].nf
    if len({ch.nf for ch in chains.values()})!=1:raise RuntimeError('Methods do not have a shared noise record length')
    for scene_rep in range(scenes):
        if kind=='shipwreck':
            phase_seed=20250821+scene_rep;coefs=amp*np.exp(1j*np.random.default_rng(phase_seed).uniform(-np.pi,np.pi,len(amp)))
        else:phase_seed=None
        truth=nm.grid_catalog(cfg,tr,ta,coefs)
        raws={m:ch.raw(truth) for m,ch in chains.items()};clean={m:chains[m].reconstruct(raws[m]) for m in METHODS}
        powers={m:float(np.mean(abs(fft.ifft(raws[m],axis=0))**2)) for m in METHODS}
        for m in METHODS:cleanrows.append({'method':m,'phase_repetition':scene_rep,'phase_seed':phase_seed,**evaluate(nm,truth,clean[m],cfg,clean['dense'])})
        scene_meta.append({'phase_repetition':scene_rep,'phase_seed':phase_seed,'raw_signal_power_by_method':powers,
            'occupied_range_cells':int(np.any(truth!=0,axis=1).sum()),'occupied_pixels':int(np.count_nonzero(truth))})
        valid=(cfg.ranges>=cfg.range_origin)&(cfg.ranges<=cfg.range_origin+cfg.range_span)
        if scene_rep==0:
            np.savez_compressed(out/'representative_clean_images.npz',ranges_m=cfg.ranges[valid],angles_deg=cfg.angles,truth_complex=truth[valid],
                catalogue_ranges_m=tr,catalogue_angles_deg=ta,catalogue_coefficients=coefs,**{m:clean[m][valid] for m in METHODS})
            panel(out,'clean_common_scale',cfg.ranges[valid],cfg.angles,{m:clean[m][valid] for m in METHODS},True,kind=='shipwreck')
            panel(out,'clean_per_image',cfg.ranges[valid],cfg.angles,{m:clean[m][valid] for m in METHODS},False,kind=='shipwreck')
        for noise_rep in range(reps):
            seed=(20261016100+noise_rep) if kind=='dense1000x50' else (20261017100+scene_rep)
            gen=np.random.default_rng(seed);noise=(gen.normal(size=(nf,128))+1j*gen.normal(size=(nf,128)))/np.sqrt(2)
            noiseFFT=fft.fft(noise,axis=0)
            noise_images={m:chains[m].reconstruct(noiseFFT[:,:chains[m].n]) for m in METHODS}
            for rule in PROTOCOL['noise_rules']:
                for snr in snrs:
                    representative={}
                    for m in METHODS:
                        sigma=np.sqrt((powers['dense'] if rule=='common_per_element_noise' else powers[m])/10**(snr/10))
                        image=clean[m]+sigma*noise_images[m]
                        rows.append({'method':m,'noise_normalization':rule,'nominal_snr_db':snr,'phase_repetition':scene_rep,'phase_seed':phase_seed,
                            'noise_repetition':noise_rep,'noise_seed':seed,'pre_frontend_noise_sigma':float(sigma),
                            'actual_raw_channel_snr_db':float(10*np.log10(powers[m]/sigma**2)),**evaluate(nm,truth,image,cfg,clean['dense'],clean[m])})
                        representative[m]=image[valid]
                    if scene_rep==noise_rep==0:
                        tag=f'noisy_{rule}_{snr}dB'
                        np.savez_compressed(out/(tag+'.npz'),ranges_m=cfg.ranges[valid],angles_deg=cfg.angles,**representative)
                        if snr==snrs[0] or kind=='shipwreck':panel(out,tag+'_common',cfg.ranges[valid],cfg.angles,representative,True,kind=='shipwreck')
            table(out/'all_trials.csv',rows)
            if kind=='dense1000x50':print(f'DENSE noise {noise_rep+1}/{reps} completed',flush=True)
        print(f'{kind} phase {scene_rep+1}/{scenes} completed',flush=True)
        table(out/'clean_trials.csv',cleanrows);dump(out/'scene_metadata.json',scene_meta)
    agg=aggregate(rows,['noise_normalization','nominal_snr_db','method']);table(out/'summary_mean_ci.csv',agg)
    table(out/'clean_summary_mean_ci.csv',aggregate(cleanrows,['method']))
    dump(out/'model_metadata.json',{'configuration':asdict(cfg),'FFT_samples':nf,'energy_metadata':{m:chains[m].meta for m in METHODS},
        'noise_free_Dense_reference':True,'separate_geometry_vs_noise_scope':PROTOCOL['dense' if kind=='dense1000x50' else 'ship']['repeat_identity'],
        'full_precision_metrics':'float64/complex128; no display gamma/TVG/smoothing','guard_background_fraction_rule':.05})
    if kind=='dense1000x50':
        for rule in PROTOCOL['noise_rules']:curves(out,'noise_'+rule,agg,rule)
    return {'rows':len(rows),'clean_rows':len(cleanrows),'shared_FFT_samples':nf,'scene_metadata':scene_meta,'shape':[int(valid.sum()),201]}

def run():
    global HERE
    p=argparse.ArgumentParser()
    p.add_argument('--suite',choices=['gates','point','dense','ship','all'],default='gates')
    p.add_argument('--out',type=Path,default=HERE)
    args=p.parse_args()
    HERE=args.out.resolve()
    if HERE.is_relative_to(PROJECT/'data') or HERE.is_relative_to(PROJECT/'code'):
        p.error('Choose an output directory outside the bundled data and code.')
    if args.suite=='ship' and not SHIP_ARRIVALS.is_file():
        p.error('Extract shipwreck-data_v10.zip at the repository root before running the shipwreck suite.')
    HERE.mkdir(parents=True,exist_ok=True)
    nm,sn,hashes=freeze(args.suite);started=time.time();gates(nm)
    requested={'point':point,'dense':lambda n,s:grid_suite(n,'dense1000x50'),'ship':lambda n,s:grid_suite(n,'shipwreck')}
    selected=['point','dense'] if args.suite=='all' else ([] if args.suite=='gates' else [args.suite])
    for name in selected:
        start=time.time();print('START',name,flush=True)
        try:
            result=requested[name](nm,sn);result.update(status='computed_pending_author_validation',elapsed_seconds=time.time()-start)
            dump(HERE/('receipt_'+name+'.json'),result);print('FINISH',name,result['elapsed_seconds'],flush=True)
        except Exception:
            error=traceback.format_exc();dump(HERE/('receipt_'+name+'.json'),{'status':'failed_preserved','traceback':error,'elapsed_seconds':time.time()-start});raise
    unchanged={p:sha(p)==h for p,h in hashes.items()}
    if not all(unchanged.values()):raise RuntimeError('Frozen source changed while running')
    dump(HERE/('run_receipt_'+args.suite+'.json'),{'selected':selected,'all_requested_computed':True,'elapsed_seconds':time.time()-started,'all_protected_hashes_unchanged':all(unchanged.values()),'sources':unchanged})

if __name__=='__main__':run()
