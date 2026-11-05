"""Uniform-protocol scene/density/noise tests; old experiment files are read-only.

The central numerical model owns Tx normalization and the common Tx/Rx FIR.
This module only creates controlled scenes, evaluates declared physical ROIs,
and saves all replications for later author inspection. No test is tuned to
produce a particular ordering. Run only after core_ready.json is present.
"""
from __future__ import annotations
import argparse, csv, hashlib, json, time, traceback, shutil
from collections import defaultdict
from dataclasses import asdict
from datetime import datetime
from functools import lru_cache
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from scipy import ndimage, optimize, stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
PROJECT = ROOT.parent
METHODS = ('dense', 'cdm', 'fdm')
LABELS = {'dense': 'Dense ULA', 'cdm': 'CDM SW-ULA', 'fdm': 'packed FDM SW-ULA'}
COLORS = {'dense': '#0072B2', 'cdm': '#D55E00', 'fdm': '#009E73'}
STYLES = {'dense': ('o', '-'), 'cdm': ('s', '--'), 'fdm': ('^', '-.')}
PROTOCOL = {
    'version': 'uniform-scene-noise-20261010-v1',
    'core': 'locked common/numerical_model.py; require core_ready.json',
    'primary_source_energy': 'equal sum of all actual simultaneous Tx-element excitation energies',
    'waveforms': 'core generation and common 129-tap Tx FIR centre crop; common continuous-record Rx FIR',
    'large_scene_range_origin_m': 30, 'large_scene_range_span_m': 2.4,
    'large_scene_angular_nodes': 201, 'FOV_deg': [-50, 50],
    'guard_range_m': .06, 'guard_angle_deg': 1.,
    'guard_sensitivity_ranges_m': [.03, .06, .12],
    'guard_sensitivity_angles_deg': [.5, 1., 2.],
    'fixed_background': 'fixed physical ROI AND 47 <= abs(angle) <= 50 deg; never trim using an image',
    'fixed_background_reference': 'mean background squared magnitude divided by mean Dense squared magnitude on the same evaluation ROI',
    'fixed_background_validity': 'report truth occupancy; source-free background interpretation unavailable if any truth cell is in fixed band',
    'density': {
        'scatter_count': [1, 10, 50, 200, 500, 1000],
        'range_layers': [1, 5, 10, 20, 50],
        'delay_extent_m': [.15, .3, .6, 1.2, 2.4],
        'dynamic_range_db': [0, 10, 20, 40],
        'additional_stress': '1000 scatterers, 50 requested layers, 0.8m extent, 40dB amplitude range',
        'independent_repetitions': 3, 'seed_base': 2026101000,
        'stress_seed_base': 2026104000,
    },
    'scene_ordering': {'synthetic_seed_base': 2026101901, 'ship_phase_seed_base': 2026102001,
                       'ship_coordinates': 'unchanged source RangeM/AngleDeg; no translation to 30m',
                       'ship_amplitude': 'source MetricAmplitude; no extra TVG/gamma or smoothing'},
    'scene_diversity': {'shapes': ['random', 'wall', 'two_clusters'], 'count': 500,
                        'layers': 50, 'extent_m': .8, 'seed_base': 2026103000, 'repetitions': 3},
    'noise': {'truth_phase_seed': 2026105001, 'calibration_seed_base': 2026105100,
              'calibration_blocks': 60, 'independent_test_blocks': 30,
              'FFT_record_samples': 8192, 'range_step_m': .01, 'angle_step_deg': .5,
              'SNR_db': [-35, -25, -15, -5],
              'noise_rules': ['common_per_element_noise', 'equal_raw_channel_snr'],
              'noise_location': 'complex raw receiver samples BEFORE common Rx FIR',
              'prefix': 'same standardized 128-channel random draw; sparse methods use first 32 channels',
              'threshold': '95th percentile (higher interpolation) of 60 unit-noise block maxima, separately by method',
              'detector': 'full physical ROI local maxima -> physical NMS -> one-to-one truth assignment',
              'NMS_range_m': .03, 'NMS_angle_deg': 2., 'match_range_m': .02, 'match_angle_deg': 1.,
              'fixed_truth': '18 grid-aligned targets: ranges30.2/30.4/30.6m, angles-40/-24/-8/8/24/40deg',
              'independence': 'same fixed geometry, 30 independent noise blocks; two rules and SNRs reuse these blocks'},
    'saved_arrays': 'truth complex64 and reconstructed magnitudes float32; all metrics computed before conversion',
}

def _dump(path, data):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False,
                    default=lambda x: x.item() if isinstance(x, np.generic) else str(x)) + '\n', encoding='utf-8')

def _csv(path, rows):
    rows = list(rows)
    if not rows: return
    keys = list(dict.fromkeys(k for row in rows for k in row))
    with Path(path).open('w', encoding='utf-8-sig', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=keys); writer.writeheader(); writer.writerows(rows)

def _sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()

@lru_cache(maxsize=1)
def _core():
    if not (ROOT/'core_ready.json').is_file():
        raise RuntimeError('Locked numerical core is not ready: core_ready.json is absent')
    gate=json.loads((ROOT/'core_ready.json').read_text(encoding='utf-8'))
    if gate.get('gate')!='passed' or gate.get('model_sha256')!=_sha(ROOT/'common'/'numerical_model.py'):
        raise RuntimeError('Core gate/hash does not match the current numerical model')
    import numerical_model as model
    return model

def _savefig(fig, path):
    fig.tight_layout()
    for ext in ['png', 'svg', 'pdf']:
        fig.savefig(str(path)+'.'+ext, dpi=300, bbox_inches='tight')
    plt.close(fig)

def _valid(cfg):
    return np.broadcast_to(((cfg.ranges >= cfg.range_origin) &
                            (cfg.ranges <= cfg.range_origin+cfg.range_span))[:, None],
                            (len(cfg.ranges), len(cfg.angles)))

def _fixed_band(cfg):
    return _valid(cfg) & np.broadcast_to(((abs(cfg.angles) >= 47.-1e-10) &
                                         (abs(cfg.angles) <= 50.+1e-10))[None, :],
                                         (len(cfg.ranges), len(cfg.angles)))

def _evaluate(truth, image, cfg, dense, guard_m=.06, guard_deg=1.):
    fixed = _fixed_band(cfg)
    q = _core().metrics(truth, image, cfg, dense=dense, guard_m=guard_m,
                         guard_deg=guard_deg, fixed_background=fixed)
    roi = _valid(cfg); p = abs(image); ref_power = float(np.mean(abs(dense)[roi]**2))
    has_truth = bool(np.any((abs(truth) > 0) & fixed)); count = int(fixed.sum())
    mean_power = float(np.mean(p[fixed]**2)) if count else None
    ratio = mean_power/ref_power if count and ref_power > 0 and not has_truth else None
    q.update({'evaluation_roi_pixels': int(roi.sum()),
              'fixed_background_source_free': not has_truth,
              'fixed_background_reference_dense_roi_mean_power': ref_power,
              'fixed_background_to_dense_mean_power_ratio': ratio,
              'fixed_background_to_dense_mean_power_db': float(10*np.log10(ratio)) if ratio and ratio > 0 else None,
              'fixed_background_definition': '47<=abs(angle)<=50deg in fixed physical range ROI',
              'actual_guard_definition': 'coordinate angular guard, discrete radial rows from requested physical halfwidth'})
    return q

def _save_scene(out, cfg, truth, images, tag, catalogue=None):
    arrays = {'ranges_m': cfg.ranges, 'angles_deg': cfg.angles,
              'truth_complex': truth.astype(np.complex64),
              'config_json': np.array(json.dumps(asdict(cfg))),
              **{m: abs(im).astype(np.float32) for m, im in images.items()}}
    if catalogue:
        arrays.update({k: np.asarray(catalogue[k]) for k in ['ranges_m', 'angles_deg', 'coefficients']
                       if k in catalogue})
        # Keep continuous catalogue coordinates distinct from image axes.
        arrays['ranges_m'] = cfg.ranges; arrays['angles_deg'] = cfg.angles
        arrays['catalogue_ranges_m'] = np.asarray(catalogue.get('ranges_m', []))
        arrays['catalogue_angles_deg'] = np.asarray(catalogue.get('angles_deg', []))
    np.savez_compressed(out/(tag+'.npz'), **arrays)
    _dump(out/(tag+'_configuration.json'), {'config': asdict(cfg), 'saved_precision': PROTOCOL['saved_arrays']})

def _image_panels(out, cfg, images, tag, common=False):
    valid = (cfg.ranges >= cfg.range_origin) & (cfg.ranges <= cfg.range_origin+cfg.range_span)
    reference = max(float(abs(im[valid]).max()) for im in images.values()) if common else None
    fig, axes = plt.subplots(1, 3, figsize=(11.7, 3.65))
    for ax, method in zip(axes, METHODS):
        mag = abs(images[method][valid]); ref = float(mag.max()) if reference is None else reference
        z = np.maximum(20*np.log10(np.maximum(mag, 1e-30)/max(ref, 1e-30)), -50)
        pcm = ax.pcolormesh(cfg.angles, cfg.ranges[valid], z, cmap='viridis', vmin=-50, vmax=0,
                            shading='auto', rasterized=True)
        ax.set(xlabel='Angle (deg)', ylabel='Range (m)', title=LABELS[method])
        fig.colorbar(pcm, ax=ax, label='Amplitude (dB)')
    _savefig(fig, out/tag)

def _aggregate(rows, group_keys, metric_keys, repeat_label):
    grouped = defaultdict(list)
    for row in rows: grouped[tuple(row[k] for k in group_keys)].append(row)
    result = []
    for key, group in grouped.items():
        item = dict(zip(group_keys, key)); item['replications'] = len(group); item['repeat_scope'] = repeat_label
        for metric in metric_keys:
            a = np.array([r[metric] for r in group if r.get(metric) is not None], dtype=float)
            item[metric] = float(a.mean()) if len(a) else None
            item[metric+'_sd'] = float(a.std(ddof=1)) if len(a) > 1 else None
            item[metric+'_ci95_halfwidth'] = float(stats.t.ppf(.975, len(a)-1)*a.std(ddof=1)/np.sqrt(len(a))) if len(a)>1 else None
            item[metric+'_available_count'] = len(a)
        result.append(item)
    return result

SCENE_METRICS = ['truth_magnitude_nmse_db', 'truth_magnitude_correlation', 'dense_magnitude_correlation',
                 'support_leakage', 'cnr_linear_magnitude', 'background_fraction',
                 'fixed_background_to_dense_mean_power_ratio', 'fixed_background_mean_power']

def _curves(rows, xkey, keys, path, xlabel, repeat_text):
    fig, axes = plt.subplots(1, len(keys), figsize=(4.0*len(keys), 3.3), squeeze=False)
    for ax, (metric, ylabel) in zip(axes[0], keys):
        for method in METHODS:
            q = sorted([r for r in rows if r['method'] == method], key=lambda r: float(r[xkey]))
            x = [float(r[xkey]) for r in q]; y = [np.nan if r.get(metric) is None else r[metric] for r in q]
            err = [0 if r.get(metric+'_ci95_halfwidth') is None else r[metric+'_ci95_halfwidth'] for r in q]
            marker, line = STYLES[method]
            ax.errorbar(x, y, yerr=err, color=COLORS[method], marker=marker, linestyle=line,
                         linewidth=1.3, capsize=2.5, label=LABELS[method])
        ax.set(xlabel=xlabel, ylabel=ylabel); ax.grid(alpha=.18); ax.legend(fontsize=8)
    _savefig(fig, path)
    _dump(Path(str(path)+'_caption.json'), {'repeat': repeat_text, 'errorbars': 'Student-t 95% intervals of the stated repetitions; unavailable metrics omitted', 'dual_encoding': 'colour+marker+linestyle'})

def _guard_rows(truth, images, cfg, context):
    result = []
    for gm in PROTOCOL['guard_sensitivity_ranges_m']:
        for ga in PROTOCOL['guard_sensitivity_angles_deg']:
            for method in METHODS:
                result.append({'scene_context': context, 'method': method,
                               **_evaluate(truth, images[method], cfg, images['dense'], gm, ga)})
    return result

def density(out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True);model=_core();cfg=model.Config(range_span=2.4,angular_nodes=201)
    conditions=[('scatter_count',n,n,50,.8,0) for n in PROTOCOL['density']['scatter_count']]
    conditions += [('range_layers',l,1000,l,.8,0) for l in PROTOCOL['density']['range_layers']]
    conditions += [('delay_extent_m',e,500,50,e,0) for e in PROTOCOL['density']['delay_extent_m']]
    conditions += [('dynamic_range_db',d,500,50,.8,d) for d in PROTOCOL['density']['dynamic_range_db']]
    conditions += [('stress_1000_50_40dB',40,1000,50,.8,40)]
    rows=[]; catalogue=[]; guards=[]
    for i,(factor,value,n,layers,extent,dynamic) in enumerate(conditions):
        for rep in range(3):
            seed=(2026104000+rep) if factor.startswith('stress') else 2026101000+100*i+rep
            truth,cat=model.scene_grid(cfg,n,layers,extent,dynamic,seed)
            images={m:model.image_grid(cfg,truth,m)[0] for m in METHODS}
            for method in METHODS:
                rows.append({'factor':factor,'factor_value':value,'scatter_count':n,'requested_layers':layers,
                             'extent_m':extent,'dynamic_range_db':dynamic,'scene_realization':rep,'scene_seed':seed,
                             'method':method,**_evaluate(truth,images[method],cfg,images['dense'])})
            catalogue.append({'factor':factor,'factor_value':value,'scene_seed':seed,'individual_reflector_power':float(np.sum(abs(cat['coefficients'])**2)),
                              'grid_coherent_power':float(np.sum(abs(truth)**2)),'occupied_pixels':int(np.count_nonzero(truth)),
                              'occupied_range_cells':int(np.any(truth!=0,axis=1).sum()),'fixed_background_truth_pixels':int(((abs(truth)>0)&_fixed_band(cfg)).sum())})
            tag=f'{factor}_{value}_rep{rep}';_save_scene(out,cfg,truth,images,tag,cat)
            if rep==0:
                _image_panels(out,cfg,images,tag)
                if (factor=='scatter_count' and n==1000) or factor.startswith('stress'):
                    _image_panels(out,cfg,images,tag+'_common',common=True)
                    guards.extend(_guard_rows(truth,images,cfg,tag))
            _csv(out/'metrics.csv',rows);_csv(out/'scene_catalogue.csv',catalogue)
            print(f'density {i+1}/{len(conditions)} {factor}={value} rep={rep+1}/3 completed',flush=True)
        _dump(out/'progress.json',{'conditions_completed':i+1,'conditions_total':len(conditions),'scene_realizations_completed':3*(i+1),'metric_rows':len(rows)})
    aggregate=_aggregate(rows,['factor','factor_value','method'],SCENE_METRICS,'3 independent generated scenes per configuration')
    _csv(out/'aggregate_metrics.csv',aggregate);_csv(out/'guard_sensitivity.csv',guards)
    for factor in list(PROTOCOL['density'])[:4]:
        q=[r for r in aggregate if r['factor']==factor]
        _curves(q,'factor_value',[('dense_magnitude_correlation','Correlation to Dense magnitude'),('support_leakage','Guard-outside power fraction'),
                                 ('fixed_background_to_dense_mean_power_ratio','Fixed-band / Dense-ROI mean power')],out/factor,factor,'3 generated scenes')
    _dump(out/'scope.json',{'primary_configurations':20,'additional_1000_50_40dB_configurations':1,'independent_scenes_per_configuration':3,
                          'metric_rows_expected':189,'metric_ROI_m':[30,32.4],'source_energy':'actual Tx-element total excitation energy=1',
                          'guard':'strict coordinate angle +/-1deg; distance +/-0.06m on matching radial grid',
                          'fixed_background':'unchanged outer angular band47..50deg; source-free checked for every scene',
                          'causal_warning':'guard-union leakage changes with protected area; fixed-band and CNR retained, no expected monotonic ranking',
                          'saved':PROTOCOL['saved_arrays']})
    return rows

@lru_cache(maxsize=1)
def _ship_data():
    path=PROJECT/'paperA_revise'/'11_A11_shipwreck_common_scene'/'results'/'A11_raytrace_arrivals.csv'
    data=np.loadtxt(path,delimiter=',',skiprows=1,usecols=(2,3,6));r,a,amp=data.T
    valid=(r>=42)&(r<=122)&(a>=-50)&(a<=50)&(amp>0)
    return path,r[valid],a[valid],amp[valid]

def scene_ordering(out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True);model=_core();rows=[];provenance=[];guards=[]
    for kind in ['synthetic1000x50','legacy_ship']:
        for rep in range(3):
            if kind=='legacy_ship':
                path,r,a,amp=_ship_data();seed=2026102001+rep;cfg=model.Config(range_origin=42,range_span=80,angular_nodes=201)
                coeff=amp*np.exp(1j*np.random.default_rng(seed).uniform(-np.pi,np.pi,len(r)))/np.sqrt(np.sum(amp**2))
                truth=model.grid_catalog(cfg,r,a,coeff);cat={'ranges_m':r,'angles_deg':a,'coefficients':coeff}
                meta={'source_csv':str(path),'source_sha256':_sha(path),'arrivals':len(r),'minimum_range_m':float(r.min()),'maximum_range_m':float(r.max()),
                      'geometry':'unchanged old RangeM/AngleDeg, orientation unchanged','repeat':'same geometry, independently generated phases',
                      'amplitude':'MetricAmplitude, no extra TVG/gamma/visual smoothing'}
            else:
                seed=2026101901+rep;cfg=model.Config(range_span=2.4,angular_nodes=201);truth,cat=model.scene_grid(cfg,1000,50,.8,0,seed)
                meta={'repeat':'independently generated geometry and phase','normalization':cat['amplitude_norm']}
            images={m:model.image_grid(cfg,truth,m)[0] for m in METHODS}
            for method in METHODS:rows.append({'scene':kind,'replication':rep,'seed':seed,'method':method,**_evaluate(truth,images[method],cfg,images['dense'])})
            tag=f'{kind}_rep{rep}';_save_scene(out,cfg,truth,images,tag,cat)
            if rep==0:
                _image_panels(out,cfg,images,tag);_image_panels(out,cfg,images,tag+'_common',common=True)
                guards.extend(_guard_rows(truth,images,cfg,tag))
            provenance.append({'scene':kind,'replication':rep,'seed':seed,**meta})
            _csv(out/'metrics.csv',rows);_dump(out/'scene_provenance.json',provenance)
            print(f'scene_ordering {kind} rep={rep+1}/3 completed',flush=True)
    _csv(out/'aggregate_metrics.csv',_aggregate(rows,['scene','method'],SCENE_METRICS,'synthetic:3 generated geometries; ship:3 phases of one geometry'))
    _csv(out/'guard_sensitivity.csv',guards)
    _dump(out/'scope.json',{'source_model':'uniform common-FIR common-energy broadband far-field angular-grid reconstruction',
                          'scene_ordering':'not predetermined; direct truth and Dense-reference metrics separately retained',
                          'ship_evidence':'geometry-derived extended-scatter stress test, not acoustic acquisition',
                          'fixed_background':'outer angular band may contain ship/seabed truth; invalid as source-free background in that case',
                          'saved':PROTOCOL['saved_arrays']})
    return rows

def scene_diversity(out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True);model=_core();cfg=model.Config(range_span=2.4,angular_nodes=201);rows=[];guards=[]
    for kind in PROTOCOL['scene_diversity']['shapes']:
        for rep in range(3):
            seed=2026103000+rep;truth,cat=model.scene_grid(cfg,500,50,.8,0,seed,shape=kind)
            images={m:model.image_grid(cfg,truth,m)[0] for m in METHODS}
            for method in METHODS:rows.append({'independent_structure':kind,'scene_realization':rep,'seed':seed,'method':method,**_evaluate(truth,images[method],cfg,images['dense'])})
            tag=f'{kind}_rep{rep}';_save_scene(out,cfg,truth,images,tag,cat)
            if rep==0:
                _image_panels(out,cfg,images,tag);guards.extend(_guard_rows(truth,images,cfg,tag))
            _csv(out/'metrics.csv',rows);print(f'scene_diversity {kind} rep={rep+1}/3 completed',flush=True)
    _csv(out/'aggregate_metrics.csv',_aggregate(rows,['independent_structure','method'],SCENE_METRICS,'3 independent random realizations within each procedural structure'))
    _csv(out/'guard_sensitivity.csv',guards)
    _dump(out/'scope.json',{'structures':['random spatial field','gently varying procedural wall','two clusters'],
                          'scatterers':500,'requested_layers':50,'individual_reflector_power':1,'repetitions':3,
                          'generators':'wall and two-cluster structures do not necessarily occupy 50 range cells; actual counts are in metrics',
                          'no_hardware_or_real_scene_claim':True})
    return rows

def guard_sensitivity(out, source_root=None):
    out=Path(out);out.mkdir(parents=True,exist_ok=True);base=out.parent if source_root is None else Path(source_root)
    rows=[]
    for suite in ['density','scene_ordering','scene_diversity']:
        path=base/suite/'guard_sensitivity.csv'
        if not path.is_file():raise RuntimeError(f'Guard sensitivity requires completed {path}')
        rows.extend(csv.DictReader(path.open(encoding='utf-8-sig')))
    _csv(out/'metrics.csv',rows)
    for context in dict.fromkeys(r['scene_context'] for r in rows):
        q=[r for r in rows if r['scene_context']==context and abs(float(r['guard_range_m'])-.06)<1e-10]
        for r in q:
            for k in ['support_leakage','cnr_linear_magnitude','background_fraction']:
                r[k]=float(r[k]) if r[k] else None
        _curves(q,'guard_angle_deg',[('support_leakage','Guard-outside power fraction'),('background_fraction','Available background fraction')],
                out/context,'Angular guard halfwidth (deg)','same saved representative reconstruction; no new independent replication')
    _dump(out/'scope.json',{'ranges_m':[.03,.06,.12],'angles_deg':[.5,1,2],'full_cross_product':True,
                          'input':'metrics computed on full-precision representative images while parent scene suites ran',
                          'not_new_independent_scenes':True,'no_selected_guard_by_outcome':True,
                          'NMSE_and_correlations':'do not depend on guard; repeated only for audit'})
    return rows

def _blind_detection(image, ranges, angles, threshold, truth_ranges, truth_angles):
    power=abs(image);dr=float(np.median(np.diff(ranges)));da=float(np.median(np.diff(angles)))
    sr=int(np.ceil(.03/dr-1e-9));sa=int(np.ceil(2./da-1e-9))
    local=ndimage.maximum_filter(power,size=(2*sr+1,2*sa+1),mode='constant',cval=-np.inf)
    candidates=np.argwhere((power>=threshold)&(power==local))
    if not len(candidates):return 0.,0,0
    order=np.argsort(power[candidates[:,0],candidates[:,1]])[::-1];selected=[]
    for index in order:
        ir,ia=candidates[index]
        if any(abs(ranges[ir]-ranges[jr])<=.03+1e-12 and abs(angles[ia]-angles[ja])<=2.+1e-12 for jr,ja in selected):continue
        selected.append((ir,ia))
    found=np.asarray(selected);r=ranges[found[:,0]];a=angles[found[:,1]]
    delta_r=abs(np.asarray(truth_ranges)[:,None]-r[None,:]);delta_a=abs(np.asarray(truth_angles)[:,None]-a[None,:])
    valid=(delta_r<=.02+1e-12)&(delta_a<=1.+1e-12)
    cost=np.sqrt((delta_r/.02)**2+delta_a**2);cost[~valid]=1e9
    it,ic=optimize.linear_sum_assignment(cost);matched=int(np.sum(valid[it,ic]))
    return matched/len(truth_ranges),len(selected)-matched,len(selected)

def common_noise(out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True);model=_core();cfg=model.Config();nf=8192
    rr=np.arange(30.,30.8001,.01);aa=np.arange(-50.,50.0001,.5)
    axis_cfg=SimpleNamespace(ranges=rr,angles=aa,dr=.01,range_origin=float(rr[0]),range_span=float(rr[-1]-rr[0]))
    tr=np.repeat([30.2,30.4,30.6],6);ta=np.tile([-40.,-24.,-8.,8.,24.,40.],3)
    coeff=np.exp(1j*np.random.default_rng(2026105001).uniform(-np.pi,np.pi,18))/np.sqrt(18)
    truth=np.zeros((len(rr),len(aa)),complex)
    for r,a,c in zip(tr,ta,coeff):truth[np.argmin(abs(rr-r)),np.argmin(abs(aa-a))]=c
    clean={};raws={};filters={};metadata={};signal_power={}
    for method in METHODS:
        raw,mf,meta=model.point_raw(cfg,method,tr,ta,coeff,geometry='farfield',nfft=nf)
        raws[method]=raw;filters[method]=mf;metadata[method]=meta
        clean[method]=model.reconstruct_points(cfg,method,raw,mf,rr,aa)
        signal_power[method]=float(np.mean(abs(np.fft.ifft(raw,axis=0))**2))
    cal=[];noise_bank={m:[] for m in METHODS}
    for rep in range(90):
        gen=np.random.default_rng(2026105100+rep)
        noise=(gen.standard_normal((nf,128))+1j*gen.standard_normal((nf,128)))/np.sqrt(2)
        nf_noise=np.fft.fft(noise,axis=0)
        for method in METHODS:
            nr=128 if method=='dense' else 32
            image=model.reconstruct_points(cfg,method,nf_noise[:,:nr],filters[method],rr,aa)
            if rep<60:cal.append({'method':method,'calibration_repetition':rep,'noise_seed':2026105100+rep,'pure_noise_max_amplitude':float(abs(image).max())})
            else:noise_bank[method].append(image.astype(np.complex64))
        if rep%5==0 or rep==89:print(f'common_noise pure-noise block {rep+1}/90 completed',flush=True)
    threshold={m:float(np.quantile([r['pure_noise_max_amplitude'] for r in cal if r['method']==m],.95,method='higher')) for m in METHODS}
    _csv(out/'pure_noise_calibration.csv',cal);rows=[];null=[];null_summary=[]
    for method in METHODS:
        count=0
        for rep,noise in enumerate(noise_bank[method]):
            hit=bool(abs(noise).max()>threshold[method]);count+=int(hit)
            null.append({'method':method,'noise_repetition':rep,'noise_seed':2026105160+rep,'pure_noise_exceeded_block_threshold':hit,
                         'unit_noise_threshold':threshold[method],'pure_noise_max_amplitude':float(abs(noise).max())})
        lo=0. if count==0 else float(stats.beta.ppf(.025,count,31-count))
        hi=1. if count==30 else float(stats.beta.ppf(.975,count+1,30-count))
        null_summary.append({'method':method,'test_blocks':30,'exceeded_blocks':count,'empirical_block_false_alarm_probability':count/30,
                             'conditional_binomial_CI95_lower':lo,'conditional_binomial_CI95_upper':hi,
                             'CI_scope':'independent test blocks conditional on the one calibrated threshold; threshold calibration uncertainty excluded'})
    for rule in PROTOCOL['noise']['noise_rules']:
        for snr in PROTOCOL['noise']['SNR_db']:
            for method in METHODS:
                sigma=np.sqrt((signal_power['dense'] if rule=='common_per_element_noise' else signal_power[method])/10**(snr/10))
                actual_snr=float(10*np.log10(signal_power[method]/sigma**2))
                for rep,noise in enumerate(noise_bank[method]):
                    image=clean[method]+sigma*noise;dense_reference=clean['dense']+sigma*noise_bank['dense'][rep]
                    # Dense reference is a same-block, same physical-sigma comparison for primary; the supplementary
                    # equal-SNR branch uses its own Dense sigma, not the sparse method's sigma.
                    if rule=='equal_raw_channel_snr':
                        dense_sigma=np.sqrt(signal_power['dense']/10**(snr/10));dense_reference=clean['dense']+dense_sigma*noise_bank['dense'][rep]
                    pd,fp,npeaks=_blind_detection(image,rr,aa,threshold[method]*sigma,tr,ta)
                    q=_evaluate(truth,image,axis_cfg,dense_reference)
                    mag=abs(image);nom=abs(clean[method]);fit=float(np.vdot(nom,mag).real/max(np.vdot(nom,nom).real,1e-30))
                    perturb=float(10*np.log10(max(float(np.sum((mag-fit*nom)**2)),1e-30)/max(float(np.sum((fit*nom)**2)),1e-30)))
                    rows.append({'method':method,'source_energy_mode':'tx_element_unit_energy','noise_normalization':rule,
                                 'nominal_dense_or_equal_input_snr_db':snr,'noise_repetition':rep,'noise_seed':2026105160+rep,
                                 'pre_frontend_noise_sigma':float(sigma),'actual_pre_frontend_raw_channel_snr_db':actual_snr,
                                 'calibrated_amplitude_threshold':threshold[method]*sigma,'detection_probability':pd,'false_detections':fp,
                                 'total_blind_detections':npeaks,'actual_total_tx_element_excitation_energy':metadata[method]['final_tx_element_excitation_energy'],
                                 'noise_induced_magnitude_NMSE_db':perturb,'fitted_clean_scale':fit,**q})
            print(f'common_noise rule={rule} SNR={snr} all3x30 cases completed',flush=True)
            _csv(out/'metrics.csv',rows)
    _csv(out/'independent_false_alarm_test.csv',null);_csv(out/'false_alarm_summary.csv',null_summary)
    aggregates=_aggregate(rows,['noise_normalization','nominal_dense_or_equal_input_snr_db','method'],
                          ['detection_probability','false_detections','noise_induced_magnitude_NMSE_db',*SCENE_METRICS],
                          '30 noise realizations of one fixed18-target geometry, reused across SNR/noise-rule branches')
    _csv(out/'aggregate_metrics.csv',aggregates)
    for rule in PROTOCOL['noise']['noise_rules']:
        _curves([r for r in aggregates if r['noise_normalization']==rule],'nominal_dense_or_equal_input_snr_db',
                [('detection_probability','Detection probability'),('false_detections','False peaks / image'),
                 ('noise_induced_magnitude_NMSE_db','Perturbation magnitude NMSE (dB)')],out/rule,
                'Nominal Dense / equal raw input SNR (dB)','30 noise realizations of one fixed scene')
    np.savez_compressed(out/'noise_realization_bank.npz',ranges_m=rr,angles_deg=aa,truth_complex=truth.astype(np.complex64),
                        truth_ranges_m=tr,truth_angles_deg=ta,coefficients=coeff,
                        **{m+'_unit_noise':np.asarray(noise_bank[m]) for m in METHODS},
                        **{m+'_clean_complex':clean[m].astype(np.complex64) for m in METHODS},
                        **{m+'_clean_raw_fft':raws[m].astype(np.complex64) for m in METHODS})
    _image_panels(out,axis_cfg,clean,'clean_18_targets',common=True)
    _dump(out/'noise_metadata.json',{'source_energy_by_method':metadata,'raw_pre_frontend_signal_power_by_method':signal_power,
                                   'threshold_unit_noise_by_method':threshold,'axes_config':{'ranges_m':[30,30.8],'range_step_m':.01,'angle_step_deg':.5}})
    _dump(out/'scope.json',{'metric_rows_expected':720,'primary_rule':'common_per_element_noise','supplementary_rule':'equal_raw_channel_snr',
                          'noise':'CN(0,sigma^2) before continuous-record common Rx FIR, shared128channel standardized prefix',
                          'threshold':'60 independent calibration blocks before looking at30test blocks; method-specific block-max95percentile',
                          'false_alarm':'image-block event, NOT per-pixel false-alarm probability; nominal5% not a guarantee',
                          'source_energy':'actual Tx-element total excitation energy1 in both rules',
                          'independence':'one geometry,30independent noise blocks; branches reuse blocks',
                          'Dense_reference_correlation':'same paired-noise Dense reconstruction, not direct scatter truth',
                          'perturbation_NMSE':'fit clean magnitude to noisy magnitude; denominator fitted clean energy; clean reference same method',
                          'all_metric_domains':'identical fixed physical range-angle evaluation ROI',
                          'no_TVG_gamma_visual_metrics':True})
    return rows

SUITES={'density':density,'scene_ordering':scene_ordering,'scene_diversity':scene_diversity,
        'guard_sensitivity':guard_sensitivity,'common_noise':common_noise}

def run_suite(name, out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    start=time.time();core_path=ROOT/'common'/'numerical_model.py';before=_sha(core_path)
    _dump(out/'protocol_before_run.json',PROTOCOL)
    snapshots=out/'code_snapshot';snapshots.mkdir(exist_ok=True)
    shutil.copy2(__file__,snapshots/'scene_noise.py');shutil.copy2(core_path,snapshots/'numerical_model.py')
    try:
        rows=SUITES[name](out)
        if _sha(core_path)!=before:raise RuntimeError('Core source changed during scene/noise run')
        receipt={'suite':name,'status':'computed_pending_author_validation','elapsed_seconds':time.time()-start,
                 'metric_rows':len(rows),'core_sha256':before,'scene_noise_sha256':_sha(__file__),
                 'artifact_files':[str(p.relative_to(out)) for p in out.rglob('*') if p.is_file()]}
        _dump(out/'receipt.json',receipt)
        return out
    except Exception:
        error=traceback.format_exc();(out/'failure.txt').write_text(error,encoding='utf-8')
        _dump(out/'receipt.json',{'suite':name,'status':'failed_preserved_diagnostics','elapsed_seconds':time.time()-start,'traceback':error})
        raise

def main():
    p=argparse.ArgumentParser();p.add_argument('--suite',choices=['density'],default='density');p.add_argument('--out',type=Path)
    args=p.parse_args();_core()
    release=ROOT.parents[1]
    out=(args.out or release/'outputs'/'density'/('run_'+datetime.now().strftime('%Y%m%d_%H%M%S_%f'))).resolve()
    if out.is_relative_to(release/'data') or out.is_relative_to(release/'code'):
        p.error('Choose an output directory outside the bundled data and code.')
    out.mkdir(parents=True,exist_ok=True);_dump(out/'protocol_before_run.json',PROTOCOL)
    if args.suite=='all':
        for name in ['density','scene_ordering','scene_diversity','guard_sensitivity','common_noise']:
            run_suite(name,out/name)
    else:run_suite(args.suite,out)
    print('scene_noise output:',str(out),flush=True)

if __name__=='__main__':main()
