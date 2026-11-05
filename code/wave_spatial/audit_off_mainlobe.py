"""Frozen full-FOV coarse-to-fine off-mainlobe audit.

A virtual-array expansion evaluates the SAME complete two-way spectral sum.
Native8192-bin characteristic functions are queried via partial CZT, with
Gaussian virtual-element integration of the exact sinc or sinc² factor.
Independent native direct-frequency checks are performed before the audit.
"""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,sys,json
import numpy as np
from scipy import signal,special
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE));import run_wave_spatial as ws
nm=ws.nm;cfg=ws.CFG;out=ws.HERE/'full_fov_off_mainlobe';out.mkdir(parents=True,exist_ok=False)
scans=sorted(set(np.arange(-50,51,1).tolist()+list(cfg.focus)));coarse=np.round(np.arange(-50,50.00001,.25),8)
protocol={'frozen_utc':datetime.now(timezone.utc).isoformat(),'scans':scans,'coarse_target_grid_deg':.25,'refinement':'All coarse local maxima outside protectedabs(theta-scan)<=3deg, plus FOV endpoints and protection-boundary-adjacent samples; each candidate +/-0.5deg refined at0.01deg, unioned and clippedtoFOV; no favourable peak threshold or topK truncation','desired_normalization':'Desired local peak refined at0.01deg within +/-3deg, clippedtoFOV','off_mainlobe_metric':'Maximum amplitude at evaluated target bearings outsideabs(theta-scan)<=3deg relative to desiredlocal peak; not first-null PSLR and not nominal-order GL-only','spectral_model':'Same8192 complete frequency bins and same all-simultaneous Tx bank; explicit Tx-Rx virtual-aperture expansion and frequency-dependent elementdirectivity','element_integral':'Dense sinc=uniform rectangular aperture Gauss16; SW sinc²=triangular virtual aperture two8-point half-interval Gauss rules; validated separately over full250–650kHz RF extent and +/-50deg','characteristic_query':'Partial native-frequency sum via CZT on fixed[-80,+80]us; primary256 times fs interpolation,128 times diagnostic; physicalabsolute frequencies fc+nu incorporatedonce insideCZT, no carrierduplicated onqueries','frozen_validation':{'Gaussian_sinc_absolute_error':1e-12,'native_vs256_normalized_mag_L2':.001,'128_vs256_off_peak_level_difference_db':.05,'angles_for_native_check':[-50.,-41.25,-33.333333,-24.25,-8.,0.,8.,24.25,33.333333,41.25,50.],'scan_checks':list(cfg.focus)},'all_results_retained':True,'original_numerical_script_sha256':hashlib.sha256((HERE/'run_wave_spatial.py').read_bytes()).hexdigest(),'audit_script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
ws.dump(out/'protocol_before_run.json',protocol)

def element_nodes(method):
 if method=='dense':
  x,w=special.roots_legendre(16);return x/2,w/2
 x,w=special.roots_legendre(8);nodes=np.r_[(x-1)/2,(x+1)/2];weights=np.r_[w/2,w/2]*(1-abs(nodes));return nodes,weights

def coefficient_virtual_array(method,scan):
 send,mf,nu,_=nm.spectra(cfg,method,8192);f=cfg.fc+nu;q=0 if method=='dense' else ws.sector(scan);coef=send*nm.frontend_response(cfg,len(f))[:,None]*np.conj(mf[:,q,None])/len(f)
 xr=nm.positions(cfg,method);wr=nm.weights(len(xr));u0=np.sin(np.deg2rad(scan));rx=wr[None,:]*np.exp(-2j*np.pi*f[:,None]*xr[None,:]*u0/cfg.c)
 if method=='dense':C=coef[:,0,None]*rx;xp=xr
 else:
  xt=nm.positions(cfg,method,True);wt=nm.weights(len(xt));focus=np.sin(np.deg2rad(cfg.focus));A=np.zeros((len(f),len(xt)),complex)
  for k in range(6):A+=coef[:,k,None]*np.exp(-2j*np.pi*f[:,None]*xt[None,:]*focus[k]/cfg.c)*wt[None,:]
  C=np.zeros((len(f),4*(len(xr)-1)+len(xt)),complex)
  for n in range(len(xr)):C[:,4*n:4*n+len(xt)]+=A*rx[:,n,None]
  xp=(np.arange(C.shape[1])-(C.shape[1]-1)/2)*(cfg.c/cfg.fc/2)
 return f,C,xp

def characteristic(method,scan,rate=256):
 f,C,xp=coefficient_virtual_array(method,scan);order=np.argsort(f);f=f[order];C=C[order];dt=1/(cfg.fs*rate);t0=-.00008;m=int(round(.00016/dt))+1;t=t0+np.arange(m)*dt;df=f[1]-f[0]
 a=np.exp(-2j*np.pi*df*t0);w=np.exp(2j*np.pi*df*dt)
 char=signal.czt(C,m=m,w=w,a=a,axis=0)*np.exp(2j*np.pi*f[0]*t[:,None])
 return char,xp,t0,dt

def query(method,ch,angles):
 char,xp,t0,dt=ch;u=np.sin(np.deg2rad(angles));nodes,weights=element_nodes(method);width=cfg.element_width_lambda*cfg.c/cfg.fc;y=np.zeros(len(u),complex)
 # Columns are virtual physical Tx+Rx positions. Absolute-frequency phase was
 # included inside char(t), so there is NO additional fc phase here.
 for col,x in enumerate(xp):
  delay=(x+width*nodes[:,None])*u[None,:]/cfg.c;index=(delay-t0)/dt;lo=np.floor(index).astype(int);fraction=index-lo
  if lo.min()<0 or lo.max()+1>=len(char):raise ValueError('Outside frozen characteristic guard')
  value=char[lo,col]*(1-fraction)+char[lo+1,col]*fraction;y+=np.sum(weights[:,None]*value,axis=0)
 return y

def refined_angles(scan,values):
 mask=abs(coarse-scan)>3+1e-12;power=abs(values);segs=np.flatnonzero(mask);p=power.copy();p[~mask]=-np.inf;peaks,_=signal.find_peaks(p)
 candidates=np.unique(np.r_[coarse[peaks[mask[peaks]]],coarse[segs[0]],coarse[segs[-1]],coarse[segs[np.argmin(abs(coarse[segs]-scan))]],-50.,50.])
 # Both protection boundaries are included independently, not just one side.
 for sign in(-1,1):
  eligible=coarse[(sign*(coarse-scan)>3+1e-12)]
  if len(eligible):candidates=np.r_[candidates,eligible[np.argmin(abs(eligible-scan))]]
 fine=np.unique(np.round(np.concatenate([np.arange(c-.5,c+.50001,.01) for c in candidates]),8));fine=fine[(fine>=-50)&(fine<=50)&(abs(fine-scan)>3+1e-12)]
 local=np.round(np.arange(max(-50,scan-3),min(50,scan+3)+.00001,.01),8)
 return candidates,fine,local

validation=[]
for method in ws.METHODS:
 nodes,w=element_nodes(method);a=np.linspace(-.65,.65,1301);actual=np.sum(w[:,None]*np.exp(2j*np.pi*nodes[:,None]*a[None,:]),axis=0);expected=np.sinc(a)**(1 if method=='dense' else 2);err=float(np.max(abs(actual-expected)));validation.append({'method':method,'validation':'element_Gaussian_integral','scan_deg':None,'normalized_mag_L2':err,'off_peak_level_difference_db':None,'pass':err<=1e-12})
 for scan in cfg.focus:
  x=np.asarray(protocol['frozen_validation']['angles_for_native_check']);direct=ws.beam(method,scan,x,1)['all_simultaneous_true_delay'];ch=characteristic(method,scan,256);y=query(method,ch,x);err=float(np.linalg.norm(y-direct)/max(np.linalg.norm(direct),1e-30));validation.append({'method':method,'validation':'native8192_vs_partialCZT256','scan_deg':scan,'normalized_mag_L2':err,'off_peak_level_difference_db':None,'pass':err<=.001})
  coarse256=query(method,ch,coarse);cand,fine,local=refined_angles(scan,coarse256);off256=query(method,ch,fine);main256=query(method,ch,local)
  ch128=characteristic(method,scan,128);off128=query(method,ch128,fine);main128=query(method,ch128,local);d256=20*np.log10(max(abs(off256).max(),1e-30)/max(abs(main256).max(),1e-30));d128=20*np.log10(max(abs(off128).max(),1e-30)/max(abs(main128).max(),1e-30));diff=float(abs(d256-d128));validation.append({'method':method,'validation':'128_vs256_off_mainlobe_peak','scan_deg':scan,'normalized_mag_L2':None,'off_peak_level_difference_db':diff,'pass':diff<=.05})
 print('off-mainlobe validation '+method,flush=True)
ws.csvout(out/'numerical_validation.csv',validation)
if not all(r['pass'] for r in validation):
 ws.dump(out/'receipt.json',{'computed':False,'frozen_validation_failed':True,'validation':validation});raise SystemExit(1)
rows=[];profiles={}
for method in ws.METHODS:
 for ii,scan in enumerate(scans):
  ch=characteristic(method,scan,256);co=query(method,ch,coarse);candidates,fine,local=refined_angles(scan,co);off=query(method,ch,fine);main=query(method,ch,local);io=int(np.argmax(abs(off)));im=int(np.argmax(abs(main)));value=float(20*np.log10(max(abs(off[io]),1e-30)/max(abs(main[im]),1e-30)))
  rows.append({'method':method,'scan_angle_deg':scan,'coarse_candidates_including_boundaries':len(candidates),'refined_off_mainlobe_samples':len(fine),'protected_halfwidth_deg':3.,'maximum_off_mainlobe_response_db':value,'off_peak_target_angle_deg':float(fine[io]),'desired_local_peak_angle_deg':float(local[im]),'desired_local_peak_amplitude':float(abs(main[im])),'frequency_bins':8192,'characteristic_query_rate':256,'coarse_target_step_deg':.25,'refined_step_deg':.01})
  if scan in cfg.focus:profiles[f'{method}_scan{scan:g}_coarse_angles_deg']=coarse;profiles[f'{method}_scan{scan:g}_coarse_response']=co;profiles[f'{method}_scan{scan:g}_fine_off_angles_deg']=fine;profiles[f'{method}_scan{scan:g}_fine_off_response']=off
  if ii%20==0:print(f'off-mainlobe {method}:{ii+1}/{len(scans)}',flush=True)
 ws.csvout(out/'all_scan_off_mainlobe.csv',rows)
summary=[]
for method in ws.METHODS:
 group=[r for r in rows if r['method']==method];r=max(group,key=lambda r:r['maximum_off_mainlobe_response_db']);summary.append({**r,'audit_scan_directions':len(group),'definition':'coarse-to-fine maximum outside+/-3deg protected desired lobe; not first-null PSLR'})
ws.csvout(out/'worst_off_mainlobe_summary.csv',summary);np.savez_compressed(out/'sector_profiles_and_refined_candidates.npz',**profiles)
ws.dump(out/'receipt.json',{'computed':True,'rows':len(rows),'validation_rows':len(validation),'all_fixed_numerical_validation_pass':all(r['pass'] for r in validation),'summary':summary,'no_waveform_or_gate_retuning':True,'old_source_unchanged':hashlib.sha256((HERE/'run_wave_spatial.py').read_bytes()).hexdigest()==protocol['original_numerical_script_sha256']})
print(json.dumps(summary,ensure_ascii=False),flush=True)
