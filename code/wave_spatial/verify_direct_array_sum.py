"""Independent exact-bin explicit-array check of the full two-way formula."""
from pathlib import Path
import sys, hashlib, json
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
import run_wave_spatial as ws
nm=ws.nm;cfg=ws.CFG
out=ws.HERE/'direct_array_sum_validation';out.mkdir(parents=True,exist_ok=False)
protocol={'formula':'Explicit weighted sums over physical Tx and Rx elements and all8192FFT bins, no array-factor closed form, no frequency quadrature. Tx and Rx frequency-dependent rectangular-element responses retained.','scan_deg':list(cfg.focus),'test_target_deg':[-50.,-41.25,-33.333333,-24.25,-8.,0.,8.,24.25,33.333333,41.25,50.],'criteria_relative_L2':1e-11,'input_code_sha256':hashlib.sha256((HERE/'run_wave_spatial.py').read_bytes()).hexdigest(),'immutable_model_sha256':nm.sha(nm.__file__)}
ws.dump(out/'protocol_before_validation.json',protocol)
theta=np.asarray(protocol['test_target_deg']);u=np.sin(np.deg2rad(theta));rows=[]
for method in ws.METHODS:
 send,mf,nu,_=nm.spectra(cfg,method,8192);f=cfg.fc+nu;xr=nm.positions(cfg,method);wr=nm.weights(len(xr));xt=nm.positions(cfg,method,True);wt=nm.weights(len(xt))
 for scan in cfg.focus:
  q=0 if method=='dense' else ws.sector(scan);u0=np.sin(np.deg2rad(scan));coeff=send*nm.frontend_response(cfg,len(f))[:,None]*np.conj(mf[:,q,None])/len(f)
  independent=np.zeros(len(theta),complex)
  for j,uj in enumerate(u):
   receive=np.exp(2j*np.pi*f[:,None]*xr[None,:]/cfg.c*(uj-u0))@wr
   ep=np.sinc(cfg.element_width_lambda/cfg.fc*f*uj)
   if method=='dense':independent[j]=np.sum(coeff[:,0]*receive*ep)
   else:
    for k,focus in enumerate(cfg.focus):
     tx=np.exp(2j*np.pi*f[:,None]*xt[None,:]/cfg.c*(uj-np.sin(np.deg2rad(focus))))@wt
     independent[j]+=np.sum(coeff[:,k]*tx*receive*ep**2)
  analytic=ws.beam(method,scan,theta,1)['all_simultaneous_true_delay']
  err=float(np.linalg.norm(analytic-independent)/max(np.linalg.norm(independent),1e-30));rows.append({'method':method,'scan_angle_deg':scan,'relative_complex_L2_error':err,'max_absolute_error':float(np.max(abs(analytic-independent))),'pass_fixed_1e_minus11':err<=1e-11})
ws.csvout(out/'explicit_element_sum_vs_closed_array_factor.csv',rows)
ws.dump(out/'receipt.json',{'rows':len(rows),'all_pass':all(r['pass_fixed_1e_minus11'] for r in rows),'maximum_relative_complex_L2_error':max(r['relative_complex_L2_error'] for r in rows),'source_unchanged':hashlib.sha256((HERE/'run_wave_spatial.py').read_bytes()).hexdigest()==protocol['input_code_sha256']})
print(json.dumps({'rows':len(rows),'all_pass':all(r['pass_fixed_1e_minus11'] for r in rows),'maximum_relative_complex_L2_error':max(r['relative_complex_L2_error'] for r in rows)},ensure_ascii=False),flush=True)
