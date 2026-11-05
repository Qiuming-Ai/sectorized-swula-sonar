"""Frozen numerical refinement only; first16x/32x weak-pair failure retained."""
from pathlib import Path
from datetime import datetime, timezone
import sys, hashlib, json
import numpy as np
from scipy import signal
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE))
import run_wave_spatial as ws
nm=ws.nm;wp=ws.wp;cfg=ws.CFG
out=ws.HERE/'continuous_query_refinement_64';out.mkdir(parents=True,exist_ok=False)
protocol={'id':'CONTINUOUS_QUERY_REFERENCE_REFINEMENT_ONLY_v01','frozen_utc':datetime.now(timezone.utc).isoformat(),'reason':'First16x/32x weak FDM distinct pair failed original fixed peak-amplitude and argmax-position gate. Do not change waveform or tolerance, retain that failure.','transmitted_waveforms':'exact original canonical same Tx primary5ms pulses; no parameter retuning','velocities':[-10.,-5.,-2.,-1.,0.,1.,2.,5.,10.],'comparisons':['32x versus64x all ordered pairs and velocities','only failed32x/64x rows, compare64x versus128x under same criteria'],'fixed_amplitude_tolerance_db':.002,'fixed_peak_equivalent_range_tolerance_mm':.5,'correlation_oversampling':16,'global_near_equal_peak_diagnostic':'Top two local maxima of canonical magnitude correlation; peak-strength gap reported. No truth-conditioned localization.','no_favorable_gate':True,'same_model_sha256':nm.sha(nm.__file__),'input_script_sha256':hashlib.sha256((HERE/'run_wave_spatial.py').read_bytes()).hexdigest()}
ws.dump(out/'protocol_before_run.json',protocol)
rows=[];summaries=[];failed=[];nfft=8192;g=nm.frontend_response(cfg,nfft)
for method in ws.METHODS:
 tx,_=nm.waveforms(cfg,method);_,mf,_,_=nm.spectra(cfg,method,nfft);nch=tx.shape[1]
 for i in range(nch):
  padded=np.pad(tx[:,i],(512,nfft-512-len(tx)));high={rate:signal.resample(padded,nfft*rate) for rate in(32,64)}
  for v in protocol['velocities']:
   beta=(cfg.c-v)/(cfg.c+v);fd=(beta-1)*cfg.fc;ys={}
   for rate in(32,64):
    echo,_=wp._continuous_tx_echo(tx[:,i],beta,fd,rate,cfg,nfft,high[rate]);ys[rate]=np.fft.fft(echo)*g
   for j in range(nch):
    lag,a,ia,ila=wp._global_and_local_correlation(ys[32],mf[:,j],cfg,os_corr=16);lag,b,ib,ilb=wp._global_and_local_correlation(ys[64],mf[:,j],cfg,os_corr=16)
    amp=abs(20*np.log10(max(abs(a[ia]),1e-30)/max(abs(b[ib]),1e-30)));bias=abs(lag[ia]-lag[ib])*cfg.c/2*1000
    maxima,_=signal.find_peaks(abs(b));order=maxima[np.argsort(abs(b[maxima]))[::-1]]
    gap=float(20*np.log10(abs(b[order[0]])/max(abs(b[order[1]]),1e-30))) if len(order)>1 else None
    row={'method':method,'Tx_i':i+1,'MF_j':j+1,'velocity_mps':v,'candidate_query_rate':32,'reference_query_rate':64,'peak_strength_difference_db':float(amp),'global_peak_equivalent_range_difference_mm':float(bias),'canonical_peak_lag_ms':float(lag[ib]*1000),'top_two_local_peak_gap_db':gap,'pass_original_criteria':bool(amp<=.002 and bias<=.5)};rows.append(row)
    if not row['pass_original_criteria']:failed.append((method,i,j,v))
  print(f'{method} Tx{i+1}:32x/64x refinement',flush=True)
for method,i,j,v in failed:
 tx,_=nm.waveforms(cfg,method);_,mf,_,_=nm.spectra(cfg,method,nfft);padded=np.pad(tx[:,i],(512,nfft-512-len(tx)));beta=(cfg.c-v)/(cfg.c+v);fd=(beta-1)*cfg.fc;cs={}
 for rate in(64,128):
  high=signal.resample(padded,nfft*rate);echo,_=wp._continuous_tx_echo(tx[:,i],beta,fd,rate,cfg,nfft,high);lag,c,ig,il=wp._global_and_local_correlation(np.fft.fft(echo)*g,mf[:,j],cfg,os_corr=16);cs[rate]=(c,ig)
 a,ia=cs[64];b,ib=cs[128];amp=abs(20*np.log10(max(abs(a[ia]),1e-30)/max(abs(b[ib]),1e-30)));bias=abs(lag[ia]-lag[ib])*cfg.c/2*1000
 maxima,_=signal.find_peaks(abs(b));order=maxima[np.argsort(abs(b[maxima]))[::-1]];gap=float(20*np.log10(abs(b[order[0]])/max(abs(b[order[1]]),1e-30))) if len(order)>1 else None
 rows.append({'method':method,'Tx_i':i+1,'MF_j':j+1,'velocity_mps':v,'candidate_query_rate':64,'reference_query_rate':128,'peak_strength_difference_db':float(amp),'global_peak_equivalent_range_difference_mm':float(bias),'canonical_peak_lag_ms':float(lag[ib]*1000),'top_two_local_peak_gap_db':gap,'pass_original_criteria':bool(amp<=.002 and bias<=.5)})
ws.csvout(out/'all_ordered_pair_reference_refinement.csv',rows)
ws.dump(out/'receipt.json',{'rows':len(rows),'stage32vs64_rows':sum(r['candidate_query_rate']==32 for r in rows),'stage32vs64_fails':sum(not r['pass_original_criteria'] and r['candidate_query_rate']==32 for r in rows),'stage64vs128_rows':sum(r['candidate_query_rate']==64 for r in rows),'stage64vs128_fails':sum(not r['pass_original_criteria'] and r['candidate_query_rate']==64 for r in rows),'maximum32vs64_amplitude_difference_db':max(r['peak_strength_difference_db'] for r in rows if r['candidate_query_rate']==32),'old_failure_preserved':True,'source_unchanged':nm.sha(nm.__file__)==protocol['same_model_sha256']})
print(json.dumps(json.loads((out/'receipt.json').read_text(encoding='utf-8')),ensure_ascii=False),flush=True)
