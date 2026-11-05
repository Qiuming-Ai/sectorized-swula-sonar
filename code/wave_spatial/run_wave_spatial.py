"""Final minimum waveform/complete two-way checks; old files are read-only.

The protocol and hashes are written before any outcome is computed. Full-band
quadrature preserves the zeroth and first discrete complex spectral moments in
every bin group and is checked against successively finer, frozen references.
"""
from __future__ import annotations
from pathlib import Path
from dataclasses import asdict
from datetime import datetime, timezone
import argparse, hashlib, json, sys, time, csv, shutil
import numpy as np
from scipy import signal
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

SOURCE=Path(__file__).resolve().parent
ROOT=SOURCE.parents[1]
OLD=ROOT/'code'
HERE=ROOT/'outputs/wave_spatial'
sys.path.insert(0,str(OLD/'common'))
import numerical_model as nm
import wave_point as wp
CFG=nm.Config()
METHODS=('dense','cdm','fdm')
COLORS={'dense':'#0072B2','cdm':'#D55E00','fdm':'#009E73'}
LABELS={'dense':'Dense LFM','cdm':'CDM-ZC','fdm':'FDM-LFM'}
STYLES={'dense':'-','cdm':'--','fdm':':'}
plt.rcParams.update({'font.family':'serif','font.size':10,'axes.labelsize':10,
 'legend.fontsize':9,'xtick.labelsize':9,'ytick.labelsize':9,'pdf.fonttype':42,
 'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':400})
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,data):Path(p).write_text(json.dumps(data,ensure_ascii=False,indent=2,default=lambda v:v.item() if isinstance(v,np.generic) else str(v))+'\n',encoding='utf-8')
def csvout(p,rows):
 rows=list(rows)
 if not rows:return
 with Path(p).open('w',encoding='utf-8-sig',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def figsave(fig,path):
 fig.tight_layout();fig.savefig(str(path)+'.pdf',bbox_inches='tight');fig.savefig(str(path)+'.png',bbox_inches='tight');plt.close(fig)
def dB(x,floor=-100):return np.maximum(20*np.log10(np.maximum(abs(x),1e-30)/max(float(np.max(abs(x))),1e-30)),floor)
def normalized(x):return abs(x)/max(float(np.max(abs(x))),1e-30)

def freeze(out):
 inputs=[Path(__file__).resolve(),OLD/'common/numerical_model.py',OLD/'common/wave_point.py']
 hashes={str(p.relative_to(ROOT)):sha(p) for p in inputs}
 protocol={'id':'RINENG_FINAL_MINIMUM_WAVE_SPATIAL_v01_20261010','frozen_utc':datetime.now(timezone.utc).isoformat(),'configuration':asdict(CFG),'source_hashes':hashes,
  'scope':'Only new wave_spatial directory written; common model/imported source and manuscript read-only.',
  'window_control':{'methods':['dense/Hamming','dense/rectangular','cdm/rectangular','fdm/Hamming','fdm/rectangular'],'all_common':'129-tap Kaiser8.6 Tx and Rx FIR, compensated64 sample Tx centre crop2000; actual simultaneous transmit-element total energy1; actual Tx+Rx unit-energy MF templates','FWHM':'full width at half maximum POWER of effective finite aperiodic template autocorrelation,32x correlation interpolation','sidelobes':'PSLR and ISLR exclude nearest local power minima around global main peak; not common equal-window optimization'},
  'complete_two_way':{'scan_directions_deg':sorted(set(np.arange(-50,51,1).tolist()+list(CFG.focus))),'target_angle_grid_deg':0.01,'target_domain_deg':[-50,50],'nfft':8192,'model':'At fixed receive scan angle theta0 and assigned matched-filter sector q, vary point-target angle theta. Sum ALL simultaneous Tx waveform contributions, true frequency-dependent TTD Tx and Rx factors, frequency-dependent Tx/Rx sinc element factors, and common Rx frontend. Dense uses the original omnidirectional single Tx reference, not invented12-element multibeam Tx.','isolated_diagnostic':'Only transmitted waveform k=q retained, same actual excitation normalization; no claim that this measures scene-level cross-waveform robustness.','approximation':'Freeze spatial Tx/Rx phase and element directivity of each transmitted component at fc (Dense/CDM) or its fk (FDM), retain exact complex spectral products and cross-waveform coherent sum.','quadrature':'Full8192FFT spectrum retained; contiguous sorted frequency groups represented at first/last frequencies with weights preserving exact group zeroth and first complex spectral moments. No PSD cutoff or omitted weak bins.','convergence_groups':[16,8,4,2,1],'fixed_convergence_criteria':{'maximum_normalized_magnitude_relative_L2':0.001,'maximum_global_peak_angle_difference_deg':0.010001,'maximum_alias_peak_level_difference_db':0.05},'convergence_scan_angles_deg':list(CFG.focus),'convergence_target_angle_step_deg':0.01,'reference':'At each stage compare candidate to next finer group; advance automatically if any check fails. Group1 is native frequency sum. Unfavourable checks retained.','grating_metric':'Visible order sin(theta_g)=sin(theta0)+p*c/(fref_q*dRx), p!=0; worst peak in +/-2deg neighbourhood clipped toFOV, relative to desired local peak within +/-3deg; a geometric neighbourhood diagnostic, not the maximum sidelobe everywhere.','angle_error':'Global output-beam peak relative to fixed receive scan bearing; local desired peak within3deg also reported; all quantities limited by0.01deg sampled grid.'},
  'ambiguity':{'finite_record_nfft':8192,'correlation_oversampling':16,'delay_scopes':['complete padded linear finite overlap record','abs(lag)<=0.6ms'],'velocities_mps':[-10,-5,-2,-1,0,1,2,5,10],'continuous_model':'beta=(c-v)/(c+v), positive receding; u(beta*t) on32x Fourier interpolant with512 guard samples, carrier shift(beta-1)fc once, unit received-pulse energy then common Rx frontend; FDM offset automatically time-scaled.','normalization':'For received channel i and matched-template j, divide by desired nominal channel-j zero-velocity MF peak; not a unit desired peak redefined after Doppler. All ordered distinct pairs retained.','Doppler_interpolation_checks':'16x vs32x received-waveform continuous-query comparison at velocities-10,0,10; fixed0.002dB peak/0.5mm peak-lag tolerances; failing results retained.'},
  'no_favorable_outcome_gate':True,'separate_spatial_and_waveform_mechanisms':True,'original_six_bands_and_codes_unchanged':True}
 dump(out/'protocol_before_run.json',protocol)
 snap=out/'code_snapshot';snap.mkdir()
 for p in inputs:
  if p.suffix=='.py':shutil.copy2(p,snap/p.name)
 return hashes

def window_control(out):
 rows=[];energy=[];arrays={};plot=[]
 for method,taper in [('dense','hamming'),('dense','rectangular'),('cdm','rectangular'),('fdm','hamming'),('fdm','rectangular')]:
  tx,meta=nm.waveforms(CFG,method,taper=taper)
  send,mf,nu,sm=nm.spectra(CFG,method,8192,taper=taper)
  excitation_energy,er,seq,active=wp._tx_energy_and_papr(CFG,method,send,nu)
  key=f'{method}_{taper}'
  energy.append({'method':method,'taper':taper,'actual_Tx_total_sample_energy':excitation_energy,'stored_Tx_samples':len(tx),'full_before_crop_samples':meta['full_convolution_samples'],'retained_pre_normalization_energy_fraction':meta['crop_retained_energy_fraction_before_normalization'],'Rx_matched_template':'unit-energy actual cropped Tx through same compensated continuous Rx FIR'})
  arrays[key+'_Tx']=tx;arrays[key+'_actual_Tx_element_sequences']=seq
  for k in range(tx.shape[1]):
   x,xmeta=wp._template_time(mf[:,k]);lag,ac=wp._autocorrelation(x,CFG,oversampling=32)
   rec={'method':method,'taper':taper,'channel':k+1,'power_FWHM_mm':nm.peak_width(lag*CFG.c/2*1000,abs(ac)**2),**wp._side_metrics(lag,ac),'effective_template_support_samples':len(x),'unit_template_sample_energy':float(np.sum(abs(x)**2)),'discarded_roundoff_energy_fraction':xmeta['discarded_roundoff_energy_fraction']}
   rows.append(rec);arrays[key+f'_ch{k+1}_template']=x;arrays[key+f'_ch{k+1}_correlation']=ac;arrays[key+f'_ch{k+1}_lag_s']=lag
   if k==0:plot.append((method,taper,lag,ac))
 csvout(out/'window_template_metrics.csv',rows);csvout(out/'window_energy_and_support.csv',energy);np.savez_compressed(out/'window_arrays.npz',**arrays)
 fig,axs=plt.subplots(1,2,figsize=(6.8,2.95))
 for method,taper,lag,ac in plot:
  ls=STYLES[method] if taper!='rectangular' or method=='cdm' else '-.'
  label=LABELS[method]+(' (no taper)' if method=='cdm' else f' ({taper})')
  for ax in axs:ax.plot(lag*CFG.c/2*1000,dB(ac),color=COLORS[method],ls=ls,lw=1.45,label=label)
 axs[0].set(xlim=(-90,90),ylim=(-45,1),xlabel='Equivalent range lag (mm)',ylabel='Correlation power (dB)')
 axs[1].set(xlim=(-900,900),ylim=(-65,1),xlabel='Equivalent range lag (mm)',ylabel='Correlation power (dB)');axs[1].legend(fontsize=8)
 for ax in axs:ax.grid(alpha=.18)
 figsave(fig,out/'window_fairness_control')
 return {'template_rows':len(rows),'energy_rows':len(energy),'maximum_Tx_energy_error':max(abs(r['actual_Tx_total_sample_energy']-1) for r in energy)}

def quadrature_nodes(f,coeff,group):
 order=np.argsort(f);ff=f[order];w=coeff[order]
 if group==1:return ff,w
 nodes=[];wts=[]
 for lo in range(0,len(ff),group):
  fs=ff[lo:lo+group];ws=w[lo:lo+group]
  if len(fs)==1:nodes.append(fs[0]);wts.append(ws[0]);continue
  frac=(fs-fs[0])/(fs[-1]-fs[0])
  nodes.extend([fs[0],fs[-1]]);wts.extend([np.sum(ws*(1-frac[:,None]),axis=0),np.sum(ws*frac[:,None],axis=0)])
 return np.asarray(nodes),np.asarray(wts)

def sector(theta):return int(np.clip(np.digitize(theta,CFG.edges[1:-1]),0,5))
def rx_pattern(method,f,target_u,scan_u):
 n=CFG.n_dense if method=='dense' else CFG.n_sparse;d=(.5 if method=='dense' else 2)*CFG.c/CFG.fc
 return nm.array_factor(2*np.pi*f[:,None]*d/CFG.c*(target_u[None,:]-scan_u),n)

def beam(method,scan_deg,angle,group):
 send,mf,nu,_=nm.spectra(CFG,method,8192);f=CFG.fc+nu
 q=0 if method=='dense' else sector(scan_deg)
 coeff=send*nm.frontend_response(CFG,len(f))[:,None]*np.conj(mf[:,q,None])/len(f)
 ff,ww=quadrature_nodes(f,coeff,group)
 u=np.sin(np.deg2rad(angle));u0=np.sin(np.deg2rad(scan_deg))
 exact=np.zeros(len(angle),complex);isolated=np.zeros_like(exact)
 for lo in range(0,len(ff),32):
  ft=ff[lo:lo+32];wt=ww[lo:lo+32]
  rx=rx_pattern(method,ft,u,u0)
  ep=np.sinc(CFG.element_width_lambda/CFG.fc*ft[:,None]*u[None,:])
  if method=='dense':
   tmp=wt[:,0,None]*rx*ep;exact+=tmp.sum(axis=0);isolated+=tmp.sum(axis=0)
  else:
   tx=nm.tx_pattern(CFG,ft,angle)
   exact+=np.sum(np.sum(wt[:,:,None]*tx,axis=1)*rx*ep,axis=0)
   isolated+=np.sum(wt[:,q,None]*tx[:,q,:]*rx*ep,axis=0)
 refs=np.array([CFG.fc]) if method=='dense' else (np.array(CFG.centres) if method=='fdm' else np.full(6,CFG.fc))
 approx=np.zeros_like(exact);approx_isolated=np.zeros_like(exact)
 for k,ref in enumerate(refs):
  rx=rx_pattern(method,np.array([ref]),u,u0)[0]
  ep=np.sinc(CFG.element_width_lambda/CFG.fc*ref*u)
  if method=='dense':tx=np.ones(len(angle))
  else:tx=nm.tx_pattern(CFG,np.array([ref]),angle)[0,k,:]
  value=coeff[:,k].sum()*rx*ep*tx;approx+=value
  if k==q:approx_isolated+=value
 return {'all_simultaneous_true_delay':exact,'all_simultaneous_centre_approx':approx,'isolated_assigned_true_delay':isolated,'isolated_assigned_centre_approx':approx_isolated}

def beam_stats(method,scan,angle,x):
 a=abs(x);i=int(np.argmax(a));inside=abs(angle-scan)<=3+1e-12;indices=np.flatnonzero(inside);il=int(indices[np.argmax(a[indices])]);ref=float(a[il])
 n=CFG.n_dense if method=='dense' else CFG.n_sparse;d=(.5 if method=='dense' else 2)*CFG.c/CFG.fc
 fq=CFG.centres[sector(scan)] if method=='fdm' else CFG.fc
 aliases=[]
 for p in range(-5,6):
  if p==0:continue
  ug=np.sin(np.deg2rad(scan))+p*CFG.c/(fq*d)
  if abs(ug)>1:continue
  ag=float(np.rad2deg(np.arcsin(ug)))
  if not CFG.fov[0]<=ag<=CFG.fov[1]:continue
  mask=abs(angle-ag)<=2+1e-12
  if mask.any():aliases.append((float(a[mask].max()),ag,p))
 worst=max(aliases,default=(0.,np.nan,0))
 dbgl=float(20*np.log10(max(worst[0],1e-30)/max(ref,1e-30))) if aliases else None
 return {'global_peak_angle_deg':float(angle[i]),'global_peak_error_deg':float(angle[i]-scan),'desired_local_peak_angle_deg':float(angle[il]),'desired_local_peak_error_deg':float(angle[il]-scan),'local_power_HPBW_deg':nm.peak_width(angle[inside],a[inside]**2),'worst_visible_alias_neighborhood_db':dbgl,'worst_alias_nominal_angle_deg':worst[1] if aliases else None,'worst_alias_order':worst[2] if aliases else None,'visible_alias_count':len(aliases),'local_peak_amplitude':ref,'global_to_local_peak_db':float(20*np.log10(max(a[i],1e-30)/max(ref,1e-30))),'power_HPBW_edge_censored':nm.peak_width(angle[inside],a[inside]**2) is None}

def compare_profiles(a,b):
 x=normalized(a);y=normalized(b);den=np.linalg.norm(y)
 return {'normalized_magnitude_relative_L2':float(np.linalg.norm(x-y)/max(den,1e-30)),'magnitude_profile_correlation':float(np.corrcoef(x,y)[0,1])}

def two_way(out):
 angle=np.round(np.arange(-50,50.00001,.01),8);conv=[];chosen=None
 for group,finer in zip((16,8,4,2),(8,4,2,1)):
  checks=[]
  for method in METHODS:
   for scan in CFG.focus:
    a=beam(method,scan,angle,group);b=beam(method,scan,angle,finer)
    for scope in ('all_simultaneous_true_delay','isolated_assigned_true_delay'):
     sa=beam_stats(method,scan,angle,a[scope]);sb=beam_stats(method,scan,angle,b[scope]);diff=compare_profiles(a[scope],b[scope])
     glerr=abs(sa['worst_visible_alias_neighborhood_db']-sb['worst_visible_alias_neighborhood_db']) if sa['worst_visible_alias_neighborhood_db'] is not None else 0.
     peakerr=abs(sa['global_peak_angle_deg']-sb['global_peak_angle_deg'])
     passed=diff['normalized_magnitude_relative_L2']<=.001 and peakerr<=.010001 and glerr<=.05
     row={'method':method,'scan_angle_deg':scan,'scope':scope,'group_bins':group,'reference_group_bins':finer,**diff,'peak_angle_difference_deg':peakerr,'alias_level_difference_db':glerr,'pass_fixed_criteria':passed};conv.append(row);checks.append(passed)
  print(f'quadrature {group} vs {finer}: {sum(checks)}/{len(checks)} pass',flush=True)
  csvout(out/'two_way_quadrature_convergence.csv',conv)
  if all(checks):chosen=group;break
 if chosen is None:chosen=1
 rows=[];comparisons=[];arrays={'target_angles_deg':angle}; scans=sorted(set(np.arange(-50,51,1).tolist()+list(CFG.focus)))
 for method in METHODS:
  for ii,scan in enumerate(scans):
   profiles=beam(method,scan,angle,chosen)
   for scope,x in profiles.items():rows.append({'method':method,'receive_scan_angle_deg':scan,'scope':scope,'target_angle_grid_deg':.01,**beam_stats(method,scan,angle,x)})
   for scope in ('all_simultaneous','isolated_assigned'):
    ex=profiles[scope+'_true_delay'];ap=profiles[scope+'_centre_approx'];se=beam_stats(method,scan,angle,ex);sa=beam_stats(method,scan,angle,ap)
    comparisons.append({'method':method,'receive_scan_angle_deg':scan,'scope':scope,**compare_profiles(ex,ap),'global_peak_angle_difference_deg':float(se['global_peak_angle_deg']-sa['global_peak_angle_deg']),'desired_local_peak_angle_difference_deg':float(se['desired_local_peak_angle_deg']-sa['desired_local_peak_angle_deg']),'full_minus_approx_alias_level_db':(se['worst_visible_alias_neighborhood_db']-sa['worst_visible_alias_neighborhood_db']) if se['worst_visible_alias_neighborhood_db'] is not None else None})
   if scan in CFG.focus or scan in(-50.,0.,50.):
    for name,x in profiles.items():arrays[f'{method}_scan{scan:g}_{name}']=x
   if ii%10==0:print(f'two-way {method} scan {ii+1}/{len(scans)}',flush=True)
  csvout(out/'two_way_all_beam_metrics.csv',rows);csvout(out/'two_way_approximation_comparison.csv',comparisons)
 np.savez_compressed(out/'two_way_representative_profiles.npz',**arrays)
 summary=[]
 for method in METHODS:
  for scope in ('all_simultaneous','isolated_assigned'):
   c=[r for r in comparisons if r['method']==method and r['scope']==scope];tr=[r for r in rows if r['method']==method and r['scope']==scope+'_true_delay'];gl=[r for r in tr if r['worst_visible_alias_neighborhood_db'] is not None]
   summary.append({'method':method,'scope':scope,'scan_directions':len(c),'minimum_magnitude_profile_correlation':min(r['magnitude_profile_correlation'] for r in c),'maximum_normalized_magnitude_relative_L2':max(r['normalized_magnitude_relative_L2'] for r in c),'maximum_absolute_global_peak_angle_difference_deg':max(abs(r['global_peak_angle_difference_deg']) for r in c),'maximum_absolute_desired_local_peak_angle_difference_deg':max(abs(r['desired_local_peak_angle_difference_deg']) for r in c),'maximum_absolute_true_delay_global_peak_error_deg':max(abs(r['global_peak_error_deg']) for r in tr),'worst_true_delay_visible_alias_neighborhood_db':max((r['worst_visible_alias_neighborhood_db'] for r in gl),default=None),'maximum_abs_alias_level_change_db':max((abs(r['full_minus_approx_alias_level_db']) for r in c if r['full_minus_approx_alias_level_db'] is not None),default=None),'target_angle_grid_deg':.01,'selected_frequency_quadrature_group_bins':chosen})
 csvout(out/'two_way_summary.csv',summary)
 fig,axs=plt.subplots(1,3,figsize=(7.1,2.7))
 scan=CFG.focus[-1]
 for ax,method in zip(axs[:2],('cdm','fdm')):
  for model,ls in [('true_delay','-'),('centre_approx','--')]:
   x=arrays[f'{method}_scan{scan:g}_all_simultaneous_{model}'];ax.plot(angle,dB(x),color=COLORS[method],ls=ls,lw=1.4,label='True delay' if model=='true_delay' else 'Centre approximation')
  ax.set(xlim=(-50,50),ylim=(-60,2),xlabel='Target bearing (deg)',ylabel='Beam power (dB)',title=LABELS[method]+f', scan {scan:g}°');ax.legend(fontsize=8);ax.grid(alpha=.15)
 for method in ('cdm','fdm'):
  c=[r for r in comparisons if r['method']==method and r['scope']=='all_simultaneous' and r['receive_scan_angle_deg'] in CFG.focus]
  axs[2].plot([r['receive_scan_angle_deg'] for r in c],[r['magnitude_profile_correlation'] for r in c],color=COLORS[method],ls=STYLES[method],marker='s' if method=='cdm' else '^',label=LABELS[method],lw=1.4)
 axs[2].set(xlabel='Receive scan bearing (deg)',ylabel='Magnitude-profile correlation',ylim=(.75,1.002));axs[2].legend(fontsize=8);axs[2].grid(alpha=.15)
 figsave(fig,out/'complete_two_way_spatial_check')
 dump(out/'two_way_scope.json',{'response_formula':'H_q(theta;theta0)=sum_f sum_k S_k(f) H_Rx(f) conj(M_q(f)) B_T,k(theta,f) B_R(theta,theta0,f) E_R(theta,f)/Nfft. B_T includesTx element directivity. Matched filter sector q determined bytheta0.','Dense':'Original omnidirectional single Tx +128 receive0.5lambda; not a different new dense Tx array.','spectrum':'All8192 complex frequency products retained via validated moment-preserving quadrature.','spatial_vs_cross':'All simultaneous waveform sum is primary; isolated assigned waveform is explicitly a diagnostic.','field_of_view':'Only +/-50deg evaluated. Endpoint HPBW may be censored, never completed with unreported outside-FOV samples.','grating':'Alias neighborhood metric relative to desired local peak; not equal to largest sidelobe anywhere.','phase':'Physical absolute-frequency TTD at f=fc+nu includedonce.','grid_resolution_deg':.01,'quadrature_group_bins':chosen})
 return {'selected_quadrature_group_bins':chosen,'convergence_rows':len(conv),'beam_rows':len(rows),'comparison_rows':len(comparisons),'summary':summary}

def ambiguity(out):
 rows=[];summary=[];arrays={};conv=[];velocities=(-10.,-5.,-2.,-1.,0.,1.,2.,5.,10.);nfft=8192
 for method in METHODS:
  tx,_=nm.waveforms(CFG,method);send,mf,nu,_=nm.spectra(CFG,method,nfft);g=nm.frontend_response(CFG,nfft);nch=tx.shape[1]
  references=np.abs(np.sum(np.fft.fft(tx,nfft,axis=0)*g[:,None]*np.conj(mf),axis=0)/nfft)
  for i in range(nch):
   padded=np.pad(tx[:,i],(512,nfft-512-len(tx)));high32=signal.resample(padded,nfft*32);high16=signal.resample(padded,nfft*16)
   for v in velocities:
    beta=(CFG.c-v)/(CFG.c+v);fd=(beta-1)*CFG.fc
    echo,_=wp._continuous_tx_echo(tx[:,i],beta,fd,32,CFG,nfft,high32);yf=np.fft.fft(echo)*g
    echo16=None
    if v in(-10.,0.,10.):echo16,_=wp._continuous_tx_echo(tx[:,i],beta,fd,16,CFG,nfft,high16);yf16=np.fft.fft(echo16)*g
    for j in range(nch):
     lag,corr,ig,il=wp._global_and_local_correlation(yf,mf[:,j],CFG,os_corr=16);a=abs(corr)/references[j]
     row={'method':method,'Tx_channel_i':i+1,'MF_channel_j':j+1,'pair_type':'auto' if i==j else 'distinct_cross','velocity_mps_positive_receding':v,'beta':beta,'carrier_shift_hz':fd,'normalized_global_peak_db':float(20*np.log10(max(a[ig],1e-30))),'global_peak_lag_ms':float(lag[ig]*1000),'normalized_local_0p6ms_peak_db':float(20*np.log10(max(a[il],1e-30))),'local_0p6ms_peak_lag_ms':float(lag[il]*1000),'normalization_reference_channel_j_nominal_MF_peak':float(references[j])};rows.append(row)
     if v==0 or (i!=j and v in(-10.,-2.,2.,10.)):arrays[f'{method}_i{i+1}_j{j+1}_v{v:g}_correlation']=corr/references[j]
     if echo16 is not None:
      lag16,c16,ig16,il16=wp._global_and_local_correlation(yf16,mf[:,j],CFG,os_corr=16)
      peakdiff=abs(20*np.log10(max(abs(c16[ig16]),1e-30)/max(abs(corr[ig]),1e-30)));biasdiff=abs(lag16[ig16]-lag[ig])*CFG.c/2*1000
      conv.append({'method':method,'Tx_channel_i':i+1,'MF_channel_j':j+1,'velocity_mps':v,'reference_query_oversampling':32,'tested_query_oversampling':16,'global_peak_difference_db':float(peakdiff),'peak_equivalent_range_difference_mm':float(biasdiff),'pass_fixed_0p002db_0p5mm':peakdiff<=.002 and biasdiff<=.5})
  for v in velocities:
   for pairtype in('auto','distinct_cross'):
    selected=[r for r in rows if r['method']==method and r['velocity_mps_positive_receding']==v and r['pair_type']==pairtype]
    if not selected:continue
    wg=max(selected,key=lambda r:r['normalized_global_peak_db']);wl=max(selected,key=lambda r:r['normalized_local_0p6ms_peak_db'])
    summary.append({'method':method,'pair_type':pairtype,'velocity_mps_positive_receding':v,'worst_global_peak_db':wg['normalized_global_peak_db'],'global_Tx_i':wg['Tx_channel_i'],'global_MF_j':wg['MF_channel_j'],'global_peak_lag_ms':wg['global_peak_lag_ms'],'worst_local_0p6ms_peak_db':wl['normalized_local_0p6ms_peak_db'],'local_Tx_i':wl['Tx_channel_i'],'local_MF_j':wl['MF_channel_j'],'local_peak_lag_ms':wl['local_0p6ms_peak_lag_ms']})
  print(f'ambiguity {method}: {nch*nch*len(velocities)} ordered pair/velocity rows',flush=True)
 arrays['lag_s']=lag;csvout(out/'finite_ambiguity_all_ordered_pairs.csv',rows);csvout(out/'finite_ambiguity_worst_pairs.csv',summary);csvout(out/'continuous_Doppler_pair_query_convergence.csv',conv);np.savez_compressed(out/'finite_ambiguity_correlations.npz',**arrays)
 fig,axs=plt.subplots(2,2,figsize=(6.8,5.0))
 rz=next(r for r in summary if r['method']=='fdm' and r['pair_type']=='distinct_cross' and r['velocity_mps_positive_receding']==0)
 key=f'fdm_i{rz["global_Tx_i"]}_j{rz["global_MF_j"]}_v0_correlation';c=arrays[key];yd=20*np.log10(np.maximum(abs(c),1e-30))
 axs[0,0].plot(lag*1000,yd,color=COLORS['fdm'],lw=1.5);axs[0,0].axvspan(-.6,.6,color='.9');axs[0,0].set(xlim=(-5.35,5.35),ylim=(-110,2),xlabel='Delay (ms)',ylabel='Normalized cross power (dB)',title='(a) Complete delay support')
 # The local worst pair need not equal the global worst pair.
 kl=f'fdm_i{rz["local_Tx_i"]}_j{rz["local_MF_j"]}_v0_correlation';cl=arrays[kl]
 axs[0,1].plot(lag*1000,20*np.log10(np.maximum(abs(cl),1e-30)),color=COLORS['fdm'],lw=1.5);axs[0,1].set(xlim=(-.6,.6),ylim=(-100,-20),xlabel='Delay (ms)',ylabel='Normalized cross power (dB)',title='(b) Local worst FDM pair')
 for method in('cdm','fdm'):
  rr=[r for r in summary if r['method']==method and r['pair_type']=='distinct_cross']
  for mode,style in(('global','-'),('local_0p6ms','--')):
   axs[1,0].plot([r['velocity_mps_positive_receding'] for r in rr],[r['worst_'+mode+'_peak_db'] for r in rr],color=COLORS[method],ls=style,marker='s' if method=='cdm' else '^',markersize=3,lw=1.3,label=LABELS[method]+(' full' if mode=='global' else ' local'))
 axs[1,0].set(xlabel='Radial velocity (m/s; receding +)',ylabel='Worst distinct-pair peak (dB)',ylim=(-100,0),title='(c) Continuous time-scaled Doppler');axs[1,0].legend(fontsize=8)
 tx,_=nm.waveforms(CFG,'fdm');f=CFG.fc+np.fft.fftshift(np.fft.fftfreq(65536,1/CFG.fs));styles=['-','--',':','-.','-','--'];palette=['#0072B2','#D55E00','#009E73','#CC79A7','#E69F00','#56B4E9']
 for k in range(6):
  p=abs(np.fft.fftshift(np.fft.fft(tx[:,k],65536)))**2;axs[1,1].plot(f/1000,10*np.log10(np.maximum(p/p.max(),1e-12)),color=palette[k],ls=styles[k],lw=1.3,label=f'{CFG.centres[k]/1000:g} kHz')
 axs[1,1].set(xlim=(395,505),ylim=(-70,2),xlabel='RF frequency (kHz)',ylabel='Normalized Tx power (dB)',title='(d) Actual FIR-shaped FDM pulses');axs[1,1].legend(fontsize=8,ncol=2)
 for ax in axs.ravel():ax.grid(alpha=.15)
 figsave(fig,out/'actual_finite_pulse_ambiguity_scope')
 return {'ordered_pair_velocity_rows':len(rows),'worst_pair_rows':len(summary),'query_convergence_rows':len(conv),'query_convergence_pass':sum(r['pass_fixed_0p002db_0p5mm'] for r in conv),'zero_velocity_distinct_pair_summary':[r for r in summary if r['pair_type']=='distinct_cross' and r['velocity_mps_positive_receding']==0]}

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--only',nargs='+',choices=['window','two_way','ambiguity']);args=parser.parse_args();selected=args.only or ['window','ambiguity','two_way']
 stamp=datetime.now().strftime('%Y%m%d_%H%M%S');out=HERE/('run_'+stamp);out.mkdir(parents=True);hashes=freeze(out);started=time.perf_counter();results={}
 for name in selected:
  print('START '+name,flush=True);sub=out/name;sub.mkdir();t=time.perf_counter()
  try:data={'window':window_control,'two_way':two_way,'ambiguity':ambiguity}[name](sub);data['elapsed_s']=time.perf_counter()-t;data['computed']=True;results[name]=data;dump(sub/'receipt.json',data);print('FINISH '+name+' '+json.dumps(data,ensure_ascii=False),flush=True)
  except Exception as exc:
   import traceback;traceback.print_exc();results[name]={'computed':False,'error':repr(exc)};dump(sub/'receipt.json',results[name])
 unchanged=all(sha(ROOT/p)==h for p,h in hashes.items())
 receipt={'run':str(out),'results':results,'elapsed_s':time.perf_counter()-started,'all_input_sources_unchanged':unchanged,'hashes':hashes,'author_validation':'pending','all_requested_computed':all(v.get('computed') for v in results.values())}
 dump(out/'receipt.json',receipt);artifacts=[{'path':str(p.relative_to(out)),'bytes':p.stat().st_size,'sha256':sha(p)} for p in out.rglob('*') if p.is_file()];dump(out/'artifact_hashes.json',artifacts);print(json.dumps(receipt,ensure_ascii=False),flush=True)
 return 0 if unchanged and receipt['all_requested_computed'] else 1
if __name__=='__main__':raise SystemExit(main())
