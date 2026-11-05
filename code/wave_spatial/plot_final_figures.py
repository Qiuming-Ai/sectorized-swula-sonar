"""Publication-size figure replot; consumes raw arrays, changes no metric."""
from pathlib import Path
import argparse,sys,csv,hashlib,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE));import run_wave_spatial as ws
P=argparse.ArgumentParser()
P.add_argument('--run',type=Path,default=ws.ROOT/'data/waveforms')
P.add_argument('--out',type=Path,default=ws.ROOT/'outputs/figures/figure09')
args=P.parse_args();run=args.run;out=args.out;out.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.family':'serif','font.size':10,'axes.labelsize':10,'axes.titlesize':10,'xtick.labelsize':9,'ytick.labelsize':9,'legend.fontsize':9,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':400})
manifest={'purpose':'Plot-only re-export; no numeric rerun or changed selection. Designed at6in for>=8pt after>=0.9 proportional insertion; use at least5.4in publication width.','actual_width_in':6.0,'minimum_text_font_pt':9,'suggested_minimum_inserted_width_in':5.4,'inputs':{},'series_encoding':'color-blind safe Okabe-Ito plus line style/marker'}
def rows(p):return list(csv.DictReader(Path(p).open(encoding='utf-8-sig')))
def save(fig,name):
 fig.savefig(out/(name+'.pdf'),bbox_inches='tight');fig.savefig(out/(name+'.png'),bbox_inches='tight');plt.close(fig)
def register(p):manifest['inputs'][str(p.relative_to(run))]=hashlib.sha256(p.read_bytes()).hexdigest()

window=run/'window/window_arrays.npz'
if window.exists():
 register(window);data=np.load(window);fig,axs=plt.subplots(1,2,figsize=(6,3.6));legend=[]
 for method,taper in [('dense','hamming'),('dense','rectangular'),('cdm','rectangular'),('fdm','hamming'),('fdm','rectangular')]:
  key=f'{method}_{taper}_ch1';lag=data[key+'_lag_s'];ac=data[key+'_correlation'];ls=ws.STYLES[method] if taper!='rectangular' or method=='cdm' else '-.';label=ws.LABELS[method]+(' (no taper)' if method=='cdm' else f' ({taper})')
  for ax in axs:line,=ax.plot(lag*ws.CFG.c/2*1000,ws.dB(ac),color=ws.COLORS[method],ls=ls,lw=1.35,label=label)
  legend.append(line)
 axs[0].set(xlim=(-90,90),ylim=(-45,1),xlabel='Equivalent range lag (mm)',ylabel='Correlation power (dB)',title='(a) Mainlobe comparison')
 axs[1].set(xlim=(-900,900),ylim=(-65,1),xlabel='Equivalent range lag (mm)',ylabel='Correlation power (dB)',title='(b) Finite-pulse sidelobes')
 for ax in axs:ax.grid(alpha=.15)
 fig.legend(handles=legend,loc='lower center',ncol=2,bbox_to_anchor=(.5,.005),frameon=False);fig.subplots_adjust(bottom=.42,wspace=.34,top=.88);save(fig,'window_fairness_control')

ambiguity=run/'ambiguity/finite_ambiguity_correlations.npz'
if ambiguity.exists():
 register(ambiguity);data=np.load(ambiguity);lag=data['lag_s'];summary=rows(run/'ambiguity/finite_ambiguity_worst_pairs.csv');rz=next(r for r in summary if r['method']=='fdm' and r['pair_type']=='distinct_cross' and float(r['velocity_mps_positive_receding'])==0)
 fig,axs=plt.subplots(2,2,figsize=(6,4.9))
 key=f'fdm_i{rz["global_Tx_i"]}_j{rz["global_MF_j"]}_v0_correlation';c=data[key];axs[0,0].plot(lag*1000,20*np.log10(np.maximum(abs(c),1e-30)),color=ws.COLORS['fdm'],lw=1.25);axs[0,0].axvspan(-.6,.6,color='.9');axs[0,0].set(xlim=(-5.35,5.35),ylim=(-110,2),xlabel='Delay (ms)',ylabel='Normalized cross power (dB)',title='(a) Full effective-pulse support')
 key=f'fdm_i{rz["local_Tx_i"]}_j{rz["local_MF_j"]}_v0_correlation';c=data[key];axs[0,1].plot(lag*1000,20*np.log10(np.maximum(abs(c),1e-30)),color=ws.COLORS['fdm'],lw=1.25);axs[0,1].set(xlim=(-.6,.6),ylim=(-105,-65),xlabel='Delay (ms)',ylabel='Normalized cross power (dB)',title='(b) Evaluated local delay interval')
 for method in('cdm','fdm'):
  rr=[r for r in summary if r['method']==method and r['pair_type']=='distinct_cross']
  for mode,ls in(('global','-'),('local_0p6ms','--')):
   axs[1,0].plot([float(r['velocity_mps_positive_receding']) for r in rr],[float(r['worst_'+mode+'_peak_db']) for r in rr],color=ws.COLORS[method],ls=ls,marker='s' if method=='cdm' else '^',markersize=3,lw=1.2,label=method.upper()+(' full' if mode=='global' else ' local'))
 axs[1,0].set(xlabel='Radial velocity (m/s; receding +)',ylabel='Worst distinct-pair peak (dB)',ylim=(-100,0),title='(c) Continuous Doppler model');axs[1,0].legend(fontsize=9,ncol=2,loc='center',bbox_to_anchor=(.5,.45),framealpha=.9,handlelength=1.4,columnspacing=.8)
 tx,_=ws.nm.waveforms(ws.CFG,'fdm');f=ws.CFG.fc+np.fft.fftshift(np.fft.fftfreq(65536,1/ws.CFG.fs));styles=['-','--',':','-.','-','--'];palette=['#0072B2','#D55E00','#009E73','#CC79A7','#E69F00','#56B4E9']
 for k in range(6):
  p=abs(np.fft.fftshift(np.fft.fft(tx[:,k],65536)))**2;axs[1,1].plot(f/1000,10*np.log10(np.maximum(p/p.max(),1e-12)),color=palette[k],ls=styles[k],lw=1.15,label=f'{ws.CFG.centres[k]/1000:g} kHz')
 axs[1,1].set(xlim=(395,505),ylim=(-60,2),xlabel='RF frequency (kHz)',ylabel='Normalized Tx power (dB)',title='(d) Actual six-band Tx spectra');axs[1,1].legend(fontsize=9,ncol=2,loc='lower center',framealpha=.85)
 for ax in axs.flat:ax.grid(alpha=.15)
 fig.tight_layout(pad=.7,h_pad=1.1,w_pad=.9);save(fig,'actual_finite_pulse_ambiguity_scope')

spatial=run/'two_way/two_way_representative_profiles.npz'
if spatial.exists():
 register(spatial);data=np.load(spatial);angle=data['target_angles_deg'];records=rows(run/'two_way/two_way_approximation_comparison.csv');fig,axs=plt.subplots(1,3,figsize=(6,2.65));scan=ws.CFG.focus[-1]
 for ax,method in zip(axs[:2],('cdm','fdm')):
  for model,ls in [('true_delay','-'),('centre_approx','--')]:
   x=data[f'{method}_scan{scan:g}_all_simultaneous_{model}'];ax.plot(angle,ws.dB(x),color=ws.COLORS[method],ls=ls,lw=1.3,label='Full model' if model=='true_delay' else 'Centre approx.')
  ax.set(xlim=(-50,50),ylim=(-60,2),xlabel='Target bearing (deg)',ylabel='Beam power (dB)',title=('(a) ' if method=='cdm' else '(b) ')+ws.LABELS[method]);ax.legend(fontsize=9,loc='lower left',handlelength=1.2,handletextpad=.45,borderpad=.35);ax.grid(alpha=.15)
 ymin=1.
 for method in('cdm','fdm'):
  r=[r for r in records if r['method']==method and r['scope']=='all_simultaneous' and float(r['receive_scan_angle_deg']) in ws.CFG.focus];v=[float(x['magnitude_profile_correlation']) for x in r];ymin=min(ymin,min(v));axs[2].plot([float(x['receive_scan_angle_deg']) for x in r],v,color=ws.COLORS[method],ls=ws.STYLES[method],marker='s' if method=='cdm' else '^',markersize=3.5,lw=1.3,label=ws.LABELS[method])
 axs[2].set(xlabel='Scan bearing (deg)',ylabel='Magnitude-profile correlation',ylim=(max(0,ymin-.005),1.001),title='(c) Six sector centres');axs[2].legend(fontsize=9,loc='lower left');axs[2].grid(alpha=.15)
 fig.tight_layout(pad=.7,w_pad=.8);save(fig,'complete_two_way_spatial_check')
ws.dump(out/'plot_only_receipt.json',{**manifest,'pdfs':[p.name for p in out.glob('*.pdf')],'raw_metrics_unchanged':True})
print(json.dumps({'out':str(out),'pdfs':[p.name for p in out.glob('*.pdf')]}),flush=True)
