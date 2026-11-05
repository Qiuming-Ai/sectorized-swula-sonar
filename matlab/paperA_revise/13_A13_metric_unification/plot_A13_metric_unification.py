from __future__ import annotations
import csv
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parent
RESULTS=ROOT/'results'; FIGURES=ROOT/'figures'; FIGURES.mkdir(exist_ok=True)
with (RESULTS/'A13_summary_mean_ci.csv').open(newline='',encoding='utf-8-sig') as f: rows=list(csv.DictReader(f))
plt.rcParams.update({'font.family':'serif','font.serif':['Times New Roman','Times','DejaVu Serif'],'font.size':8,'axes.labelsize':8,'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':7,'axes.linewidth':0.7,'savefig.dpi':600})
methods=[('dense','Dense ULA','#222222','o'),('cdm','CDM SW-ULA','#D55E00','s'),('fdm_current','packed FDM SW-ULA','#0072B2','^')]
fig,axes=plt.subplots(2,2,figsize=(7.15,4.65))
panels=[('point18','MeanCNR','CI95CNR','Point-scene CNR','CNR'),('point18','MeanNoiseNMSEdB','CI95NoiseNMSEdB','Point-scene noise sensitivity','Noise-induced NMSE (dB)'),('dense1000x50','MeanCorrelation2DToDenseReference','CI95Correlation2DToDenseReference','Dense-scene agreement with dense reference','Correlation to dense reference'),('dense1000x50','MeanSupportLeakage','CI95SupportLeakage','Dense-scene support leakage','Energy fraction')]
for idx,(scenario,yfield,efield,title,ylabel) in enumerate(panels):
 ax=axes.flat[idx]
 for method,label,color,marker in methods:
  block=sorted([r for r in rows if r['Scenario']==scenario and r['Method']==method],key=lambda r:float(r['InputSNRDb']))
  x=np.array([float(r['InputSNRDb']) for r in block]); y=np.array([float(r[yfield]) for r in block]); e=np.array([float(r[efield]) for r in block])
  ax.errorbar(x,y,yerr=e,color=color,marker=marker,linewidth=1,markersize=3.5,capsize=2,label=label)
 ax.set(xlabel='Raw-channel input SNR (dB)',ylabel=ylabel,title=title); ax.grid(True,color='#D9D9D9',linewidth=.45); ax.spines[['top','right']].set_visible(False); ax.text(.02,.96,f'({chr(97+idx)})',transform=ax.transAxes,va='top',fontweight='bold')
axes[0,0].legend(frameon=False); fig.tight_layout(pad=.8,w_pad=1,h_pad=1)
fig.savefig(FIGURES/'A13_metric_unification.pdf',bbox_inches='tight'); fig.savefig(FIGURES/'A13_metric_unification.png',bbox_inches='tight')
print(FIGURES/'A13_metric_unification.pdf')
