from __future__ import annotations
import csv
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
ROOT=Path(__file__).resolve().parent;RESULTS=ROOT/'results';FIGURES=ROOT/'figures';FIGURES.mkdir(exist_ok=True)
def read(n):
 with (RESULTS/n).open(newline='',encoding='utf-8-sig') as f:return list(csv.DictReader(f))
m=read('A15_wideband_metrics.csv');p=read('A15_wideband_profiles.csv')
plt.rcParams.update({'font.family':'serif','font.serif':['Times New Roman','Times','DejaVu Serif'],'font.size':8,'axes.labelsize':8,'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':7,'axes.linewidth':.7,'savefig.dpi':600})
fig,axes=plt.subplots(1,3,figsize=(7.15,2.55));colors={'cdm':'#D55E00','fdm_current':'#0072B2'}
for ax,method,target,title in [(axes[0],'cdm',41.25,'CDM edge sector'),(axes[1],'fdm_current',41.25,'packed FDM edge sector')]:
 b=sorted([r for r in p if r['Method']==method and abs(float(r['TargetAngleDeg'])-target)<1e-9],key=lambda r:float(r['ScanAngleDeg']));x=np.array([float(r['ScanAngleDeg']) for r in b]);y0=np.maximum(-55,np.array([float(r['CentreFrequencyMagnitudeDb']) for r in b]));yw=np.maximum(-55,np.array([float(r['WaveformWeightedMagnitudeDb']) for r in b]));ax.plot(x,y0,color=colors[method],linestyle='--',label='centre-frequency');ax.plot(x,yw,color='#222222',linestyle='-',label='waveform-weighted delay');ax.set(xlabel='Scan angle (deg)',ylabel='Normalized magnitude (dB)',xlim=(20,50),ylim=(-55,2),title=title);ax.legend(frameon=True,facecolor='white',edgecolor='none',loc='upper left',bbox_to_anchor=(.03,.84))
for method,label,marker in [('cdm','CDM','s'),('fdm_current','packed FDM','^')]:
 b=sorted([r for r in m if r['Method']==method],key=lambda r:float(r['TargetAngleDeg']));x=np.array([float(r['TargetAngleDeg']) for r in b]);y=np.array([float(r['Correlation']) for r in b]);axes[2].plot(x,y,color=colors[method],marker=marker,label=label)
axes[2].set(xlabel='Target/sector angle (deg)',ylabel='Profile correlation',ylim=(0.97,1.001),title='Approximation agreement');axes[2].legend(frameon=False)
for i,ax in enumerate(axes):ax.grid(True,color='#D9D9D9',linewidth=.45);ax.spines[['top','right']].set_visible(False);ax.text(.02,.96,f'({chr(97+i)})',transform=ax.transAxes,va='top',fontweight='bold')
fig.tight_layout(pad=.7,w_pad=.9);fig.savefig(FIGURES/'A15_wideband_spatial_audit.pdf',bbox_inches='tight');fig.savefig(FIGURES/'A15_wideband_spatial_audit.png',bbox_inches='tight');print(FIGURES/'A15_wideband_spatial_audit.pdf')
