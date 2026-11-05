from __future__ import annotations
import csv
import argparse
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
ROOT=Path(__file__).resolve().parents[3]
parser=argparse.ArgumentParser(description='Plot historical MATLAB A16 outputs, not the v08 Figure 14 statistics.')
parser.add_argument('--results',type=Path,default=ROOT/'outputs/matlab/A16_paper')
RESULTS=parser.parse_args().results;FIGURES=RESULTS/'figures'
def read(n):
 with (RESULTS/n).open(newline='',encoding='utf-8-sig') as f:return list(csv.DictReader(f))
raw=read('A16_all_trials.csv');sens=read('A16_leakage_sensitivity.csv')
FIGURES.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'serif','font.serif':['Times New Roman','Times','DejaVu Serif'],'font.size':8,'axes.labelsize':8,'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':7,'axes.linewidth':.7,'savefig.dpi':600})
methods=[('dense','Dense ULA','#222222'),('cdm','CDM SW-ULA','#D55E00'),('fdm_current','packed FDM SW-ULA','#0072B2')]
fig,axes=plt.subplots(1,3,figsize=(7.15,2.55));metrics=[('Correlation2DToDense','2D correlation to dense'),('RangeProfileCorrelationToDense','Range-profile correlation'),('AngleProfileCorrelationToDense','Angle-profile correlation')]
for ax,(field,title) in zip(axes[:2],metrics[:2]):
 data=[np.array([float(r[field]) for r in raw if r['Method']==m]) for m,_,_ in methods];bp=ax.boxplot(data,positions=np.arange(3),widths=.55,patch_artist=True,showfliers=False);[patch.set_facecolor(methods[i][2]+'55') for i,patch in enumerate(bp['boxes'])];ax.set_xticks(np.arange(3),[x[1].replace(' SW-ULA','') for x in methods],rotation=15);ax.set_ylabel(title);ax.set_title(title)
for method,label,color in methods:
 b=[r for r in sens if r['Method']==method];guards=sorted({(int(float(r['GuardRangeBins'])),int(float(r['GuardAngleBins']))) for r in b});x=np.arange(len(guards));mean=[];ci=[]
 for g in guards:
  v=np.array([float(r['SupportLeakage']) for r in b if int(float(r['GuardRangeBins']))==g[0] and int(float(r['GuardAngleBins']))==g[1]]);mean.append(v.mean());ci.append(1.96*v.std(ddof=1)/np.sqrt(len(v)))
 axes[2].errorbar(x,mean,yerr=ci,color=color,marker='o',markersize=3,label=label)
axes[2].set_xticks(np.arange(len(guards)),[f'{g[0]}/{g[1]}' for g in guards]);axes[2].set(xlabel='Guard radius (range/angle bins)',ylabel='Support leakage',title='Mask sensitivity');axes[2].legend(frameon=False)
for i,ax in enumerate(axes):ax.grid(True,color='#D9D9D9',linewidth=.45);ax.spines[['top','right']].set_visible(False);ax.text(.02,.96,f'({chr(97+i)})',transform=ax.transAxes,va='top',fontweight='bold')
fig.tight_layout(pad=.7,w_pad=.9);fig.savefig(FIGURES/'A16_shipwreck_phase_ensemble.pdf',bbox_inches='tight');fig.savefig(FIGURES/'A16_shipwreck_phase_ensemble.png',bbox_inches='tight');print(FIGURES/'A16_shipwreck_phase_ensemble.pdf')
