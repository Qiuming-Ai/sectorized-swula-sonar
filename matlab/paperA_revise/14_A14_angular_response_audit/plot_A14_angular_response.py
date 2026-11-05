from __future__ import annotations
import csv
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
ROOT=Path(__file__).resolve().parent;RESULTS=ROOT/'results';FIGURES=ROOT/'figures';FIGURES.mkdir(exist_ok=True)
with (RESULTS/'A14_angular_profiles.csv').open(newline='',encoding='utf-8-sig') as f: rows=list(csv.DictReader(f))
plt.rcParams.update({'font.family':'serif','font.serif':['Times New Roman','Times','DejaVu Serif'],'font.size':8,'axes.labelsize':8,'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':7,'axes.linewidth':.7,'savefig.dpi':600})
fig,axes=plt.subplots(1,2,figsize=(7.15,2.6));styles={'dense':('#222222','-'),'sw_ula':('#0072B2','--'),'cdm':('#D55E00','--'),'fdm_current':('#0072B2','-.')}
def block(kind,method,target=8):return sorted([r for r in rows if r['ResponseType']==kind and r['Method']==method and abs(float(r['TargetAngleDeg'])-target)<1e-9],key=lambda r:float(r['ScanAngleDeg']))
for method,label in [('dense','Dense ULA'),('sw_ula','SW-ULA (CDM/FDM common Rx)')]:
 b=block('pure_receive',method);x=np.array([float(r['ScanAngleDeg']) for r in b]);y=np.maximum(-55,np.array([float(r['MagnitudeDb']) for r in b]));c,ls=styles[method];axes[0].plot(x,y,color=c,linestyle=ls,linewidth=1.1,label=label)
for method,label in [('dense','Dense ULA'),('cdm','CDM effective PSF'),('fdm_current','packed FDM effective PSF')]:
 b=block('effective_angular_psf',method);x=np.array([float(r['ScanAngleDeg']) for r in b]);y=np.maximum(-55,np.array([float(r['MagnitudeDb']) for r in b]));c,ls=styles[method];axes[1].plot(x,y,color=c,linestyle=ls,linewidth=1.1,label=label)
for i,(ax,title) in enumerate(zip(axes,['Pure receive response','Transmit–receive effective angular PSF'])):
 ax.set(xlabel='Scan angle (deg)',ylabel='Normalized magnitude (dB)',xlim=(-50,50),ylim=(-55,2),title=title);ax.grid(True,color='#D9D9D9',linewidth=.45);ax.spines[['top','right']].set_visible(False);ax.legend(frameon=True,facecolor='white',edgecolor='none',loc='upper left',bbox_to_anchor=(.06,.88));ax.text(.02,.96,f'({chr(97+i)})',transform=ax.transAxes,va='top',fontweight='bold')
fig.tight_layout(pad=.7,w_pad=1);fig.savefig(FIGURES/'A14_angular_response.pdf',bbox_inches='tight');fig.savefig(FIGURES/'A14_angular_response.png',bbox_inches='tight');print(FIGURES/'A14_angular_response.pdf')
