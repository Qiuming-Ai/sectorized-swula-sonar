"""Redraw manuscript Figures 10 and 11 from bundled, read-only data.

Run from any working directory. Outputs are written under outputs/figures;
no echo simulation, metric recomputation, or shipwreck data access occurs.
The shipwreck figures are outside this public replotter's scope.
"""
from pathlib import Path
import csv,json,hashlib,argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

SCRIPT=Path(__file__).resolve()
ROOT=SCRIPT.parents[2]
DATA=ROOT/'data'
OUT=ROOT/'outputs'/'figures'
POINT_MAIN=DATA/'point30m/point30m_main_minus15dB_arrays.npz'
POINT_CLEAN=DATA/'point30m/clean_images.npz'
POINT_STATS=DATA/'point30m/summary_mean_ci.csv'
DENSE_STATS=DATA/'dense1000x50/summary_mean_ci.csv'
INPUTS={'point':(POINT_MAIN,POINT_CLEAN),'noise':(POINT_STATS,DENSE_STATS)}
METHODS=('dense','cdm','fdm')
LABELS={'dense':'Dense ULA','cdm':'CDM SW-ULA','fdm':'FDM SW-ULA'}
COLORS={'dense':'#0072B2','cdm':'#D55E00','fdm':'#009E73'}
STYLES={'dense':('-', 'o'),'cdm':('--','s'),'fdm':(':','^')}
plt.rcParams.update({'font.family':'serif','font.serif':['Times New Roman','DejaVu Serif'],
 'font.size':10,'axes.labelsize':10,'xtick.labelsize':9.6,'ytick.labelsize':9.6,'legend.fontsize':9.6,
 'pdf.fonttype':42,'ps.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def output_path(name):
    path=OUT/name
    if path.resolve().parent!=OUT:
        raise RuntimeError(f'Refusing redirected figure output: {path}')
    return path
def write_json(name,data):
    output_path(name).write_text(json.dumps(data,indent=2)+'\n',encoding='utf-8')
def save(fig,name):
    for ext in ('pdf','png','svg'):fig.savefig(output_path(name+'.'+ext),dpi=600,bbox_inches='tight')
    plt.close(fig)
def rows(path):
    with Path(path).open(encoding='utf-8-sig') as f:return list(csv.DictReader(f))
def level(x,ref=None):
    return np.maximum(20*np.log10(np.maximum(abs(x),1e-30)/(float(abs(x).max()) if ref is None else ref)),-50)

def point():
    axis_keys=('ranges_m','angles_deg','true_ranges_m','true_angles_deg')
    with np.load(POINT_MAIN,allow_pickle=False) as main, np.load(POINT_CLEAN,allow_pickle=False) as quiet:
        if any(not np.array_equal(main[k],quiet[k]) for k in axis_keys):
            raise ValueError('Saved noisy and clean point grids or target coordinates differ')
        r,a,tr,ta=(main[k] for k in axis_keys)
        noisy={m:main[m] for m in METHODS};clean={m:quiet[m] for m in METHODS}
    if any(z.shape!=(len(r),len(a)) or not np.isfinite(z).all()
           for images in (noisy,clean) for z in images.values()):
        raise ValueError('Saved point image does not match the declared finite range-angle grid')
    for prefix,images in [('point30m_clean',clean),('point30m',noisy)]:
        ref=max(abs(images[m]).max() for m in METHODS)
        for common in (False,True):
            fig,axes=plt.subplots(2,3,figsize=(7.15,4.9),sharex='row',sharey='row')
            for j,m in enumerate(METHODS):
                z=level(images[m],ref if common else None)
                axes[0,j].pcolormesh(a,r,z,shading='auto',cmap='cividis',vmin=-50,vmax=0,rasterized=True)
                axes[0,j].scatter(ta,tr,s=12,facecolor='none',edgecolor='white',linewidth=.6)
                axes[0,j].add_patch(Rectangle((-2.5,30.065),5,.115,fill=False,edgecolor='white',ls='--',lw=.7))
                ir=(r>=30.065)&(r<=30.18);ia=(a>=-2.5)&(a<=2.5)
                im=axes[1,j].pcolormesh(a[ia],r[ir],z[np.ix_(ir,ia)],shading='auto',cmap='cividis',vmin=-50,vmax=0,rasterized=True)
                selected=(tr>=30.065)&(tr<=30.18)&(abs(ta)<=2.5)
                axes[1,j].scatter(ta[selected],tr[selected],s=15,facecolor='none',edgecolor='white',linewidth=.7)
                axes[0,j].set_title(f'({chr(97+j)}) {LABELS[m]}',fontsize=10)
                axes[0,j].set_xlim(-50,50);axes[0,j].set_ylim(30,30.8)
                axes[1,j].set(xlabel='Angle (deg)',xlim=(-2.5,2.5),ylim=(30.065,30.18))
                axes[1,j].text(.03,.04,'5° × 0.115 m',transform=axes[1,j].transAxes,color='white',fontsize=9.6)
            axes[0,0].set_ylabel('Range (m)');axes[1,0].set_ylabel('Range (m)')
            fig.subplots_adjust(left=.075,right=.89,bottom=.10,top=.95,wspace=.18,hspace=.18)
            cb=fig.add_axes([.91,.15,.017,.72]);fig.colorbar(im,cax=cb,label='Amplitude (dB)')
            save(fig,prefix+('_common_scale' if common else '_per_image'))
    write_json('point_main_display_protocol.json',{'noise_normalization':'common_per_element_noise',
        'nominal_Dense_raw_input_snr_db':-15,'noise_repetition':0,'noise_seed':20261015200,
        'source_array':POINT_MAIN.relative_to(ROOT).as_posix(),
        'common_scale':'one scalar maximum across three numerical output images',
        'display_range_db':[-50,0],'zoom_deg_m':[5,.115],
        'clean_outputs':'point30m_clean_* and point30m_axis_profiles are supplementary'})
    fig,axes=plt.subplots(1,2,figsize=(7.15,2.55))
    for m in METHODS:
        ls,mark=STYLES[m]
        axes[0].plot(a,level(clean[m][np.argmin(abs(r-30.22))]),ls=ls,color=COLORS[m],label=LABELS[m])
        axes[1].plot((r-30.22)*1000,level(clean[m][:,np.argmin(abs(a))]),ls=ls,color=COLORS[m],label=LABELS[m])
    axes[0].set(xlabel='Angle (deg), at r = 30.22 m',ylabel='Normalized magnitude (dB)',xlim=(-50,50),ylim=(-50,2))
    axes[1].set(xlabel='Range offset from 30.22 m (mm)',ylabel='Normalized magnitude (dB)',xlim=(-160,420),ylim=(-50,2))
    for ax in axes:ax.legend();ax.grid(alpha=.15)
    fig.tight_layout();save(fig,'point30m_axis_profiles')

def noise():
    pointrows=rows(POINT_STATS);dense=rows(DENSE_STATS)
    for rule in ('common_per_element_noise','equal_raw_channel_snr'):
        fig,axes=plt.subplots(2,2,figsize=(7.15,5.6))
        keys=[('cnr_linear_magnitude','Point-scene CNR'),('noise_induced_magnitude_NMSE_db','Magnitude perturbation NMSE (dB)'),
              ('dense_magnitude_correlation','Noise-free Dense correlation'),('truth_magnitude_correlation','Scatterer-support magnitude correlation')]
        for panel,ax in enumerate(axes.flat):
            metric,label=keys[panel]
            for scene_index,data in enumerate((pointrows,dense)):
                if panel==0 and scene_index==1:continue
                for m in METHODS:
                    q=sorted([r for r in data if r['method']==m and r['noise_normalization']==rule],key=lambda r:float(r['nominal_snr_db']))
                    ls,mark=STYLES[m]
                    ax.errorbar([float(r['nominal_snr_db']) for r in q],[float(r[metric]) for r in q],
                        yerr=[float(r[metric+'_ci95_halfwidth'] or 0) for r in q],
                        ls='-' if scene_index==0 else '--',marker=mark,markersize=4,
                        markerfacecolor=COLORS[m] if scene_index==0 else 'white',
                        color=COLORS[m],capsize=2,linewidth=1.1,label=LABELS[m] if scene_index==0 else None)
            ax.set(xlabel='Nominal raw input SNR (dB)',ylabel=label);ax.set_xticks([-15,-10,-5,0]);ax.grid(alpha=.15)
            ax.text(.03,.95,f'({chr(97+panel)})',transform=ax.transAxes,va='top',fontsize=10)
            if panel==0:ax.legend(fontsize=9.6)
        axes[0,1].plot([],[],color='.35',ls='-',marker='o',label='18 points')
        axes[0,1].plot([],[],color='.35',ls='--',marker='o',markerfacecolor='white',label='1000 points / 50 layers')
        handles,labels=axes[0,1].get_legend_handles_labels()
        scene_keys=[(h,l) for h,l in zip(handles,labels) if l in ('18 points','1000 points / 50 layers')]
        axes[0,1].legend([h for h,l in scene_keys],[l for h,l in scene_keys],fontsize=9.6)
        fig.tight_layout();save(fig,'channel_mc_'+rule)

def main():
    ship_status='For Figure 13, extract the shipwreck data asset and run code/scenes/plot_shipwreck.py.'
    parser=argparse.ArgumentParser(description=__doc__,epilog=ship_status)
    parser.add_argument('--only',nargs='+',choices=['point','noise','ship'],
                        help='Select point (Figure 10) and/or noise (Figure 11); default: both.')
    args=parser.parse_args();selected=list(dict.fromkeys(args.only or ['point','noise']))
    if 'ship' in selected:parser.error(ship_status)
    required=[path for group in selected for path in INPUTS[group]]
    missing=[path.relative_to(ROOT).as_posix() for path in required if not path.is_file()]
    if missing:parser.error('Missing bundled input(s): '+', '.join(missing))
    if OUT.resolve()!=OUT:raise RuntimeError(f'Refusing redirected figure output directory: {OUT}')
    OUT.mkdir(parents=True,exist_ok=True)
    hashes={path.relative_to(ROOT).as_posix():sha(path) for path in required}
    if 'point' in selected:point()
    if 'noise' in selected:noise()
    if any(sha(ROOT/name)!=digest for name,digest in hashes.items()):
        raise RuntimeError('A bundled input changed while plotting')
    main_files={'point':'point30m_common_scale.pdf','noise':'channel_mc_common_per_element_noise.pdf'}
    write_json('figure_protocol_and_audit.json',{'selected':selected,
        'protocol':'saved arrays and statistics only; no echo or metric recomputation',
        'plotter_sha256':sha(SCRIPT),'input_sha256':hashes,'inputs_unchanged':True,
        'main_outputs':{group:main_files[group] for group in selected},'shipwreck_status':ship_status,
        'minimum_font_pt':9.6,'PNGs_dpi':600,'heatmaps':'rasterized data, vector axes/labels in PDF and SVG'})
    print('Replotted '+', '.join(selected)+' into '+str(OUT))
if __name__=='__main__':main()
