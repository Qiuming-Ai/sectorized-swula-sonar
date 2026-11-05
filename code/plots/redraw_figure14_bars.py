"""Compact three-row Figure 14 bars from unchanged saved statistics.

No acoustic experiment or reconstruction is rerun.  The root agent manages the
unique authoring marker and replacement of the canonical main-package asset.
"""
from pathlib import Path
import csv
import json
import hashlib
import numpy as np
from scipy.stats import t
import fitz
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager

ROOT=Path(__file__).resolve().parents[2]
WORK=ROOT/'outputs/figures/figure14'
WORK.mkdir(parents=True,exist_ok=True)
SOURCE=ROOT/'data/shipwreck'
AUDIT=SOURCE/'FIGURE14_SOURCE_AUDIT.json'
PDF=WORK/'figure14_ship_statistics.pdf'
PNG=WORK/'figure14_ship_statistics_preview.png'
METHODS=('dense','cdm','fdm')
LABELS=('Dense','CDM','FDM')
COLORS=('#0072B2','#D55E00','#009E73')
METRICS=('dense_magnitude_correlation','angle_profile_correlation_to_clean_Dense','support_leakage')
YLABELS=('2-D Dense\ncorrelation','Angular-profile\ncorrelation','Support leakage')
LIMITS=((0,1.05),(-.4,1.05),(0,.4))
TICKS=((0,.5,1),(-.4,0,.5,1),(0,.2,.4))
WIDTH=.52
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()

def selected_csv(path):
    with path.open(encoding='utf-8-sig',newline='') as stream:
        return [r for r in csv.DictReader(stream)if r['noise_normalization']=='common_per_element_noise'and float(r['nominal_snr_db'])==10]

prior=json.loads(AUDIT.read_text(encoding='utf-8'))
assert prior['summary_sha256']==sha(SOURCE/'summary_mean_ci.csv')
assert prior['trials_sha256']==sha(SOURCE/'all_trials.csv')
summary=selected_csv(SOURCE/'summary_mean_ci.csv')
trials=selected_csv(SOURCE/'all_trials.csv')
assert len(summary)==3 and len(trials)==90
data={m:next(r for r in summary if r['method']==m)for m in METHODS}
statistics=[]
for method in METHODS:
    assert int(data[method]['replications'])==30
    raw=[r for r in trials if r['method']==method]
    assert len(raw)==30
    for metric in METRICS:
        mean=float(data[method][metric]);half=float(data[method][metric+'_ci95_halfwidth'])
        values=np.array([float(r[metric])for r in raw])
        assert abs(mean-float(values.mean()))<1e-12
        assert abs(half-float(t.ppf(.975,29)*values.std(ddof=1)/np.sqrt(30)))<1e-12
        expected=next(r for r in prior['statistics']if r['method']==method and r['metric']==metric)
        assert mean==expected['mean']and half==expected['ci95_halfwidth']
        statistics.append({'method':method,'metric':metric,'n':30,'mean':mean,'ci95_halfwidth':half})
font_path=font_manager.findfont(font_manager.FontProperties(family=['Times New Roman','DejaVu Serif']))
plt.rcParams.update({
    'font.family':'serif','font.serif':['Times New Roman','DejaVu Serif'],'font.size':9,
    'axes.labelsize':9,'xtick.labelsize':9,'ytick.labelsize':9,
    'axes.unicode_minus':False,'axes.linewidth':.65,
    'axes.spines.top':False,'axes.spines.right':False,
    'pdf.fonttype':42,'ps.fonttype':42,'savefig.facecolor':'white',
})
fig,axes=plt.subplots(3,1,figsize=(85/25.4,82/25.4))
fig.subplots_adjust(left=.245,right=.97,bottom=.065,top=.95,hspace=.45)
expected_vectors=[]
for panel,(ax,metric,ylabel,ylim,ticks)in enumerate(zip(axes,METRICS,YLABELS,LIMITS,TICKS)):
    means=np.array([float(data[m][metric])for m in METHODS])
    cis=np.array([float(data[m][metric+'_ci95_halfwidth'])for m in METHODS])
    containers=ax.bar(range(3),means,bottom=0,width=WIDTH,color=COLORS,
        edgecolor='#303030',linewidth=.6,zorder=3)
    ax.errorbar(range(3),means,yerr=cis,fmt='none',color='#202020',
        capsize=2.6,capthick=.7,elinewidth=.7,zorder=4)
    assert all(rect.get_y()==0 for rect in containers)
    assert all(rect.get_height()==means[i]for i,rect in enumerate(containers))
    ax.set(xlim=(-.55,2.55),ylim=ylim)
    ax.set_xticks(range(3),LABELS)
    ax.set_yticks(ticks)
    ax.set_ylabel(ylabel,labelpad=3)
    ax.set_title(f'({chr(97+panel)})',loc='left',fontsize=9,pad=2)
    ax.tick_params(direction='out',length=2.5,width=.65,pad=2.2)
    ax.grid(axis='y',alpha=.18,linewidth=.5,zorder=0)
    if panel==1:ax.axhline(0,color='#606060',linewidth=.65,zorder=2)
    for j,(method,mean,ci)in enumerate(zip(METHODS,means,cis)):
        # Convert Matplotlib's data transform to PDF coordinates (origin top left).
        # This independent vector-level QA detects wrong heights/baselines at export.
        x0,y0=ax.transData.transform((j-WIDTH/2,0))*72/fig.dpi
        x1,y1=ax.transData.transform((j+WIDTH/2,mean))*72/fig.dpi
        xc,yl=ax.transData.transform((j,mean-ci))*72/fig.dpi
        _,yu=ax.transData.transform((j,mean+ci))*72/fig.dpi
        height=fig.get_figheight()*72
        expected_vectors.append({'panel':panel,'method':method,'metric':metric,'mean':float(mean),
            'ci95_halfwidth':float(ci),'color_rgb':list(matplotlib.colors.to_rgb(COLORS[j])),
            'bar_rect_pt':[float(x0),float(min(height-y0,height-y1)),float(x1),float(max(height-y0,height-y1))],
            'ci_rect_pt':[float(xc),float(height-yu),float(xc),float(height-yl)]})
fig.savefig(PDF,format='pdf',metadata={'Title':'Figure 14: Shipwreck image statistics',
    'Subject':'Unchanged means and 95% confidence intervals shown as zero-baseline bars',
    'Creator':'Reproducible Matplotlib bar plot from frozen summary and trial CSV'})
plt.close(fig)

doc=fitz.open(PDF);page=doc[0]
page.get_pixmap(matrix=fitz.Matrix(4,4),alpha=False).save(PNG)
spans=[s for b in page.get_text('dict')['blocks']if 'lines'in b for l in b['lines']for s in l['spans']]
assert len(doc)==1 and len(page.get_images())==0
assert min(s['size']for s in spans)>=8
assert all(page.rect.contains(fitz.Rect(s['bbox']))for s in spans)
assert all(page.get_text().splitlines().count(label)==3 for label in LABELS)
drawings=page.get_drawings()
checks=[]
for expected in expected_vectors:
    bars=[d for d in drawings if d['fill']is not None and np.max(abs(np.array(d['fill'])-expected['color_rgb']))<1e-6]
    bar=min(bars,key=lambda d:np.max(abs(np.array(d['rect'])-expected['bar_rect_pt'])))
    bar_error=float(np.max(abs(np.array(bar['rect'])-expected['bar_rect_pt'])))
    vertical=[d for d in drawings if d['type']=='s'and len(d['items'])==1 and d['items'][0][0]=='l'and abs(d['rect'].x1-d['rect'].x0)<1e-5]
    ci=min(vertical,key=lambda d:np.max(abs(np.array(d['rect'])-expected['ci_rect_pt'])))
    ci_error=float(np.max(abs(np.array(ci['rect'])-expected['ci_rect_pt'])))
    assert bar_error<1e-4 and ci_error<1e-4,(expected,bar_error,ci_error)
    checks.append({'method':expected['method'],'metric':expected['metric'],
                   'zero_baseline_mean_vector_error_pt':bar_error,'ci_endpoint_vector_error_pt':ci_error})
result={
    'status':'authored_vector_checked_pending_visual_review',
    'input_summary':str(SOURCE/'summary_mean_ci.csv'),'input_summary_sha256':sha(SOURCE/'summary_mean_ci.csv'),
    'input_trials':str(SOURCE/'all_trials.csv'),'input_trials_sha256':sha(SOURCE/'all_trials.csv'),
    'statistics':statistics,'statistics_exactly_equal_to_previous_validated_figure':True,
    'recalculation_scope':'saved trial mean and Student-t CI audit only; no scene, image or experiment rerun',
    'canvas_mm':[85,82],'panels':3,'bar_baseline_all_panels':0,'bar_width_categorical_units':WIDTH,
    'method_order':list(METHODS),'axis_limits':[list(x)for x in LIMITS],
    'zero_reference_line_in_negative_panel':True,'negative_cdm_angular_result_preserved':True,
    'bars_and_ci_verified_from_saved_pdf_vectors':checks,
    'minimum_pdf_font_pt':min(s['size']for s in spans),'maximum_pdf_font_pt':max(s['size']for s in spans),
    'embedded_font_file':font_path,'embedded_pdf_fonts':page.get_fonts(full=True),
    'vector_path_count':len(drawings),'bitmap_image_count':len(page.get_images()),'text_bounds_inside_canvas':True,
    'no_added_value_annotations':True,'tiny_ci_not_artificially_enlarged':True,
    'output_pdf':str(PDF),'preview_png':str(PNG),'pdf_sha256':sha(PDF),'png_sha256':sha(PNG),
}
(WORK/'FIGURE14_BAR_QA.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(result,ensure_ascii=False,indent=2))
