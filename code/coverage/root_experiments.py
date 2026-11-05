"""Spatial audits and rerun sensitivity tests; all use the locked common model."""
from numerical_model import *
import legacy_experiments as legacy
from dataclasses import replace

def coverage(out):
    legacy.sector_search(out)
    cfg=Config();rows=[]
    for astep,fstep in [(.01,1000),(.005,500)]:
        theta=np.unique(np.r_[np.arange(-50,50+astep/2,astep),cfg.edges]);worst_cover=1.;worst_gl=0.;atc=atg=None
        for f in np.arange(400000,500001,fstep):
            tx=abs(tx_pattern(cfg,[f],theta)[0])
            for k in range(6):
                in_sector=(theta>=cfg.edges[k])&(theta<=cfg.edges[k+1]);t=theta[in_sector];response=tx[k,in_sector];cov=response/tx[k].max();idx=int(cov.argmin())
                if cov[idx]<worst_cover:worst_cover=float(cov[idx]);atc=(float(t[idx]),int(f),k+1)
                u=np.sin(np.deg2rad(t))
                for order in [-4,-3,-2,-1,1,2,3,4]:
                    ug=u+order*cfg.fc/(2*f);visible=(abs(ug)<=1)&(ug>=np.sin(np.deg2rad(-50)))&(ug<=np.sin(np.deg2rad(50)))
                    if not visible.any():continue
                    alias=np.rad2deg(np.arcsin(ug[visible]));ar=abs(tx_pattern(cfg,[f],alias,centres=[cfg.focus[k]])[0,0]);ratio=ar/np.maximum(response[visible],1e-30);idx=int(ratio.argmax())
                    if ratio[idx]>worst_gl:worst_gl=float(ratio[idx]);atg=(float(t[visible][idx]),float(alias[idx]),int(f),order,k+1)
        cdb=float(20*np.log10(worst_cover));gdb=float(20*np.log10(worst_gl))
        rows.append({'angular_step_deg':astep,'frequency_step_hz':fstep,'worst_coverage_db':cdb,'worst_grating_ratio_db':gdb,'coverage_margin_to_minus6p6_db':cdb+6.6,'grating_margin_to_minus24_db':-24-gdb,'coverage_angle_deg':atc[0],'coverage_frequency_hz':atc[1],'coverage_sector':atc[2],'grating_target_angle_deg':atg[0],'grating_alias_angle_deg':atg[1],'grating_frequency_hz':atg[2],'grating_order':atg[3],'grating_sector':atg[4],'nominal_feasible':bool(cdb>=-6.6 and gdb<=-24),'scope':'nominal RF400-500kHz, fixed implemented six sectors; not a calibration robustness claim'})
        print('Final coverage audit',astep,fstep,cdb,gdb,flush=True)
    write_csv(out/'grid_convergence_audit.csv',rows)

def boundary_2d(out):
    cfg=Config();targets=np.unique(np.r_[cfg.edges,np.clip(cfg.edges-.13,-50,50),np.clip(cfg.edges+.13,-50,50)])
    rr=np.arange(30,30.8001,.005);aa=np.arange(-50,50.0001,.1);rows=[];saved={}
    for a in targets:
        for m in METHODS:
            raw,mf,meta=point_raw(cfg,m,[30.5],[a],geometry='spherical');im=reconstruct_points(cfg,m,raw,mf,rr,aa)
            idx=np.unravel_index(abs(im).argmax(),im.shape);threshold=.1*float(abs(im).max());pd,fp,count=legacy._blind_detection(im,rr,aa,threshold,np.array([30.5]),np.array([a]))
            rows.append({'method':m,'truth_angle_deg':float(a),'truth_range_m':30.5,'range_error_mm':float((rr[idx[0]]-30.5)*1000),'angle_error_deg':float(aa[idx[1]]-a),'peak_amplitude':float(abs(im[idx])),'blind_detection_probability':pd,'false_peaks':fp,'total_blind_peaks':count,'threshold_amplitude':threshold,'threshold_scope':'noise-free relative-threshold boundary diagnostic, not calibrated FA'})
            if a in [-50,0,50]:saved[m+'_'+str(int(a))]=abs(im).astype(np.float32)
        print('Boundary 2D at30m',round(float(a),5),flush=True)
    write_csv(out/'metrics.csv',rows);np.savez_compressed(out/'representative_images.npz',ranges_m=rr,angles_deg=aa,**saved)
    legacy.curves(rows,'truth_angle_deg','false_peaks',out/'boundary_false_peaks','True angle (deg)','False peaks per noise-free image')

REGISTRY={'coverage':coverage,'beam_audit':legacy.beam_audit,'boundary_2d':boundary_2d,'calibration':legacy.calibration,'nonideal':legacy.nonideal,'wideband':legacy.wideband,'assignment':legacy.assignment,'zc_parameters':legacy.zc_parameters}
