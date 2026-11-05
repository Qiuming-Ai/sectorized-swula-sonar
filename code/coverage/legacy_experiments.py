"""Scientific experiments, separated from wording and reviewer certification."""
from dataclasses import asdict,replace
from pathlib import Path
import itertools,json,time,csv
import numpy as np
from scipy import signal,ndimage
import matplotlib.pyplot as plt
from numerical_model import *

# Explicit 50m diagnostic origin for the inherited sensitivity test coordinates.
# Main point/noise experiments are rebuilt separately at 30m. Old files remain untouched.
_BaseConfig=Config
def Config(*args,**kwargs):
    kwargs.setdefault('range_origin',50.0)
    return _BaseConfig(*args,**kwargs)

def curves(rows,x,y,path,xlabel,ylabel,group='method'):
    fig,ax=plt.subplots(figsize=(5,3.6));styles=['o-','s--','^-.','d:']
    groups=list(dict.fromkeys(r[group] for r in rows))
    for i,g in enumerate(groups):
        rr=[r for r in rows if r[group]==g and r.get(y) is not None]
        if not rr:continue
        rr=sorted(rr,key=lambda r:r[x]);ax.plot([r[x] for r in rr],[r[y] for r in rr],styles[i%4],label=LABELS.get(g,g),color=COLORS.get(g,None),markersize=4)
    ax.set(xlabel=xlabel,ylabel=ylabel);ax.grid(alpha=.2);ax.legend();savefig(fig,path)

def waveform(out):
    cfg=Config();rows=[];spectral=[];nfft=65536;freq=np.fft.fftshift(np.fft.fftfreq(nfft,1/cfg.fs));lag_store={}
    fig,axs=plt.subplots(1,2,figsize=(10,3.7))
    for method in METHODS:
        for taper in ['rectangular','hamming']:
            for shape in [False,True]:
                w,meta=waveforms(cfg,method,taper,shape);tag=f'{method}_{taper}_'+('FIR' if shape else 'original_sampling')
                p=np.abs(np.fft.fftshift(np.fft.fft(w,nfft,axis=0),axes=0))**2
                c=signal.correlate(w[:,0],w[:,0],mode='full',method='fft')
                ci=signal.resample(c,8*len(c));lag=(np.arange(len(ci))/8-(len(w)-1))/cfg.fs
                width=peak_width(lag*cfg.c/2,np.abs(ci)**2)
                cross=[];cross_window=[]
                for i,j in itertools.combinations(range(w.shape[1]),2):
                    cc=signal.correlate(w[:,i],w[:,j],mode='full',method='fft');ll=signal.correlation_lags(len(w),len(w))/cfg.fs
                    cross.append(float(np.max(np.abs(cc))))
                    cross_window.append(float(np.max(np.abs(cc[np.abs(ll)<=.0006]))))
                for k in range(w.shape[1]):
                    s=p[:,k];cdf=np.cumsum(s)/s.sum();lo=freq[np.searchsorted(cdf,.005)];hi=freq[min(np.searchsorted(cdf,.995),nfft-1)]
                    spectral.append({'configuration':tag,'channel':k+1,'outside_nominal_100khz_fraction':float(s[np.abs(freq)>50000].sum()/s.sum()),'equal_tail_99pct_bandwidth_hz':float(hi-lo),'lower_99pct_hz':float(lo),'upper_99pct_hz':float(hi),'support_ms':len(w)/cfg.fs*1000})
                rows.append({'configuration':tag,'method':method,'taper':taper,'shaping':shape,'support_ms':len(w)/cfg.fs*1000,'power_fwhm_mm':None if width is None else width*1000,'max_pair_cross_full_db':None if not cross else float(20*np.log10(max(cross))),'max_pair_cross_window_db':None if not cross_window else float(20*np.log10(max(cross_window))),'outside_nominal_band_worst':max(x['outside_nominal_100khz_fraction'] for x in spectral if x['configuration']==tag)})
                lag_store[tag]=np.column_stack([lag,ci.real,ci.imag])
                if taper==('rectangular' if method=='cdm' else 'hamming') and not shape:
                    ls={'dense':'-','cdm':'--','fdm':'-.'}[method]
                    axs[0].plot(freq/1000,10*np.log10(np.maximum(p[:,0]/p[:,0].max(),1e-12)),label=LABELS[method],color=COLORS[method],ls=ls)
                    ix=np.abs(lag*cfg.c/2)<=.15;axs[1].plot(lag[ix]*cfg.c/2*1000,db(ci[ix]),label=LABELS[method],color=COLORS[method],ls=ls)
    axs[0].set(xlabel='Baseband frequency (kHz)',ylabel='Relative power (dB)',ylim=(-80,3));axs[0].axvline(-50,color='k',ls=':');axs[0].axvline(50,color='k',ls=':')
    axs[1].set(xlabel='Range lag (mm)',ylabel='Matched-filter amplitude (dB)',ylim=(-60,3));axs[1].legend();savefig(fig,out/'original_waveforms')
    write_csv(out/'metrics.csv',rows);write_csv(out/'all_channels_spectrum.csv',spectral);np.savez_compressed(out/'range_correlation_raw.npz',**lag_store)
    # Fractional-delay and interpolation diagnostics are measurements, not labels.
    cfgp=Config(interpolation_oversampling=8);w,_=waveforms(cfgp,'dense');base=np.fft.fft(w[:,0],8192)
    checks=[]
    for delay in [0.,.125,.37,.75,1.33]:
        shifted=np.fft.ifft(base*np.exp(-2j*np.pi*np.fft.fftfreq(8192)*delay));mf=np.fft.ifft(np.fft.fft(shifted)*base.conj());fine=signal.resample(mf,8192*8)
        index=delay*8;i=int(np.floor(index));v=fine[i]*(1-index+i)+fine[i+1]*(index-i)
        exact=np.sum(np.abs(base)**2)/8192
        checks.append({'delay_samples':delay,'interpolated_peak_magnitude':float(abs(v)),'exact_spectral_peak':float(exact),'relative_peak_error':float(abs(abs(v)-exact)/exact)})
    write_csv(out/'fractional_delay_checks.csv',checks)
    dump(out/'scope.json',{'waveform_rows':len(rows),'all_channel_rows':len(spectral),'FIR_control':'same 129-tap FIR applied to every method; full transients retained, pulse support is longer and explicitly reported','bandwidth_definition':'equal-tail 99% of discrete sampled spectrum; not analogue infinite-band chip spectrum','full_delay_cross':'all finite aperiodic lags','window_delay_cross':'abs(tau)<=0.6ms','FWHM':'half maximum of power, Fourier-interpolated correlation','configuration':asdict(cfg)})

def nearfield(out):
    rows=[];store={};angles=np.linspace(-50,50,401)
    for distance in [5.,15.,30.,60.]:
        cfg=Config(range_origin=distance-.5,interpolation_oversampling=8)
        rr=np.linspace(distance-.12,distance+.12,121)
        for target_angle in [0.,40.]:
            for method in METHODS:
                raw,mf,meta=point_raw(cfg,method,[distance],[target_angle],geometry='spherical')
                reference=None
                for focus in ['spherical','farfield']:
                    im=reconstruct_points(cfg,method,raw,mf,rr,angles,focus)
                    ii=np.unravel_index(np.abs(im).argmax(),im.shape);peak=float(abs(im[ii]))
                    if focus=='spherical':reference=peak
                    local=np.abs(angles-target_angle)<=8
                    width=peak_width(angles[local],np.max(np.abs(im[:,local])**2,axis=0))
                    outside=np.abs(angles-target_angle)>8
                    ghost=float(20*np.log10(max(np.max(np.abs(im[:,outside])),1e-30)/peak))
                    rows.append({'method':method,'distance_m':distance,'target_angle_deg':target_angle,'reconstruction':focus,'angle_estimate_deg':float(angles[ii[1]]),'range_estimate_m':float(rr[ii[0]]),'angle_error_deg':float(angles[ii[1]]-target_angle),'range_error_mm':float((rr[ii[0]]-distance)*1000),'peak':peak,'loss_to_spherical_db':float(20*np.log10(peak/reference)),'local_power_hpbw_deg':width,'outside_8deg_peak_db':ghost})
                    if distance==5 and target_angle==0:store[method+'_'+focus]=np.abs(im).astype(np.float32)
    write_csv(out/'metrics.csv',rows);np.savez_compressed(out/'five_metre_images.npz',angles=angles,ranges=np.linspace(4.88,5.12,121),**store)
    select=[x for x in rows if x['target_angle_deg']==0 and x['reconstruction']=='farfield']
    curves(select,'distance_m','loss_to_spherical_db',out/'nearfield_loss','Target range (m)','Loss relative to spherical focusing (dB)')
    dump(out/'model_scope.json',{'generator':'exact spherical Tx/target/Rx distances and Tx waveform TTD','focusing':'spherical Rx delay vs plane-wave Rx delay; Tx remains the implemented far-field-steered aperture','angle_grid_step_deg':.25,'range_grid_step_m':.002,'carrier_space_phase':'counted once in generation; complex-baseband delay retained','no_strength_tuning':True})

def offgrid(out):
    rows=[]
    for step in [1.,.25,.1,.05]:
        aa=np.arange(-50,50+step/2,step)
        for true_angle in [-41.37,-17.23,0.37,24.13,47.41]:
            cfg=Config();true_range=50.5037;rr=np.arange(50.46,50.55,.001)
            for method in METHODS:
                raw,mf,_=point_raw(cfg,method,[true_range],[true_angle],geometry='farfield')
                im=reconstruct_points(cfg,method,raw,mf,rr,aa)
                ii=np.unravel_index(np.abs(im).argmax(),im.shape)
                rows.append({'method':method,'angular_step_deg':step,'truth_angle_deg':true_angle,'truth_range_m':true_range,'estimate_angle_deg':float(aa[ii[1]]),'estimate_range_m':float(rr[ii[0]]),'absolute_angle_error_deg':float(abs(aa[ii[1]]-true_angle)),'absolute_range_error_mm':float(abs(rr[ii[0]]-true_range)*1000),'peak':float(np.abs(im[ii]))})
    write_csv(out/'metrics.csv',rows)
    grouped=[]
    for method in METHODS:
        for step in [1.,.25,.1,.05]:
            q=[x for x in rows if x['method']==method and x['angular_step_deg']==step]
            grouped.append({'method':method,'angular_step_deg':step,'max_angle_error_deg':max(x['absolute_angle_error_deg'] for x in q)})
    curves(grouped,'angular_step_deg','max_angle_error_deg',out/'grid_convergence','Reconstruction angle spacing (deg)','Maximum angular error (deg)')

def boundaries(out):
    cfg=Config();rows=[];border=cfg.edges
    truth=np.unique(np.concatenate([border,np.clip(border-.13,-50,50),np.clip(border+.13,-50,50)]))
    aa=np.arange(-50,50.0001,.1);rr=np.arange(50.45,50.55,.002)
    for a in truth:
        for method in METHODS:
            raw,mf,_=point_raw(cfg,method,[50.5],[a],geometry='farfield')
            im=reconstruct_points(cfg,method,raw,mf,rr,aa);profile=np.max(np.abs(im),axis=0)
            ii=np.unravel_index(np.abs(im).argmax(),im.shape)
            candidates,_=signal.find_peaks(profile,height=.1*profile.max(),distance=max(1,int(1/.1)))
            candidates=np.unique(np.r_[candidates,0 if profile[0]>=profile[1] else [],len(profile)-1 if profile[-1]>=profile[-2] else []]).astype(int)
            candidates=candidates[profile[candidates]>=.1*profile.max()]
            peaks=aa[candidates]
            rows.append({'method':method,'truth_angle_deg':float(a),'estimate_angle_deg':float(aa[ii[1]]),'angle_error_deg':float(aa[ii[1]]-a),'peak_amplitude':float(profile.max()),'peaks_above_minus20db':int(len(peaks)),'matched_peaks_within_1deg':int(np.sum(abs(peaks-a)<=1)),'false_peaks_beyond_1deg':int(np.sum(abs(peaks-a)>1)),'detector_note':'diagnostic relative-threshold peaks; calibrated false alarm in R1-23'})
    write_csv(out/'metrics.csv',rows);curves(rows,'truth_angle_deg','peak_amplitude',out/'boundary_gain','Target angle (deg)','Reconstructed peak amplitude')

def zc_parameters(out):
    rows=[]
    sets=[('implemented',Config().roots),('consecutive',(1,2,3,4,5,6)),('separated',(17,23,31,37,41,43))]
    for l in [251,499,751]:
        cfg=replace(Config(),code_length=l,duration=l/100000.)
        for name,roots in sets:
            w,meta=waveforms(cfg,'cdm',roots=roots)
            peak=0.;aper_auto=0.;periodic_residual=0.
            for i,j in itertools.combinations(range(6),2):peak=max(peak,float(np.abs(signal.correlate(w[:,i],w[:,j],method='fft')).max()))
            chip=np.arange(l);code=np.exp(-1j*np.pi*roots[0]*chip*(chip+1)/l)
            cyc=np.fft.ifft(np.abs(np.fft.fft(code))**2)/l;periodic_residual=float(np.max(abs(cyc[1:])))
            auto=signal.correlate(code,code)/l;lags=signal.correlation_lags(l,l);aper_auto=float(np.max(abs(auto[lags!=0])))
            rows.append({'code_length':l,'root_set':name,'roots':';'.join(map(str,roots)),'chip_rate_hz':cfg.chip_rate,'duration_ms':cfg.duration*1000,'max_aperiodic_cross_amplitude_db':20*np.log10(peak),'chip_aperiodic_autosidelobe_db':20*np.log10(max(aper_auto,1e-30)),'periodic_autosidelobe_max_amplitude':periodic_residual,'fixed_parameters':'chip rate/fs; duration changes with L, no truncation'})
    write_csv(out/'metrics.csv',rows);curves(rows,'code_length','max_aperiodic_cross_amplitude_db',out/'root_length_sensitivity','ZC length (chips)','Worst finite-pulse cross-correlation (dB)',group='root_set')

def beam_audit(out):
    cfg=Config();rows=[];coverage=[]
    for f in np.arange(400000,500001,5000):
        aa=np.unique(np.r_[np.arange(-50,50.0001,.05),cfg.edges]);tx=abs(tx_pattern(cfg,[f],aa)[0]);tx/=np.max(tx,axis=1,keepdims=True)
        for k in range(6):
            ix=(aa>=cfg.edges[k])&(aa<=cfg.edges[k+1]);coverage.append({'frequency_hz':f,'sector':k+1,'minimum_coverage_db':float(20*np.log10(tx[k,ix].min())),'worst_angle_deg':float(aa[ix][np.argmin(tx[k,ix])])})
        for scan in np.arange(-50,50.0001,1):
            k=int(np.clip(np.digitize(scan,cfg.edges[1:-1]),0,5));at=np.arange(max(-50,scan-8),min(50,scan+8)+.005,.01)
            tr=abs(tx_pattern(cfg,[f],at)[0,k]);u0=np.sin(np.deg2rad(scan));ut=np.sin(np.deg2rad(at))
            for method in ['dense','cdm']:
                n=cfg.n_dense if method=='dense' else cfg.n_sparse;d=(.5 if method=='dense' else 2)*cfg.c/cfg.fc
                rx=abs(array_factor(2*np.pi*f*d/cfg.c*(ut-u0),n))*np.sinc(cfg.element_width_lambda/cfg.fc*f*ut)
                two=rx if method=='dense' else rx*tr
                width=peak_width(at,two**2)
                alias=[];inside=[]
                if method!='dense':
                    desired=float(abs(tx_pattern(cfg,[f],[scan])[0,k,0]));
                    for q in range(-4,5):
                        if not q:continue
                        ug=u0+q*cfg.c/(f*d)
                        if abs(ug)>1:continue
                        ag=float(np.rad2deg(np.arcsin(ug)));v=float(abs(tx_pattern(cfg,[f],[ag])[0,k,0]))/max(desired,1e-30)
                        ratio=20*np.log10(max(v,1e-30));alias.append(ratio)
                        if -50<=ag<=50:inside.append(ratio)
                rows.append({'method':method,'frequency_hz':f,'scan_angle_deg':scan,'assigned_sector':k+1,'local_power_hpbw_deg':width,'local_peak_angle_shift_deg':float(at[np.argmax(two)]-scan),'tx_only_alias_ratio_in_fov_db':max(inside) if inside else None,'tx_only_alias_ratio_all_physical_visible_db':max(alias) if alias else None,'hpbw_boundary_censored':width is None})
    write_csv(out/'metrics.csv',rows);write_csv(out/'sector_coverage.csv',coverage)
    q=[x for x in rows if x['frequency_hz']==450000];curves(q,'scan_angle_deg','local_power_hpbw_deg',out/'full_fov_hpbw','Receive steering angle (deg)','Two-way local power HPBW (deg)')
    summary=[]
    for method in ['dense','cdm']:
        q=[x for x in rows if x['method']==method];valid=[x for x in q if x['local_power_hpbw_deg'] is not None];worst=max(valid,key=lambda x:x['local_power_hpbw_deg'])
        summary.append({'method':method,'worst_uncensored_hpbw_deg':worst['local_power_hpbw_deg'],'at_angle_deg':worst['scan_angle_deg'],'at_frequency_hz':worst['frequency_hz'],'censored_boundary_rows':sum(x['hpbw_boundary_censored'] for x in q),'worst_alias_in_fov_db':max((x['tx_only_alias_ratio_in_fov_db'] for x in q if x['tx_only_alias_ratio_in_fov_db'] is not None),default=None)})
    write_csv(out/'extrema.csv',summary)
    dump(out/'scope.json',{'primary_domain_deg':[-50,50],'frequency_step_hz':5000,'fine_pattern_step_deg':.01,'coverage_step_deg':.05,'all_boundaries_inserted':True,'grating_definition':'Tx amplitude at receive alias / Tx amplitude at intended incident angle; in-FOV primary, all physically visible orders separately recorded','HPBW_at_endpoints':'censored when half-power crossing lies outside evaluation FOV, not silently extrapolated'})

def _locate_point(cfg,method,raw,mf,rx_offset=0.):
    rr=np.arange(49.25,51.7501,.01);aa=np.arange(-50,50.0001,1.)
    im=reconstruct_points(cfg,method,raw,mf,rr,aa,rx_steering_deg=rx_offset);ii=np.unravel_index(abs(im).argmax(),im.shape)
    refined_r=np.arange(rr[ii[0]]-.012,rr[ii[0]]+.0121,.001)
    refined_a=np.arange(max(-50,aa[ii[1]]-1),min(50,aa[ii[1]]+1)+.001,.1)
    fine=reconstruct_points(cfg,method,raw,mf,refined_r,refined_a,rx_steering_deg=rx_offset);jj=np.unravel_index(abs(fine).argmax(),fine.shape)
    outside=abs(aa-aa[ii[1]])>8;ghost=np.max(abs(im[:,outside]))/max(float(abs(fine[jj])),1e-30)
    return float(refined_r[jj[0]]),float(refined_a[jj[1]]),float(abs(fine[jj])),float(20*np.log10(max(ghost,1e-30)))

def calibration(out):
    definitions={'sync_sd_us':[0,.1,.25,.5,1.], 'sound_speed_fraction':[0,.001,.005,.01,.02], 'tx_steering_deg':[0,.1,.25,.5,1.], 'rx_steering_deg':[0,.1,.25,.5,1.], 'position_lambda':[0,.005,.01,.025,.05], 'gain_sd_db':[0,.1,.5,1.,2.], 'phase_sd_deg':[0,1,3,5,10]}
    random={'sync_sd_us','position_lambda','gain_sd_db','phase_sd_deg'};rows=[]
    for factor,values in definitions.items():
        for rep in range(3 if factor in random else 1):
            cfg=replace(Config(),seed=2026100950+rep,interpolation_oversampling=8)
            for method in METHODS:
                baseline=None
                for value in values:
                    err={} if factor=='rx_steering_deg' else {factor:value}
                    raw,mf,_=point_raw(cfg,method,[50.5],[20.],geometry='farfield',errors=err)
                    r,a,pk,gh=_locate_point(cfg,method,raw,mf,value if factor=='rx_steering_deg' else 0)
                    if value==0:baseline=pk
                    rows.append({'factor':factor,'factor_value':value,'realization':rep,'method':method,'angle_error_deg':a-20.,'range_error_mm':(r-50.5)*1000,'peak_loss_db':20*np.log10(pk/baseline),'outside_8deg_peak_db':gh,'true_generator':'perturbed forward; nominal reconstruction','error_scope':'synthetic pressure-test range, not measured tolerances'})
            print('calibration',factor,'realization',rep,flush=True)
    write_csv(out/'metrics.csv',rows)
    for factor in definitions:
        curves([r for r in rows if r['factor']==factor and r['realization']==0],'factor_value','peak_loss_db',out/('sensitivity_'+factor),factor,'Peak change (dB)')
    dump(out/'error_definitions.json',{'ranges':definitions,'random_errors':'zero-mean independent Gaussian per element; deterministic common-prefix standardized draws; 3 seeds','sound_speed':'positive uniform speed mismatch; fixed gate uses nominal c','rx_steering':'actual scan weights offset; image labels remain nominal','localization':'blind global coarse search followed by local refinement around detected maximum, not truth','range_step_refined_m':.001,'angle_step_refined_deg':.1})

def nonideal(out):
    cfg=Config();rows=[];f=cfg.fc+np.fft.fftfreq(8192,1/cfg.fs)
    for method in METHODS:
        clean,mf,_=point_raw(cfg,method,[50.5],[20.],geometry='farfield');baseline=_locate_point(cfg,method,clean,mf)[2]
        xn=positions(cfg,method)
        for dr in [.03,.15,.6]:
            for level in [-30,-20,-10,0]:
                utrue=np.sin(np.deg2rad(20));upath=np.sin(np.deg2rad(-20))
                extra=10**(level/20)*np.exp(1j*np.pi/3)*np.exp(-2j*np.pi*f[:,None]*(2*dr/cfg.c-xn[None,:]*(upath-utrue)/cfg.c))
                ep_ratio=np.sinc(cfg.element_width_lambda/cfg.fc*f*upath)/np.sinc(cfg.element_width_lambda/cfg.fc*f*utrue)
                raw=clean*(1+extra*ep_ratio[:,None]);r,a,pk,gh=_locate_point(cfg,method,raw,mf)
                rows.append({'method':method,'factor':'additional_path','value':level,'extra_equivalent_range_m':dr,'range_error_mm':(r-50.5)*1000,'angle_error_deg':a-20,'peak_change_db':20*np.log10(pk/baseline),'outside_8deg_peak_db':gh})
        for factor,values in {'transducer_ripple_db':[0,1,2,3,6],'transducer_phase_deg':[0,5,10,20,40],'element_width_fraction':[0,.1,.25,.5]}.items():
            for v in values:
                raw,mf,_=point_raw(cfg,method,[50.5],[20.],geometry='farfield',errors={factor:v});r,a,pk,gh=_locate_point(cfg,method,raw,mf)
                rows.append({'method':method,'factor':factor,'value':v,'extra_equivalent_range_m':0,'range_error_mm':(r-50.5)*1000,'angle_error_deg':a-20,'peak_change_db':20*np.log10(pk/baseline),'outside_8deg_peak_db':gh})
    write_csv(out/'metrics.csv',rows)
    for factor in ['transducer_ripple_db','transducer_phase_deg','element_width_fraction']:
        curves([r for r in rows if r['factor']==factor],'value','peak_change_db',out/factor,factor,'Peak change (dB)')
    dump(out/'scope.json',{'multipath':'same transmitted waveform mixture from the physical target; second Rx arrival has prescribed delay and DOA -20deg, phase pi/3; not a second newly illuminated physical object','transducer':'synthetic sinusoidal passband magnitude/phase ripple or width perturbation; nominal matched filter','not_modelled':'full reverberation, water-tank propagation or measured transducers'})

def wideband(out):
    cfg=Config();rows=[];checks=[]
    for method in METHODS:
        w,_=waveforms(cfg,method);nf=8192;s=np.fft.fft(w,nf,axis=0)*frontend_response(cfg,nf)[:,None];freq=cfg.fc+np.fft.fftfreq(nf,1/cfg.fs)
        n=cfg.n_dense if method=='dense' else cfg.n_sparse;d=(.5 if method=='dense' else 2)*cfg.c/cfg.fc
        xn=positions(cfg,method);wr=weights(n);char=[]
        for k in range(s.shape[1]):
            ps=abs(s[:,k])**2;ps/=ps.sum()
            # Characteristic function of the FULL FFT-bin spectral energy.
            # Equivalent to integrating the explicit broadband element sum;
            # 128x Fourier interpolation avoids billions of repeated sines.
            char.append(signal.resample(np.fft.ifft(ps)*nf,nf*128))
        for target in np.arange(-50,50.0001,1):
            k=0 if method=='dense' else int(np.clip(np.digitize(target,cfg.edges[1:-1]),0,5))
            ps=abs(s[:,k])**2;ps/=ps.sum()
            # Frequency integration is independent of the channel-centre approximation.
            scan=np.arange(-50,50.0001,.1);delta=np.sin(np.deg2rad(target))-np.sin(np.deg2rad(scan))
            tau=xn[:,None]*delta[None,:]/cfg.c;idx=tau*cfg.fs*128;low=np.floor(idx).astype(int);fraction=idx-low;cs=char[k]
            val=cs[low%len(cs)]*(1-fraction)+cs[(low+1)%len(cs)]*fraction
            exact=np.sum(wr[:,None]*val*np.exp(2j*np.pi*cfg.fc*tau),axis=0)
            fk=cfg.fc if method!='fdm' else cfg.centres[k]
            central=array_factor(2*np.pi*d/cfg.c*fk*delta,n)
            if target in [0,40]:
                ix=np.array([0,130,350,500,600,850,1000]);direct=np.sum(ps[:,None]*array_factor(2*np.pi*d/cfg.c*freq[:,None]*delta[ix][None,:],n),axis=0)
                checks.append({'method':method,'target_deg':target,'max_full_bin_integration_difference':float(np.max(abs(direct-exact[ix]))),'fixed_tolerance':.0001})
            local=abs(scan-target)<=8
            ae=float(scan[local][np.argmax(abs(exact[local]))]);ac=float(scan[local][np.argmax(abs(central[local]))])
            false=~local
            rows.append({'method':method,'target_angle_deg':float(target),'max_local_angle_difference_deg':abs(ae-ac),'exact_outside_main_db':float(20*np.log10(max(np.max(abs(exact[false])),1e-30)/max(np.max(abs(exact[local])),1e-30))),'centre_outside_main_db':float(20*np.log10(max(np.max(abs(central[false])),1e-30)/max(np.max(abs(central[local])),1e-30))),'max_normalized_pattern_error':float(np.max(abs(exact/exact[local].max()-central/central[local].max()))),'frequency_model':'energy-weighted physical fc+nu; centre reference uses channel RF centre'})
    write_csv(out/'metrics.csv',rows);curves(rows,'target_angle_deg','exact_outside_main_db',out/'wideband_alias','Target angle (deg)','Receive-only broadband residual outside main lobe (dB)')
    extrema=[]
    for method in METHODS:
        q=[r for r in rows if r['method']==method]
        extrema.append({'method':method,'max_angle_difference_deg':max(r['max_local_angle_difference_deg'] for r in q),'worst_exact_residual_db':max(r['exact_outside_main_db'] for r in q),'worst_centre_residual_db':max(r['centre_outside_main_db'] for r in q),'max_pattern_error':max(r['max_normalized_pattern_error'] for r in q)})
    write_csv(out/'extrema.csv',extrema)
    write_csv(out/'full_bin_quadrature_checks.csv',checks)
    dump(out/'scope.json',{'component':'receive-only broadband approximation diagnostic, separated from two-way grating audit and complete waveform image','maximum_angle_error':'max over 101 continuous evaluation angles within ±50deg; 0.1deg reconstruction grid; report 0 as below this numerical resolution, not exact zero','waveform_implementation':'new uniform Tx frontend pulse; CDM IFFT+resample, no temporal taper; energy-weighted receive-only test, not a full noise reconstruction','integration':'all 8192 physical-frequency bins; characteristic-function form with 128x Fourier interpolation; separately checked against full-bin direct quadrature'})
    if any(r['max_full_bin_integration_difference']>r['fixed_tolerance'] for r in checks):raise RuntimeError('Broadband quadrature numerical check failed')

def range_tradeoff(out):
    cfg=Config();rows=[];allprofiles={};distances=[.005,.01,.02,.04,.08,.12]
    for shaping in [False,True]:
        for method in METHODS:
            w,_=waveforms(cfg,method,shaping=shaping);auto=signal.correlate(w[:,0],w[:,0],method='fft')
            auto=signal.resample(auto,len(auto)*8);lags=(np.arange(len(auto))/8-(len(w)-1))/cfg.fs*cfg.c/2
            for separation in distances:
                x=np.arange(-.2,.2001,.0005)
                def interp(t):return np.interp(t,lags,auto.real,left=0,right=0)+1j*np.interp(t,lags,auto.imag,left=0,right=0)
                y=interp(x-separation/2)+interp(x+separation/2)*np.exp(1j*.7)
                pp=abs(y)**2;peaks,_=signal.find_peaks(pp,height=.25*pp.max(),distance=4)
                chosen=x[peaks];matched=all(np.any(abs(chosen-t)<max(.002,separation/4)) for t in [-separation/2,separation/2])
                valley=float(pp[np.argmin(abs(x))]/pp.max())
                rows.append({'method':method,'shaping':shaping,'target_separation_mm':separation*1000,'two_targets_resolved_under_declared_rule':bool(matched and valley<.8),'central_valley_power_ratio':valley,'detected_peaks':len(peaks),'phase_difference_rad':.7})
                if separation in [.01,.08] and not shaping:allprofiles[f'{method}_{separation}']=np.column_stack([x,abs(y)])
    write_csv(out/'metrics.csv',rows);np.savez_compressed(out/'two_target_profiles.npz',**allprofiles)
    fig,axs=plt.subplots(1,2,figsize=(9,3.7))
    for ax,sep in zip(axs,[.01,.08]):
        for i,method in enumerate(METHODS):
            p=allprofiles[f'{method}_{sep}'];ax.plot(p[:,0]*1000,db(p[:,1]),['-','--','-.'][i],label=LABELS[method],color=COLORS[method])
        ax.set(xlabel='Range offset (mm)',ylabel='Amplitude (dB)',title=f'Separation {sep*1000:g} mm',ylim=(-40,3),xlim=(-120,120));ax.legend()
    savefig(fig,out/'range_resolution_tradeoff')
    dump(out/'scope.json',{'resolution_rule':'blind local peaks above -6dB power; each within max(2mm, separation/4) and central power valley <0.8; rule fixed before results','model':'two coherent point echoes with phase difference 0.7rad; not a universal Rayleigh criterion','fair_taper_spectrum_controls':'separately reported in R1-12'})

def doppler(out):
    cfg=Config();rows=[]
    for method in METHODS:
        w,_=waveforms(cfg,method);n=len(w)
        for v in [-10.,-5.,-2.,-1.,0.,1.,2.,5.,10.]:
            beta=(cfg.c-v)/(cfg.c+v);fd=(beta-1)*cfg.fc
            for k in range(w.shape[1]):
                for model in ['frequency_shift','time_scaling_and_shift']:
                    if model=='frequency_shift':x=w[:,k].copy()
                    else:x=signal.resample(w[:,k],int(round(n/beta)))
                    t=np.arange(len(x))/cfg.fs;x=x*np.exp(2j*np.pi*fd*t);x/=np.linalg.norm(x)
                    corr=signal.correlate(x,w[:,k],method='fft');lag=signal.correlation_lags(len(x),n);i=int(abs(corr).argmax())
                    rows.append({'method':method,'channel':k+1,'radial_velocity_mps':v,'model':model,'doppler_hz':fd,'time_scale':beta,'peak_loss_db':float(20*np.log10(max(abs(corr[i]),1e-30))),'range_bias_mm':float(lag[i]/cfg.fs*cfg.c/2*1000),'actual_scaled_samples':len(x)})
    write_csv(out/'metrics.csv',rows)
    summary=[]
    for method in METHODS:
        for v in [-10.,-5.,-2.,-1.,0.,1.,2.,5.,10.]:
            rr=[r for r in rows if r['method']==method and r['radial_velocity_mps']==v and r['model']=='time_scaling_and_shift']
            summary.append({'method':method,'radial_velocity_mps':v,'worst_peak_loss_db':min(r['peak_loss_db'] for r in rr),'worst_absolute_range_bias_mm':max(abs(r['range_bias_mm']) for r in rr)})
    curves(summary,'radial_velocity_mps','worst_peak_loss_db',out/'expanded_doppler','Relative radial velocity (m/s)','Worst-channel peak loss (dB)')
    dump(out/'scope.json',{'velocity_sign':'positive receding; received carrier shift is negative','beta':'(c-v)/(c+v), ideal monostatic constant radial velocity','tested_domain_mps':[-10,10],'application_claim':'pressure-test domain, not verified coverage of all USV/AUV scenarios','comparison':'frequency-only approximation vs resampled pulse time scaling plus carrier shift; amplitude normalized to isolate mismatch'})

def _save_images(out,cfg,scene,images,name='raw_images'):
    np.savez_compressed(out/(name+'.npz'),truth_real=scene.real.astype(np.float32),truth_imag=scene.imag.astype(np.float32),angles_deg=cfg.angles,ranges_m=cfg.ranges,**{m:abs(im).astype(np.float32) for m,im in images.items()})
    dump(out/(name+'_configuration.json'),asdict(cfg))

def truth_metrics(out):
    cfg=Config(angular_nodes=201);truth,cat=scene_grid(cfg,50,5,extent=.6,seed=2026100981)
    images={m:image_grid(cfg,truth,m)[0] for m in METHODS};rows=[]
    for m in METHODS:rows.append({'method':m,**metrics(truth,images[m],cfg,images['dense'])})
    # Validate metric conventions without an image solver or dense reference.
    checks=[]
    for label,test in [('identity',truth),('global_gain',3*truth),('zero_image',np.zeros_like(truth)),('off_support_energy',np.roll(truth,80,axis=0))]:
        q=metrics(truth,test,cfg);checks.append({'test':label,**q})
    write_csv(out/'metrics.csv',rows);write_csv(out/'metric_definition_checks.csv',checks);_save_images(out,cfg,truth,images);image_panel(cfg,images,out/'direct_truth_image_comparison')
    dump(out/'metric_definitions.json',{'input':'unprocessed complex scatter coefficient array; magnitude for direct truth metric; no gamma/smoothing','NMSE':'10log10(sum((a*abs(I)-abs(truth))^2)/sum(abs(truth)^2)); a>=0 minimizes linear-magnitude residual','correlation':'zero-mean Pearson correlation of linear magnitudes on fixed evaluated grid','leakage':'sum(abs(I)^2 outside fixed physical guard)/sum(abs(I)^2 over evaluated image), discrete sampled image energy, not acoustic radiation energy','CNR':'(mean magnitude on declared support - mean background magnitude)/std background magnitude','guard':'0.06m and 1deg, fixed before outcome','background_validity':'leakage/CNR unavailable when unguarded background fraction <5%; denominator retained in CSV','dense_relative':'separate reference-reconstruction similarity, not ground truth precision','grid':'uniform sin(theta); all methods use the same cells and evaluation domain'})

def density(out):
    cfg=Config(angular_nodes=201,range_span=2.4);rows=[];catalogue=[]
    cond=[('scatter_count',n,n,50,.8,0) for n in [1,10,50,200,500,1000]]
    cond += [('range_layers',l,1000,l,.8,0) for l in [1,5,10,20,50]]
    cond += [('delay_extent_m',e,500,50,e,0) for e in [.15,.3,.6,1.2,2.4]]
    cond += [('dynamic_range_db',d,500,50,.8,d) for d in [0,10,20,40]]
    for i,(factor,value,n,layers,extent,dynamic) in enumerate(cond):
        for rep in range(3):
            seed=2026101000+100*i+rep;scene,cat=scene_grid(cfg,n,layers,extent,dynamic,seed)
            images={m:image_grid(cfg,scene,m)[0] for m in METHODS}
            for m in METHODS:rows.append({'factor':factor,'factor_value':value,'scatter_count':n,'requested_layers':layers,'extent_m':extent,'dynamic_range_db':dynamic,'scene_realization':rep,'scene_seed':seed,'method':m,**metrics(scene,images[m],cfg,images['dense'])})
            catalogue.append({'factor':factor,'value':value,'scene_seed':seed,'individual_reflector_power':float(np.sum(abs(cat['coefficients'])**2)),'grid_coherent_power':float(np.sum(abs(scene)**2)),'occupied_pixels':int(np.count_nonzero(scene)),'occupied_range_cells':int(np.count_nonzero(np.any(scene!=0,axis=1)))})
            if factor=='scatter_count' and n==1000 and rep==0:
                _save_images(out,cfg,scene,images,'1000_scatter_50_layers');image_panel(cfg,images,out/'1000_scatter_50_layers')
        print('density',factor,value,flush=True)
    write_csv(out/'metrics.csv',rows);write_csv(out/'scene_catalogue.csv',catalogue)
    from scipy.stats import t as student_t
    aggregate=[]
    for factor in dict.fromkeys(r['factor'] for r in rows):
        for value in dict.fromkeys(r['factor_value'] for r in rows if r['factor']==factor):
            for method in METHODS:
                q=[r for r in rows if r['factor']==factor and r['factor_value']==value and r['method']==method]
                rr={'factor':factor,'factor_value':value,'method':method,'independent_scene_count':len(q),'CI_scope':'three independent generated scenes at fixed factor level'}
                for metric in ['truth_magnitude_nmse_db','dense_magnitude_correlation','support_leakage']:
                    a=np.array([r[metric] for r in q if r[metric] is not None]);rr[metric]=float(a.mean()) if len(a) else None
                    rr[metric+'_ci95_halfwidth']=float(student_t.ppf(.975,len(a)-1)*a.std(ddof=1)/np.sqrt(len(a))) if len(a)>1 else None
                aggregate.append(rr)
        for metric in ['dense_magnitude_correlation','support_leakage','truth_magnitude_nmse_db']:
            curves([r for r in aggregate if r['factor']==factor],'factor_value',metric,out/(factor+'_'+metric),factor,metric)
    write_csv(out/'aggregate_metrics.csv',aggregate)
    dump(out/'scope.json',{'generator':'common broadband far-field angular-grid channel model; not old reduced_model CSV','receive_counts':{'dense':128,'sparse':32},'fixed_physical_evaluation_range_m':[50,52.4],'angular_nodes':201,'maximum_angle_spacing_deg':float(np.max(np.diff(cfg.angles))),'reflector_energy':'sum individual reflectivity powers fixed at 1; grid collisions add coherently and actual grid power is recorded','paired_methods':True,'replications':'3 independent generated geometries per factor level; not 30 phase replicates','density_axes':'scatter count, actual occupied range cells, and delay extent kept distinct','test_conclusion':'not predetermined; absent valid metric background is retained as unavailable'})

def _ship_scene(seed=2026102001,nodes=201):
    path=PROJECT/'paperA_revise'/'11_A11_shipwreck_common_scene'/'results'/'A11_raytrace_arrivals.csv'
    data=np.loadtxt(path,delimiter=',',skiprows=1,usecols=(2,3,6))
    r,a,amp=data.T;valid=(r>=42)&(r<=122)&(a>=-50)&(a<=50)&(amp>0)
    r,a,amp=r[valid],a[valid],amp[valid]
    cfg=Config(range_origin=42,range_span=80,angular_nodes=nodes)
    phase=np.random.default_rng(seed).uniform(-np.pi,np.pi,len(r));coeff=amp*np.exp(1j*phase)/np.sqrt(np.sum(amp**2))
    scene=grid_catalog(cfg,r,a,coeff)
    return cfg,scene,{'source_csv':str(path),'source_sha256':sha(path),'valid_arrivals':len(r),'phase_seed':seed,'normalization':'individual geometry-derived coefficient powers sum to 1','geometry':'unchanged legacy transformed ship/seabed arrival cloud; existing angle/range columns used directly','unprocessed_amplitude':'MetricAmplitude from source; no TVG/gamma/smoothing','source_min_range_m':float(r.min()),'source_max_range_m':float(r.max())}

def scene_ordering(out):
    rows=[];provenance=[]
    for kind in ['synthetic1000x50','legacy_ship']:
        for shape in [False,True]:
            for rep in range(3):
                if kind=='legacy_ship':cfg,truth,meta=_ship_scene(2026102001+rep)
                else:
                    cfg=Config(range_span=2.4,angular_nodes=201);truth,cat=scene_grid(cfg,1000,50,extent=.8,seed=2026101901+rep)
                    meta={'phase_and_scene_seed':2026101901+rep,'normalization':cat['amplitude_norm']}
                images={m:image_grid(cfg,truth,m,shaping=shape)[0] for m in METHODS}
                for m in METHODS:rows.append({'scene':kind,'common_FIR_control':shape,'replication':rep,'method':m,**metrics(truth,images[m],cfg,images['dense'])})
                if rep==0:
                    tag=kind+('_FIR' if shape else '_original');_save_images(out,cfg,truth,images,tag);image_panel(cfg,images,out/tag)
                provenance.append({'scene':kind,'FIR':shape,'replication':rep,**meta})
                print('scene_ordering',kind,'FIR',shape,'rep',rep,flush=True)
    write_csv(out/'metrics.csv',rows);dump(out/'scene_provenance.json',provenance)
    dump(out/'scope.json',{'ranking':'all rankings retained; no assumption that new independent model reproduces old narrowband-model values','factor_controls':'density/delay/amplitude controls in R2-11; same passband FIR applied to every method here','phase_repetitions':'ship: 3 phases of one unchanged geometry; synthetic: 3 generated scenes, reported separately','STL_interpretation':'geometry-derived extended-scatter stress test, not new acoustic ray-tracing validation'})

def scene_diversity(out):
    cfg=Config(range_span=2.4,angular_nodes=201);rows=[]
    for kind in ['random','wall','two_clusters']:
        for rep in range(3):
            truth,cat=scene_grid(cfg,500,50,extent=.8,seed=2026103000+rep,shape=kind)
            images={m:image_grid(cfg,truth,m)[0] for m in METHODS}
            for m in METHODS:rows.append({'independent_structure':kind,'scene_realization':rep,'method':m,**metrics(truth,images[m],cfg,images['dense'])})
            if rep==0:_save_images(out,cfg,truth,images,kind);image_panel(cfg,images,out/kind)
    write_csv(out/'metrics.csv',rows)
    dump(out/'scope.json',{'independent_structures':['random spatial field','procedural gently varying wall','two separate clusters'],'same_count':500,'same_total_individual_reflectivity_power':1,'geometry_or_density_alternative':'both density configurations and these limited independent structures are provided; not claimed as multiple real shipwrecks','further_STL_geometry':'not necessary to rerun large CAD/MEX chain for this OR route'})

def common_display(out):
    from executors import run_bundle
    base=run_bundle('truth_metrics');saved=np.load(base/'raw_images.npz');cfg=Config(angular_nodes=201)
    images={m:saved[m] for m in METHODS};truth=saved['truth_real']+1j*saved['truth_imag']
    image_panel(cfg,images,out/'individual_peak_normalization',common=False);image_panel(cfg,images,out/'common_peak_reference',common=True)
    ref=max(float(abs(im).max()) for im in images.values());rows=[{'method':m,'own_peak':float(abs(im).max()),'common_reference':ref,'peak_relative_to_common_db':float(20*np.log10(abs(im).max()/ref)),'display_floor_db':-50} for m,im in images.items()]
    write_csv(out/'metrics.csv',rows);dump(out/'input_reference.json',{'shared_input':str(base),'source_npz_sha256':sha(base/'raw_images.npz'),'reference':'one scalar: maximum peak across the paired reconstructions, not separate method references','no_gamma_for_metrics':True})

def assignment(out):
    cfg=Config(angular_nodes=201,range_span=1.2);w,_=waveforms(cfg,'fdm');pair=np.zeros((6,6))
    for i,j in itertools.combinations(range(6),2):pair[i,j]=pair[j,i]=float(np.max(abs(signal.correlate(w[:,i],w[:,j],method='fft')))**2)
    overlap=np.zeros((6,6,6))
    for i in range(6):
        a=np.linspace(cfg.edges[i],cfg.edges[i+1],25)
        for j in range(6):
            for slot in range(6):overlap[i,j,slot]=np.mean(abs(tx_pattern(cfg,[cfg.centres[slot]],a)[0,j])**2)
    candidates=[]
    for perm in itertools.permutations(range(6)):
        # Explicit design proxy only: off-sector illumination times full-delay
        # waveform cross-energy peak. This is NOT global image-error optimum.
        score=sum(overlap[i,j,perm[j]]*pair[perm[i],perm[j]] for i in range(6) for j in range(6) if i!=j)
        candidates.append({'permutation':';'.join(map(str,perm)),'proxy_score':float(score)})
    candidates.sort(key=lambda x:x['proxy_score']);write_csv(out/'all_720_proxy_candidates.csv',candidates)
    orders={'implemented':tuple(range(6)),'ascending_frequency':tuple(np.argsort(cfg.centres)),'descending_frequency':tuple(np.argsort(cfg.centres)[::-1]),'proxy_best':tuple(map(int,candidates[0]['permutation'].split(';'))),'proxy_worst':tuple(map(int,candidates[-1]['permutation'].split(';')))}
    truth,_=scene_grid(cfg,200,20,extent=.8,seed=2026104101);dense=image_grid(cfg,truth,'dense')[0];rows=[]
    for name,perm in orders.items():
        test=replace(cfg,centres=tuple(cfg.centres[i] for i in perm));im,_=image_grid(test,truth,'fdm');rows.append({'candidate':name,'sector_frequency_centres_khz':';'.join(str(test.centres[i]/1000) for i in range(6)),**metrics(truth,im,cfg,dense)})
    write_csv(out/'metrics.csv',rows);dump(out/'selection_scope.json',{'original_assignment':'implemented permutation retained; no original optimization provenance invented','new_proxy':'sum off-sector illumination weight * squared full-delay cross-correlation peak; computed for all 720 permutations','representative_image_candidates':list(orders),'caveat':'proxy overlap weights use nominal geometry/frequency illumination; proxy rank is not proof of best full imaging design'})

def sector_search(out):
    cfg=Config();freq=np.arange(400000,500001,5000);cache={};candidate_rows=[]
    def scores(candidates,edges,step=.25,frequencies=freq):
        targets=np.unique(np.r_[np.arange(edges[0],edges[1]+step/2,step),edges])
        cov=np.full(len(candidates),np.inf);gl=np.full(len(candidates),-np.inf)
        for f in frequencies:
            response=abs(tx_pattern(cfg,[f],targets,centres=candidates)[0])
            # Global mainbeam is audited on a dense FOV grid at final stage.
            # Search normalization samples the local mainbeam within ±1deg;
            # Hamming sidelobes cannot dominate this mainbeam in this FOV.
            near=candidates[:,None]+np.arange(-1,1.001,.01)[None,:];near=np.clip(near,-50,50)
            q=2*np.pi*f/(2*cfg.fc)*(np.sin(np.deg2rad(near))-np.sin(np.deg2rad(candidates[:,None])))
            normpeak=np.max(abs(array_factor(q,cfg.n_tx)*np.sinc(cfg.element_width_lambda*f/cfg.fc*np.sin(np.deg2rad(near)))),axis=1)
            cov=np.minimum(cov,np.min(response/normpeak[:,None],axis=1))
            u=np.sin(np.deg2rad(targets))
            for order in range(-4,5):
                if not order:continue
                ug=u+order*cfg.fc/(2*f);valid=(abs(ug)<=1)&(ug>=np.sin(np.deg2rad(-50)))&(ug<=np.sin(np.deg2rad(50)))
                if not np.any(valid):continue
                alias=np.rad2deg(np.arcsin(ug[valid]));ar=abs(tx_pattern(cfg,[f],alias,centres=candidates)[0]);ratio=ar/np.maximum(response[:,valid],1e-30)
                gl=np.maximum(gl,np.max(ratio,axis=1))
        return 20*np.log10(cov),20*np.log10(np.maximum(gl,1e-30))
    for n in range(4,9):
        edges=np.linspace(-50,50,n+1);sectors=[]
        for k in range(n):
            candidates=np.unique(np.r_[np.arange(edges[k],edges[k+1]+.025,.05),edges[k+1]])
            co,gr=scores(candidates,edges[k:k+2]);sectors.append((candidates,co,gr))
            candidate_rows.extend({'sector_count':n,'sector':k+1,'steering_deg':float(a),'search_min_coverage_db':float(c),'search_worst_in_fov_grating_db':float(g)} for a,c,g in zip(candidates,co,gr))
        cache[n]=(edges,sectors);print('joint sector candidates',n,flush=True)
    write_csv(out/'all_joint_candidates.csv',candidate_rows)
    threshold_rows=[];selected=[]
    for cov_limit,gl_limit in itertools.product([-6.,-6.3,-6.6,-7.,-8.],[-18.,-21.,-24.,-27.,-30.]):
        feasible=None;centres=[];margins=[]
        for n in range(4,9):
            edges,sectors=cache[n];choices=[];ok=True
            for candidates,co,gr in sectors:
                good=(co>=cov_limit)&(gr<=gl_limit)
                if not good.any():ok=False;break
                indices=np.where(good)[0];pick=indices[np.argmax(co[indices])];choices.append(float(candidates[pick]))
            if ok:feasible=n;centres=choices;break
        threshold_rows.append({'coverage_threshold_db':cov_limit,'grating_threshold_db':gl_limit,'first_search_grid_feasible_count':feasible,'steering_angles_deg':';'.join(map(str,centres)) if centres else '', 'search_scope':'equal-width sectors, all 0.05deg steering candidates, 0.25deg targets + edges, 5kHz RF grid','not_global_continuous_optimum':True})
        if cov_limit==-6.6 and gl_limit==-24 and feasible is not None:selected=centres
    write_csv(out/'metrics.csv',threshold_rows)
    # Independently audit actual implemented and selected nominal-threshold
    # candidate on the requested fine angle/frequency grid, not search scores.
    audits=[]
    for name,centres in [('implemented',list(cfg.focus)),('joint_selected',selected)]:
        if not centres:continue
        n=len(centres);edges=np.linspace(-50,50,n+1);full=np.unique(np.r_[np.arange(-50,50.0001,.01),edges]);worst_cov=1.;worst_gl=0.;atc=None;atg=None
        for f in np.arange(400000,500001,1000):
            alltx=abs(tx_pattern(cfg,[f],full,centres=centres)[0])
            for k in range(n):
                valid=(full>=edges[k])&(full<=edges[k+1]);target=full[valid];response=alltx[k,valid];coverage=response/alltx[k].max();ii=int(coverage.argmin())
                if coverage[ii]<worst_cov:worst_cov=float(coverage[ii]);atc=(float(target[ii]),int(f))
                u=np.sin(np.deg2rad(target))
                for order in range(-4,5):
                    if not order:continue
                    ug=u+order*cfg.fc/(2*f);visible=(abs(ug)<=1)&(ug>=np.sin(np.deg2rad(-50)))&(ug<=np.sin(np.deg2rad(50)))
                    if not visible.any():continue
                    alias=np.rad2deg(np.arcsin(ug[visible]));ar=abs(tx_pattern(cfg,[f],alias,centres=[centres[k]])[0,0]);ratio=ar/np.maximum(response[visible],1e-30);ii=int(ratio.argmax())
                    if ratio[ii]>worst_gl:worst_gl=float(ratio[ii]);atg=(float(target[visible][ii]),float(alias[ii]),int(f),order)
        cd=20*np.log10(worst_cov);gd=20*np.log10(worst_gl)
        audits.append({'design':name,'sector_count':n,'steering_angles_deg':';'.join(map(str,centres)),'dense_worst_coverage_db':float(cd),'coverage_angle_deg':atc[0],'coverage_frequency_hz':atc[1],'dense_worst_in_fov_grating_db':float(gd),'grating_target_deg':atg[0],'grating_incident_angle_deg':atg[1],'grating_frequency_hz':atg[2],'grating_order':atg[3],'dense_nominal_constraints_satisfied':bool(cd>=-6.6 and gd<=-24)})
    write_csv(out/'dense_final_audit.csv',audits)
    fig,ax=plt.subplots(figsize=(5,3.5));matrix=np.array([[next(r['first_search_grid_feasible_count'] or np.nan for r in threshold_rows if r['coverage_threshold_db']==c and r['grating_threshold_db']==g) for g in [-18,-21,-24,-27,-30]] for c in [-6.,-6.3,-6.6,-7.,-8.]])
    h=ax.imshow(matrix,vmin=4,vmax=8,cmap='viridis',aspect='auto');ax.set_xticks(range(5),[-18,-21,-24,-27,-30]);ax.set_yticks(range(5),[-6.,-6.3,-6.6,-7.,-8.]);ax.set(xlabel='Grating constraint (dB)',ylabel='Coverage constraint (dB)')
    for i in range(5):
        for j in range(5):ax.text(j,i,'none' if not np.isfinite(matrix[i,j]) else str(int(matrix[i,j])),ha='center',va='center',color='white' if not np.isfinite(matrix[i,j]) or matrix[i,j]<7 else 'black')
    fig.colorbar(h,ax=ax,label='First feasible candidate count');savefig(fig,out/'threshold_selection_map')
    dump(out/'search_proof_scope.json',{'joint_test':True,'candidate_counts':[4,5,6,7,8],'steering_grid_deg':.05,'search_target_grid_deg':.25,'search_frequency_grid_hz':freq.tolist(),'final_target_grid_deg':.01,'final_frequency_spacing_hz':1000,'factorization':'for the fixed equal-width partition, per-sector Tx constraints depend on its own steering candidate, so exhaustive per-sector feasibility is exhaustive within this finite design family','threshold_reuse':'reuse every candidate constraint score, then search all candidates under each threshold, not just one coverage-best angle','normalization_caveat':'search peak sampled locally at 0.01deg around steering; final uses full dense FOV maximum; final feasibility is reported separately, not assumed','continuous_minimum':'not proven; finer steering candidates or nonuniform boundaries not covered'})

def _blind_detection(image,rr,aa,threshold,truth_r,truth_a):
    power=abs(image);sr=max(1,int(np.ceil(.03/(rr[1]-rr[0]))));sa=max(1,int(np.ceil(2/(aa[1]-aa[0]))))
    local=ndimage.maximum_filter(power,size=(2*sr+1,2*sa+1),mode='constant',cval=-np.inf)
    ix=np.argwhere((power==local)&(power>=threshold));ix=ix[np.argsort(power[ix[:,0],ix[:,1]])[::-1]] if len(ix) else ix
    assigned=set();tp=0;false=0
    for i,j in ix:
        distance=((truth_r-rr[i])/.02)**2+((truth_a-aa[j])/1.)**2;eligible=[k for k in np.argsort(distance) if k not in assigned and abs(truth_r[k]-rr[i])<=.02 and abs(truth_a[k]-aa[j])<=1.]
        if eligible:assigned.add(int(eligible[0]));tp+=1
        else:false+=1
    return tp/len(truth_r),false,len(ix)

def common_noise(out):
    cfg=Config();rr=np.arange(50,50.8001,.01);aa=np.arange(-50,50.0001,.5);nf=8192
    truth_r=np.repeat([50.2,50.4,50.6],6);truth_a=np.tile([-40.,-24.,-8.,8.,24.,40.],3)
    coeff=np.exp(1j*np.random.default_rng(2026105001).uniform(-np.pi,np.pi,18))/np.sqrt(18)
    clean={};raws={};filters={};metadata={};signal_power={}
    for m in METHODS:
        raw,mf,meta=point_raw(cfg,m,truth_r,truth_a,coeff,geometry='farfield');raws[m]=raw;filters[m]=mf;metadata[m]=meta
        clean[m]=reconstruct_points(cfg,m,raw,mf,rr,aa);signal_power[m]=float(np.mean(abs(np.fft.ifft(raw,axis=0))**2))
    cal=[];testnoise={m:[] for m in METHODS};calibration_count=60;test_count=30
    for rep in range(calibration_count+test_count):
        rng=np.random.default_rng(2026105100+rep);z=(rng.standard_normal((nf,128))+1j*rng.standard_normal((nf,128)))/np.sqrt(2);zf=np.fft.fft(z,axis=0)
        for m in METHODS:
            nr=128 if m=='dense' else 32
            ni=reconstruct_points(cfg,m,zf[:,:nr],filters[m],rr,aa)
            if rep<calibration_count:cal.append({'method':m,'calibration_repetition':rep,'pure_noise_max_amplitude':float(abs(ni).max())})
            else:testnoise[m].append(ni.astype(np.complex64))
        if rep%15==0:print('noise pure-noise/calibration blocks',rep+1,'/',calibration_count+test_count,flush=True)
    thresholds={m:float(np.quantile([r['pure_noise_max_amplitude'] for r in cal if r['method']==m],.95,method='higher')) for m in METHODS}
    write_csv(out/'pure_noise_calibration.csv',cal)
    rows=[];null_rows=[]
    for energy in ['bank_unit_energy','tx_element_unit_energy']:
        for normal in ['equal_raw_channel_snr','common_per_element_noise']:
            for snr in [-35,-25,-15,-5]:
                for m in METHODS:
                    gain=1 if energy=='bank_unit_energy' else 1/np.sqrt(metadata[m]['unscaled_tx_element_excitation_energy'])
                    sig_power=signal_power[m]*gain**2
                    sigma=np.sqrt((sig_power if normal=='equal_raw_channel_snr' else signal_power['dense'])/10**(snr/10))
                    for rep,noise in enumerate(testnoise[m]):
                        im=clean[m]*gain+noise*sigma;threshold=thresholds[m]*sigma
                        pd,fp,npeaks=_blind_detection(im,rr,aa,threshold,truth_r,truth_a)
                        rows.append({'method':m,'source_energy_mode':energy,'noise_normalization':normal,'nominal_dense_or_equal_input_snr_db':snr,'noise_repetition':rep,'noise_sigma':float(sigma),'actual_raw_channel_snr_db':float(10*np.log10(sig_power/sigma**2)),'calibrated_amplitude_threshold':float(threshold),'detection_probability':float(pd),'false_detections':fp,'total_blind_detections':npeaks,'input_reflector_power':1.,'total_tx_element_excitation_energy':float(metadata[m]['unscaled_tx_element_excitation_energy'] if energy=='bank_unit_energy' else 1.)})
                        if energy=='bank_unit_energy' and normal=='common_per_element_noise' and snr==-15:null_rows.append({'method':m,'test_repetition':rep,'pure_noise_exceeded_calibrated_block_threshold':bool(abs(noise).max()>thresholds[m])})
    write_csv(out/'metrics.csv',rows);write_csv(out/'independent_false_alarm_test.csv',null_rows)
    aggregate=[]
    for energy,normal,snr,m in itertools.product(['bank_unit_energy','tx_element_unit_energy'],['equal_raw_channel_snr','common_per_element_noise'],[-35,-25,-15,-5],METHODS):
        q=[r for r in rows if r['source_energy_mode']==energy and r['noise_normalization']==normal and r['nominal_dense_or_equal_input_snr_db']==snr and r['method']==m];pd=np.array([r['detection_probability'] for r in q]);fp=np.array([r['false_detections'] for r in q])
        aggregate.append({'method':m,'source_energy_mode':energy,'noise_normalization':normal,'snr_db':snr,'mean_detection_probability':float(pd.mean()),'sd_detection_probability':float(pd.std(ddof=1)),'mean_false_detections':float(fp.mean()),'noise_trials':len(q),'CI_scope':'noise variation of one fixed 18-target scene'})
    write_csv(out/'aggregate_metrics.csv',aggregate)
    for normal in ['equal_raw_channel_snr','common_per_element_noise']:
        curves([r for r in aggregate if r['noise_normalization']==normal and r['source_energy_mode']=='tx_element_unit_energy'],'snr_db','mean_detection_probability',out/normal,'Nominal dense / equal input SNR (dB)','Detection probability')
    np.savez_compressed(out/'fixed_clean_images.npz',angles_deg=aa,ranges_m=rr,truth_ranges=truth_r,truth_angles=truth_a,**{m:abs(im).astype(np.float32) for m,im in clean.items()})
    dump(out/'scope.json',{'threshold':'95th percentile of per-block pure-noise maximum from 60 independent calibration blocks; not relative to test image peak','test_trials':30,'calibration_test_seeds_disjoint':True,'blind_search':'full image -> physical NMS -> one-to-one truth matching; truth not used to generate candidates','noise':'independent CN(0,sigma^2) per raw receiver sample/channel; common standardized prefix across methods; same 400kHz processed rate','two_SNR_questions':'mechanism at equal raw-channel SNR vs physical channel-noise sensitivity at common per-element variance','source_energy_controls':'original unit waveform-bank energy and supplementary equal sum-of-Tx-element excitation energy; exact energy includes simultaneous waveform cross terms','dense_Tx':'ideal broad point transmitter; sparse uses stated 12-element sectorized Tx; transmit gain is explicit, not hardware link-budget calibration','dimensions':{'fft_record_samples':8192,'range_reconstruction_spacing_m':.01,'angular_reconstruction_spacing_deg':.5},'range_resolution_accuracy':'fine off-grid localization is separately tested in R1-22; noise scene deliberately uses grid-aligned distinct targets','transducer_band':'flat numerical processed band; no measured frontend calibration'})

def competitor(out):
    """A real 32-channel nonuniform-array wideband OMP reconstruction.

    The dictionary is range-angle joint, the data use continuous coordinates
    and spherical receive paths. K is NOT supplied to OMP. No truth-centered
    covariance, no exact source-count MUSIC, and no point-only range oracle.
    """
    cfg=Config(interpolation_oversampling=8);nf=8192
    rng=np.random.default_rng(2026106001);span=np.ptp(positions(cfg,'fdm'))
    xpos=np.sort(np.r_[-span/2,rng.uniform(-span/2,span/2,30),span/2])
    rr=np.arange(50.3,50.7001,.01);aa=np.arange(-50,50.001,.5)
    nu=np.fft.fftfreq(nf,1/cfg.fs);possible=np.where(abs(nu)<=45000)[0]
    bins=np.sort(np.random.default_rng(2026106002).choice(possible,41,replace=False))
    sf=np.fft.fft(waveforms(cfg,'dense')[0][:,0],nf)[bins];f=cfg.fc+nu[bins]
    rgrid,agrid=np.meshgrid(rr,aa,indexing='ij');u=np.sin(np.deg2rad(agrid.ravel()))
    delay=(2*(rgrid.ravel()[None,:]-cfg.gate_range)-xpos[:,None]*u[None,:])/cfg.c
    ep=np.sinc(cfg.element_width_lambda/cfg.fc*f[:,None]*u[None,:])
    dictionary=(sf[:,None,None]*ep[:,None,:]*np.exp(-2j*np.pi*f[:,None,None]*delay[None,:,:])).reshape(len(bins)*32,-1)
    norms=np.linalg.norm(dictionary,axis=0);dictionary/=norms
    ah=np.ascontiguousarray(dictionary.conj().T)
    def data(ranges,angles,coeff,seed):
        y=np.zeros((len(bins),32),complex)
        for r,a,c in zip(ranges,angles,coeff):
            xx=r*np.sin(np.deg2rad(a));zz=r*np.cos(np.deg2rad(a));dist=np.hypot(xx-xpos,zz)
            y+=c*sf[:,None]*np.sinc(cfg.element_width_lambda/cfg.fc*f[:,None]*(xx-xpos)[None,:]/dist[None,:])*np.exp(-2j*np.pi*f[:,None]*(r+dist[None,:]-2*cfg.gate_range)/cfg.c)
        # 20 dB per selected frequency-channel sample: a declared synthetic
        # noise normalization, NOT measured channel noise or identical raw SNR.
        sigma=float(np.linalg.norm(y)/np.sqrt(y.size)/10)
        gen=np.random.default_rng(seed);noise=(gen.standard_normal(y.shape)+1j*gen.standard_normal(y.shape))/np.sqrt(2)
        return (y+sigma*noise).ravel(),sigma
    def solve(y,sigma):
        residual=y.copy();selected=[];steps=[];coef=np.empty(0);start=time.perf_counter()
        # Discrepancy rule fixed before test data; no target count as input.
        tolerance=1.3*sigma*np.sqrt(len(y))
        for it in range(32):
            score=abs(ah@residual);score[selected]=0;idx=int(score.argmax());selected.append(idx)
            coef=np.linalg.lstsq(dictionary[:,selected],y,rcond=None)[0];residual=y-dictionary[:,selected]@coef
            steps.append({'iteration':it+1,'selected_dictionary_index':idx,'residual_norm':float(np.linalg.norm(residual)),'stop_tolerance':tolerance})
            if np.linalg.norm(residual)<=tolerance:break
        im=np.zeros((len(rr),len(aa)),float);im.ravel()[selected]=abs(coef/norms[selected])
        return im,steps,time.perf_counter()-start
    scenarios={'independent_training':(np.array([50.36,50.64]),np.array([-18.,27.])),
               'offgrid_test':(np.array([50.3537,50.4213,50.4967,50.5731,50.6459]),np.array([-39.37,-17.23,.37,24.13,46.41]))}
    rows=[];log=[];saved={};truth_table=[]
    for scene,(tr,ta) in scenarios.items():
        coeff=np.exp(1j*np.random.default_rng(2026106100+(scene=='offgrid_test')).uniform(-np.pi,np.pi,len(tr)))/np.sqrt(len(tr))
        for rep in range(3):
            y,sigma=data(tr,ta,coeff,2026106200+100*(scene=='offgrid_test')+rep);im,steps,elapsed=solve(y,sigma)
            pd,fp,npeaks=_blind_detection(im,rr,aa,.05*im.max(),tr,ta)
            rows.append({'scene':scene,'noise_replication':rep,'method':'32-channel random array + joint wideband OMP','detection_probability':pd,'false_peaks':fp,'blind_peak_count':npeaks,'omp_iterations':len(steps),'solver_seconds_local_CPU':elapsed,'sigma_frequency_samples':sigma,'dictionary_bytes':dictionary.nbytes,'conjugate_dictionary_bytes':ah.nbytes})
            log.extend({'scene':scene,'noise_replication':rep,**s} for s in steps)
            if rep==0:saved[scene]=im
        truth_table.extend({'scene':scene,'range_m':float(r),'angle_deg':float(a),'amplitude':float(abs(c)),'phase_rad':float(np.angle(c))} for r,a,c in zip(tr,ta,coeff))
    # Principal architectures keep their implemented Tx banks. This supplement
    # isolates a meaningful competing receiver/reconstructor, not a strict
    # waveform/Tx/noise-matched proof of architectural superiority.
    tr,ta=scenarios['offgrid_test'];coeff=np.exp(1j*np.random.default_rng(2026106101).uniform(-np.pi,np.pi,5))/np.sqrt(5)
    baseline={}
    for m in METHODS:
        raw,mf,_=point_raw(cfg,m,tr,ta,coeff,geometry='spherical');baseline[m]=reconstruct_points(cfg,m,raw,mf,rr,aa)
        pd,fp,npk=_blind_detection(baseline[m],rr,aa,.05*abs(baseline[m]).max(),tr,ta)
        rows.append({'scene':'offgrid_test_noise_free_reference','noise_replication':0,'method':m,'detection_probability':pd,'false_peaks':fp,'blind_peak_count':npk,'omp_iterations':None,'solver_seconds_local_CPU':None,'sigma_frequency_samples':0.,'dictionary_bytes':None,'conjugate_dictionary_bytes':None})
    write_csv(out/'metrics.csv',rows);write_csv(out/'omp_residual_trace.csv',log);write_csv(out/'continuous_truth.csv',truth_table)
    np.savez_compressed(out/'raw_results.npz',ranges_m=rr,angles_deg=aa,random_positions_m=xpos,selected_baseband_hz=nu[bins],**saved,**{m:abs(im).astype(np.float32) for m,im in baseline.items()})
    fig,axs=plt.subplots(1,4,figsize=(14,3.5))
    for ax,(label,im) in zip(axs,list(baseline.items())+[('Random + OMP',saved['offgrid_test'])]):
        z=ax.pcolormesh(aa,rr,db(im,floor=-40),cmap='viridis',vmin=-40,vmax=0,shading='auto',rasterized=True)
        ax.scatter(ta,tr,facecolors='none',edgecolors='white',s=24);ax.set(title=LABELS.get(label,label),xlabel='Angle (deg)',ylabel='Range (m)');fig.colorbar(z,ax=ax,label='Amplitude (dB)')
    savefig(fig,out/'supplementary_reconstructions')
    dump(out/'scope.json',{'candidate':'32 physical nonuniform receive channels over same 211.7mm aperture, endpoints fixed; ideal broad point Tx, full-band Hamming LFM','appropriate_algorithm':'joint wideband range-angle OMP from actual frequency-channel measurements','source_count_given_to_solver':False,'stop_rule':'residual <= 1.3 sigma sqrt(M), cap 32 atoms; fixed before independent training/test','dictionary_grid':{'range_step_m':.01,'angle_step_deg':.5,'atoms':len(norms),'frequency_bins':41,'frequency_bin_selection_seed':2026106002},'generator':'continuous off-grid truth, spherical receive distances; dictionary uses far-field delays','training':'separate two-target scene; no test-truth parameter adjustment','test_scope':'small sparse off-grid scene, not diffuse/shipwreck recovery','reference_comparison':'principal three configurations plotted noise-free; competitive OMP uses 20dB selected-frequency noise. Not a controlled superiority ranking across unequal transmit/noise conditions','detection_scope':'5%-peak diagnostic only in this supplement; independently calibrated primary threshold test is R1-23','timing_scope':'local CPU supplementary solver only; NOT Jetson or end-to-end throughput'})

def model_validation(out):
    """Fixed numerical correctness checks, before interpreting experiments."""
    rows=[];cfg=Config(angular_nodes=101,interpolation_oversampling=8)
    q=np.linspace(-15,15,501);n=128;explicit=np.sum(weights(n)[:,None]*np.exp(1j*(np.arange(n)-(n-1)/2)[:,None]*q[None,:]),axis=0)
    err=float(np.max(abs(explicit-array_factor(q,n))))
    rows.append({'test':'closed_Hamming_AF_vs_explicit_element_sum','error':err,'tolerance':1e-10,'passed':err<=1e-10})
    r=cfg.ranges[np.argmin(abs(cfg.ranges-50.5))];a=cfg.angles[63];truth=grid_catalog(cfg,[r],[a],[1])
    ir=np.where((cfg.ranges>=50.3)&(cfg.ranges<=50.7))[0];ia=np.arange(0,101,2)
    for m in METHODS:
        im,_=image_grid(cfg,truth,m);raw,mf,meta=point_raw(cfg,m,[r],[a],geometry='farfield')
        direct=reconstruct_points(cfg,m,raw,mf,cfg.ranges[ir],cfg.angles[ia]);err=float(np.max(abs(im[ir[:,None],ia]-direct))/np.max(abs(im)))
        rows.append({'test':m+'_angular_convolution_vs_independent_element_DAS','error':err,'tolerance':.002,'passed':err<=.002})
        fine=reconstruct_points(replace(cfg,interpolation_oversampling=16),m,raw,mf,cfg.ranges[ir],cfg.angles[ia]);err=float(np.max(abs(fine-direct))/max(np.max(abs(fine)),1e-30))
        rows.append({'test':m+'_Fourier_interpolation_8_vs_16','error':err,'tolerance':.002,'passed':err<=.002})
    raw,mf,_=point_raw(cfg,'dense',[50.5037],[13.37],geometry='spherical')
    rr=np.arange(50.4937,50.5138,.001);aa=np.arange(12.37,14.371,.05);im=reconstruct_points(cfg,'dense',raw,mf,rr,aa,focus_model='spherical');i=np.unravel_index(abs(im).argmax(),im.shape)
    er=float(abs(rr[i[0]]-50.5037));ea=float(abs(aa[i[1]]-13.37))
    rows.extend([{'test':'spherical_point_range_carrier_counted_once','error':er,'tolerance':.0011,'passed':er<=.0011},{'test':'spherical_point_angle_carrier_counted_once','error':ea,'tolerance':.0501,'passed':ea<=.0501}])
    w,_=waveforms(cfg,'dense');time_energy=float(np.sum(abs(w)**2));fft_energy=float(np.sum(abs(np.fft.fft(w,8192,axis=0))**2)/8192);err=abs(time_energy-fft_energy)
    rows.append({'test':'prototype_Parseval_energy','error':err,'tolerance':1e-12,'passed':err<=1e-12})
    write_csv(out/'metrics.csv',rows);dump(out/'validation.json',{'all_passed':all(r['passed'] for r in rows),'tests':rows,'thresholds':'fixed numerical tolerances, not favorable performance gates','bounded_scope':'verifies this software model, not measured sonar or all possible scenes','separate_discretization':'large-scene grid model and exact continuous-point channel model cross-checked at grid-aligned point'})
    if not all(r['passed'] for r in rows):raise RuntimeError('Model correctness gate failed; preserve diagnostics and stop interpreting downstream results')
