"""Independent, auditable simulation utilities. No paper/source files are edited.

Complex baseband is relative to fc. Exact paths use exp(-j*2*pi*(fc+nu)*tau)
ONCE. TTD is an additional physical Tx excitation delay, not a carrier phase
coefficient substituted for an envelope delay. Point echo generation is exact
spherical; large scenes use a separately labelled far-field angular-grid model.
"""
from __future__ import annotations
from dataclasses import dataclass, asdict, replace
from functools import lru_cache
import csv, hashlib, json, time
from pathlib import Path
import numpy as np
from scipy import signal, ndimage, fft
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
PROJECT=ROOT.parent
METHODS=('dense','cdm','fdm')
LABELS={'dense':'Dense ULA','cdm':'CDM SW-ULA','fdm':'packed FDM SW-ULA'}
COLORS={'dense':'#0072B2','cdm':'#D55E00','fdm':'#009E73'}
plt.rcParams.update({'font.size':10,'axes.labelsize':10,'legend.fontsize':9,'figure.dpi':110,'savefig.dpi':300,'axes.spines.top':False,'axes.spines.right':False})

@dataclass(frozen=True)
class Config:
    c:float=1500.0
    fc:float=450000.0
    fs:float=400000.0
    bandwidth:float=100000.0
    duration:float=0.005
    n_tx:int=12
    n_dense:int=128
    n_sparse:int=32
    focus:tuple=(-41.25,-24.25,-8.,8.,24.25,41.25)
    centres:tuple=(426000.,442000.,490000.,474000.,458000.,412000.)
    roots:tuple=(5,7,1,3,11,9)
    code_length:int=499
    chip_rate:float=100000.0
    subband:float=20000.0
    element_width_lambda:float=0.45
    fov:tuple=(-50.,50.)
    range_origin:float=30.0
    range_span:float=0.8
    leading_guard_m:float=1.5
    angular_nodes:int=201
    interpolation_oversampling:int=8
    seed:int=2026100901
    energy_mode:str='tx_element_unit_energy'
    tx_full:bool=False
    rx_frontend:bool=True
    @property
    def dr(self):return self.c/(2*self.fs)
    @property
    def gate_range(self):return self.range_origin-self.leading_guard_m
    @property
    def edges(self):return np.linspace(*self.fov,len(self.focus)+1)
    @property
    def u(self):return np.linspace(np.sin(np.deg2rad(self.fov[0])),np.sin(np.deg2rad(self.fov[1])),self.angular_nodes)
    @property
    def angles(self):return np.rad2deg(np.arcsin(self.u))
    @property
    def ranges(self):return self.gate_range+np.arange(int(np.ceil((self.range_span+2*self.leading_guard_m)/self.dr))+1)*self.dr

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def dump(path,obj):
    Path(path).parent.mkdir(parents=True,exist_ok=True)
    Path(path).write_text(json.dumps(obj,ensure_ascii=False,indent=2,allow_nan=False,default=lambda x:x.item() if isinstance(x,np.generic) else str(x))+'\n',encoding='utf-8')
def write_csv(path,rows):
    rows=list(rows)
    if not rows:return
    Path(path).parent.mkdir(parents=True,exist_ok=True)
    with Path(path).open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def savefig(fig,path):
    fig.tight_layout();fig.savefig(str(path)+'.png',bbox_inches='tight');fig.savefig(str(path)+'.svg',bbox_inches='tight');plt.close(fig)
def db(x,ref=None,floor=-80):
    a=np.abs(x);r=np.max(a) if ref is None else ref
    return np.maximum(20*np.log10(np.maximum(a,1e-30)/max(float(r),1e-30)),floor)

def positions(cfg,method,tx=False):
    n=cfg.n_tx if tx else (cfg.n_dense if method=='dense' else cfg.n_sparse)
    d=(0.5 if tx or method=='dense' else 2.0)*cfg.c/cfg.fc
    return (np.arange(n)-(n-1)/2)*d
def weights(n):
    w=np.hamming(n);return w/w.sum()

def dirichlet(q,n):
    q=np.asarray(q);den=np.sin(q/2);out=np.empty_like(q,dtype=float);ok=np.abs(den)>1e-10
    np.divide(np.sin(n*q/2),den,out=out,where=ok)
    out[~ok]=n*np.cos(n*q[~ok]/2)/np.cos(q[~ok]/2)
    return out
def array_factor(q,n):
    """Signed centred Hamming array response, including correct 2pi limits."""
    if n==1:return np.ones_like(q)
    a=2*np.pi/(n-1)
    return (.54*dirichlet(q,n)+.23*(dirichlet(q-a,n)+dirichlet(q+a,n)))/np.hamming(n).sum()
def tx_pattern(cfg,freq,angles,centres=None,steering='ttd'):
    f=np.asarray(freq).reshape(-1,1,1);u=np.sin(np.deg2rad(np.asarray(angles))).reshape(1,1,-1)
    u0=np.sin(np.deg2rad(cfg.focus if centres is None else centres)).reshape(1,-1,1)
    fs=f if steering=='ttd' else cfg.fc
    q=2*np.pi*(cfg.c/cfg.fc/2)/cfg.c*(f*u-fs*u0)
    return array_factor(q,cfg.n_tx)*np.sinc(cfg.element_width_lambda/cfg.fc*f*u)

@lru_cache(maxsize=32)
def frontend_taps(cfg):
    return signal.firwin(129,cfg.bandwidth/2,fs=cfg.fs,window=('kaiser',8.6))

@lru_cache(maxsize=32)
def frontend_response(cfg,nfft):
    """Known 64-sample group delay compensated; linear-record guards required."""
    if not cfg.rx_frontend:return np.ones(nfft,dtype=complex)
    nu=np.fft.fftfreq(nfft,1/cfg.fs)
    return np.fft.fft(frontend_taps(cfg),nfft)*np.exp(2j*np.pi*nu*64/cfg.fs)

def waveforms(cfg,method,taper=None,shaping=True,length=None,roots=None,tx_full=None):
    n=int(round(cfg.duration*cfg.fs));t=np.arange(n)/cfg.fs;tau=t-t.mean()
    if method=='dense':
        x=np.exp(1j*np.pi*(cfg.bandwidth/cfg.duration)*tau**2)[:,None]
    elif method=='fdm':
        off=np.asarray(cfg.centres)-cfg.fc
        x=np.exp(1j*np.pi*(cfg.subband/cfg.duration)*tau[:,None]**2+2j*np.pi*tau[:,None]*off)
    else:
        l=cfg.code_length if length is None else length;rt=cfg.roots if roots is None else roots
        chip=np.arange(l)[:,None]
        if any(np.gcd(int(r),l)!=1 for r in rt):raise ValueError('ZC root must be coprime with L')
        frequency_zc=np.exp(-1j*np.pi*chip*(chip+1)*np.asarray(rt)[None,:]/l)
        chips=np.fft.ifft(frequency_zc,axis=0)
        chips=chips/np.sqrt(np.mean(np.abs(chips)**2,axis=0,keepdims=True))
        spc=int(round(cfg.fs/cfg.chip_rate))
        if abs(spc-cfg.fs/cfg.chip_rate)>1e-12:raise ValueError('integer resampling factor required')
        x=signal.resample_poly(chips,spc,1,axis=0,window=('kaiser',5.0))
        x=np.pad(x,((0,max(0,n-l*spc)),(0,0)))[:n]
    actual_taper=('rectangular' if method=='cdm' else 'hamming') if taper is None else taper
    if actual_taper=='hamming':x=x*np.hamming(n)[:,None]
    elif actual_taper!='rectangular':raise ValueError(actual_taper)
    retained_fraction=1.0
    full_length=len(x)
    use_full=cfg.tx_full if tx_full is None else tx_full
    if shaping:
        h=frontend_taps(cfg)
        complete=np.stack([signal.fftconvolve(x[:,k],h,mode='full') for k in range(x.shape[1])],axis=1)
        full_length=len(complete)
        cropped=complete[64:64+n]
        retained_fraction=float(np.sum(np.abs(cropped)**2)/np.sum(np.abs(complete)**2))
        x=complete if use_full else cropped
    x=x/np.sqrt(np.sum(np.abs(x)**2,axis=0,keepdims=True))
    return x,{'method':method,'taper':actual_taper,'shaping':('common 129-tap Kaiser FIR; full sensitivity record' if use_full else 'common 129-tap Kaiser FIR; 64-sample delay compensation and common centre crop') if shaping else 'no common Tx FIR (diagnostic only)','nominal_duration_s':cfg.duration,'stored_support_s':len(x)/cfg.fs,'full_convolution_samples':full_length,'crop_retained_energy_fraction_before_normalization':retained_fraction,'common_FIR_group_delay_samples':64 if shaping else 0,'energy_per_prototype':np.sum(np.abs(x)**2,axis=0).tolist(),'code_length':cfg.code_length if length is None else length,'roots':list(cfg.roots if roots is None else roots),'sampling_hz':cfg.fs,'cdm_generation':'frequency-domain ZC -> IFFT -> Kaiser-5 polyphase resample; no temporal taper' if method=='cdm' else None}

def spectra(cfg,method,nfft,taper=None,shaping=True):
    w,meta=waveforms(cfg,method,taper,shaping)
    s=np.fft.fft(w,n=nfft,axis=0);send=s/np.sqrt(s.shape[1]);nu=np.fft.fftfreq(nfft,1/cfg.fs);f=cfg.fc+nu
    if method=='dense':energy=1.0
    else:
        xm=positions(cfg,method,True);wt=weights(cfg.n_tx)
        exc=send[:,None,:]*np.exp(-2j*np.pi*f[:,None,None]*xm[None,:,None]*np.sin(np.deg2rad(cfg.focus))[None,None,:]/cfg.c)
        energy=float(np.sum(np.abs(np.sum(exc,axis=2)*wt[None,:])**2)/nfft)
    if cfg.energy_mode=='tx_element_unit_energy':send=send/np.sqrt(energy)
    elif cfg.energy_mode!='bank_unit_energy':raise ValueError(cfg.energy_mode)
    q=s*frontend_response(cfg,nfft)[:,None]
    q_energy=np.sum(np.abs(q)**2,axis=0)/nfft
    matched=q/np.sqrt(q_energy[None,:])
    meta.update({'energy_mode':cfg.energy_mode,'unscaled_tx_element_excitation_energy':energy,'applied_send_gain':1/np.sqrt(energy) if cfg.energy_mode=='tx_element_unit_energy' else 1.0,'final_tx_element_excitation_energy':1.0 if cfg.energy_mode=='tx_element_unit_energy' else energy,'matched_filter':'unit-energy actual Tx pulse through common Rx frontend','rx_frontend':'same delay-compensated 129-tap FIR on guarded continuous record' if cfg.rx_frontend else 'disabled diagnostic','rx_template_energy_before_unit_normalization':q_energy.tolist()})
    return send,matched,nu,meta

def point_raw(cfg,method,ranges,angles,amplitudes=None,geometry='spherical',errors=None,shaping=True,nfft=None):
    ranges=np.atleast_1d(ranges);angles=np.atleast_1d(angles);p=len(ranges)
    coeff=np.ones(p,dtype=complex) if amplitudes is None else np.asarray(amplitudes,dtype=complex)
    err={} if errors is None else errors
    nfft=8192 if nfft is None else nfft
    send,matched,nu,meta=spectra(cfg,method,nfft,shaping=shaping);f=cfg.fc+nu
    xn=positions(cfg,method);xm=np.array([0.]) if method=='dense' else positions(cfg,method,True)
    nr=len(xn);rng=np.random.default_rng(cfg.seed)
    eps={k:rng.standard_normal(cfg.n_dense)[:nr] for k in ['position','gain','phase','sync']}
    xn_true=xn+float(err.get('position_lambda',0))*cfg.c/cfg.fc*eps['position']
    cg=cfg.c*(1+float(err.get('sound_speed_fraction',0)))
    g=10**(float(err.get('gain_sd_db',0))*eps['gain']/20)*np.exp(1j*np.deg2rad(float(err.get('phase_sd_deg',0)))*eps['phase'])
    sync=float(err.get('sync_sd_us',0))*1e-6*eps['sync']
    width=cfg.element_width_lambda*(1+float(err.get('element_width_fraction',0)))
    out=np.zeros((nfft,nr),complex)
    for rr,aa,ap in zip(ranges,angles,coeff):
        u=np.sin(np.deg2rad(aa));xx=rr*u;zz=rr*np.cos(np.deg2rad(aa))
        dt=np.sqrt((xx-xm)**2+zz**2) if geometry=='spherical' else rr-xm*u
        dr=np.sqrt((xx-xn_true)**2+zz**2) if geometry=='spherical' else rr-xn_true*u
        tr=(dr/cg-cfg.gate_range/cfg.c)+sync
        rx=np.exp(-2j*np.pi*f[:,None]*tr[None,:])*g[None,:]
        arrival_u=(xx-xn_true)/dr if geometry=='spherical' else np.full(nr,u)
        rx*=np.sinc(width/cfg.fc*f[:,None]*arrival_u[None,:])
        if method=='dense':illum=send[:,0]*np.exp(-2j*np.pi*f*(dt[0]/cg-cfg.gate_range/cfg.c))
        else:
            tt=dt/cg-cfg.gate_range/cfg.c
            wt=weights(cfg.n_tx)
            alltx=np.zeros((nfft,len(cfg.focus)),complex)
            for k,focus in enumerate(cfg.focus):
                steering=xm*np.sin(np.deg2rad(focus+float(err.get('tx_steering_deg',0))))/cfg.c
                phase=np.exp(-2j*np.pi*f[:,None]*(tt[None,:]+steering[None,:]))
                depart_u=(xx-xm)/dt if geometry=='spherical' else np.full(len(xm),u)
                ep=np.sinc(width/cfg.fc*f[:,None]*depart_u[None,:])
                alltx[:,k]=np.sum(phase*ep*wt[None,:],axis=1)
            illum=np.sum(send*alltx,axis=1)
        ripple=float(err.get('transducer_ripple_db',0));phase_r=float(err.get('transducer_phase_deg',0))
        response=10**((ripple*np.sin(2*np.pi*nu/cfg.bandwidth))/20)*np.exp(1j*np.deg2rad(phase_r)*np.sin(2*np.pi*nu/cfg.bandwidth))
        out+=ap*illum[:,None]*rx*response[:,None]
    return out,matched,meta

def reconstruct_points(cfg,method,raw_fft,matched,ranges,angles,focus_model='farfield',rx_steering_deg=0):
    ranges=np.asarray(ranges);angles=np.asarray(angles);xn=positions(cfg,method);wr=weights(len(xn));nfft=len(raw_fft)
    u=np.sin(np.deg2rad(angles+rx_steering_deg));zz=ranges[:,None]*np.cos(np.deg2rad(angles))[None,:];xx=ranges[:,None]*u[None,:]
    os=cfg.interpolation_oversampling;image=np.zeros((len(ranges),len(angles)),complex)
    sec=np.clip(np.digitize(angles,cfg.edges[1:-1]),0,len(cfg.focus)-1)
    for k in range(matched.shape[1]):
        ix=np.arange(len(angles)) if method=='dense' else np.where(sec==k)[0]
        if not len(ix):continue
        mf=np.fft.ifft(raw_fft*frontend_response(cfg,nfft)[:,None]*np.conj(matched[:,k,None]),axis=0)
        # Fourier interpolation of a zero-padded linear matched-filter record.
        fine=signal.resample(mf,nfft*os,axis=0)
        scan=np.zeros((len(ranges),len(ix)),complex)
        for n,pos in enumerate(xn):
            if focus_model=='spherical':tau=(ranges[:,None]+np.sqrt((xx[:,ix]-pos)**2+zz[:,ix]**2))/cfg.c-2*cfg.gate_range/cfg.c
            else:tau=2*(ranges[:,None]-cfg.gate_range)/cfg.c-pos*u[ix][None,:]/cfg.c
            idx=tau*cfg.fs*os
            if np.any(idx<0) or np.any(idx>=len(fine)-1):raise ValueError('Scan outside guarded linear record')
            low=np.floor(idx).astype(int);frac=idx-low
            val=fine[low,n]*(1-frac)+fine[low+1,n]*frac
            scan+=wr[n]*val*np.exp(2j*np.pi*cfg.fc*tau)
        image[:,ix]=scan
    return image

def scene_grid(cfg,n,layers,extent=None,dynamic_db=0,seed=None,shape='random'):
    rng=np.random.default_rng(cfg.seed if seed is None else seed)
    extent=cfg.range_span if extent is None else extent
    # Continuous angles; physically fixed region and explicit finite layer support.
    if shape=='wall':a=rng.uniform(-30,30,n);r=cfg.range_origin+extent*(.5+.12*np.sin(np.deg2rad(a*2)))
    elif shape=='two_clusters':a=rng.choice([-20.,20.],n)+rng.normal(0,4,n);r=cfg.range_origin+rng.choice([.2,.8],n)*extent
    else:a=rng.uniform(-45,45,n);r=cfg.range_origin+np.linspace(.1,.9,layers)[rng.integers(layers,size=n)]*extent
    amp=10**(-rng.uniform(0,dynamic_db,n)/20);amp=amp/np.sqrt(np.sum(amp**2))
    phase=rng.uniform(-np.pi,np.pi,n);coeff=amp*np.exp(1j*phase)
    return grid_catalog(cfg,r,a,coeff),{'ranges_m':r,'angles_deg':a,'coefficients':coeff,'seed':cfg.seed if seed is None else seed,'count':n,'layers':layers,'extent_m':extent,'dynamic_range_db':dynamic_db,'shape':shape,'amplitude_norm':'sum of individual reflectivity powers = 1'}

def grid_catalog(cfg,r,a,coeff):
    grid=np.zeros((len(cfg.ranges),len(cfg.u)),complex)
    ri=np.rint((np.asarray(r)-cfg.gate_range)/cfg.dr).astype(int)
    ai=np.rint((np.sin(np.deg2rad(a))-cfg.u[0])/(cfg.u[1]-cfg.u[0])).astype(int)
    valid=(ri>=0)&(ri<len(grid))&(ai>=0)&(ai<grid.shape[1])
    np.add.at(grid,(ri[valid],ai[valid]),np.asarray(coeff)[valid])
    return grid

def image_grid(cfg,scene,method,shaping=True,energy_mode=None):
    """Broadband far-field DAS on a uniform sin(theta) grid.

    Linear angular convolution is algebraically the explicit uniform-array sum,
    not a receive grating mask. Only the assigned Tx sector output is selected.
    Frequency and angular convolutions are padded. Grid discretization is kept
    separate from exact continuous-coordinate point validation.
    """
    if energy_mode is not None:cfg=replace(cfg,energy_mode=energy_mode)
    nr,na=scene.shape;nw=len(waveforms(cfg,method,shaping=shaping)[0]);nf=fft.next_fast_len(nr+2*nw-2)
    send,matched,nu,meta=spectra(cfg,method,nf,shaping=shaping)
    f=cfg.fc+nu;carrier=np.exp(-4j*np.pi*cfg.fc*(cfg.ranges-cfg.gate_range)/cfg.c)
    sc=np.fft.fft(scene*carrier[:,None],n=nf,axis=0)
    y=np.zeros_like(sc);sec=np.clip(np.digitize(cfg.angles,cfg.edges[1:-1]),0,5)
    lag=np.arange(-na+1,na)*(cfg.u[1]-cfg.u[0]);spatial_n=fft.next_fast_len(3*na-2)
    n=cfg.n_dense if method=='dense' else cfg.n_sparse;d=(.5 if method=='dense' else 2)*cfg.c/cfg.fc
    for lo in range(0,nf,128):
        sl=slice(lo,min(nf,lo+128));ff=f[sl]
        if method=='dense':illum=send[sl,0,None]*np.ones((1,na))
        else:illum=np.sum(send[sl,:,None]*tx_pattern(cfg,ff,cfg.angles),axis=1)
        incident=sc[sl]*illum*np.sinc(cfg.element_width_lambda/cfg.fc*ff[:,None]*cfg.u[None,:])*frontend_response(cfg,nf)[sl,None]
        kernel=array_factor(2*np.pi*ff[:,None]*d/cfg.c*lag[None,:],n)
        scan=np.fft.ifft(np.fft.fft(incident,n=spatial_n,axis=1)*np.fft.fft(kernel,n=spatial_n,axis=1),axis=1)[:,na-1:na-1+na]
        if method=='dense':y[sl]=scan*np.conj(matched[sl,0,None])
        else:y[sl]=scan*np.conj(matched[sl][:,sec])
    im=np.fft.ifft(y,axis=0)[:nr]
    im*=np.exp(4j*np.pi*cfg.fc*(cfg.ranges-cfg.gate_range)/cfg.c)[:,None]
    return im,meta

def physical_guard(support,cfg,guard_m=.06,guard_deg=1.):
    """Exact rectangular coordinate guard, independent of angular indexing."""
    rows=int(np.floor((guard_m+1e-12)/cfg.dr))
    radial=ndimage.maximum_filter1d(support.astype(np.uint8),size=2*rows+1,axis=0,mode='constant')>0
    guard=np.zeros_like(radial)
    for j,a in enumerate(cfg.angles):
        ix=np.where(abs(cfg.angles-a)<=guard_deg+1e-12)[0]
        guard[:,j]=np.any(radial[:,ix],axis=1)
    return guard

def metrics(truth,image,cfg,dense=None,guard_m=.06,guard_deg=1.,fixed_background=None):
    t=np.abs(truth);p=np.abs(image)
    support=t>0
    guard=physical_guard(support,cfg,guard_m,guard_deg)
    # Use the declared physical evaluation range, not padded empty rows.
    evaluated=(cfg.ranges>=cfg.range_origin)&(cfg.ranges<=cfg.range_origin+cfg.range_span)
    valid=np.broadcast_to(evaluated[:,None],p.shape);background=valid&~guard
    alpha=np.vdot(p[valid],t[valid]).real/max(np.vdot(p[valid],p[valid]).real,1e-30)
    err=np.sum((alpha*p[valid]-t[valid])**2)/max(np.sum(t[valid]**2),1e-30)
    energy=np.sum(p[valid]**2);leak=np.sum(p[background]**2)/max(energy,1e-30)
    def rho(a,b):
        a=np.asarray(a).ravel();b=np.asarray(b).ravel();a=a-a.mean();b=b-b.mean();den=np.linalg.norm(a)*np.linalg.norm(b)
        return float(a@b/den) if den>0 else None
    bg_count=int(background.sum());bg_fraction=bg_count/int(valid.sum())
    core=support&valid
    cnr=(float(p[core].mean())-float(p[background].mean()))/max(float(p[background].std()),1e-30) if bg_count and core.any() else None
    fixed=valid&np.broadcast_to((abs(cfg.angles)>=47)[None,:],p.shape) if fixed_background is None else valid&fixed_background
    return {'truth_magnitude_nmse_db':float(10*np.log10(max(err,1e-30))),'truth_scale':float(alpha),'truth_magnitude_correlation':rho(t[valid],p[valid]),'dense_magnitude_correlation':rho(np.abs(dense)[valid],p[valid]) if dense is not None else None,'support_leakage':float(leak) if bg_fraction>=.05 else None,'cnr_linear_magnitude':cnr if bg_fraction>=.05 else None,'background_fraction':bg_fraction,'background_pixels':bg_count,'occupied_pixels':int(support.sum()),'occupied_range_cells':int(np.any(support,axis=1).sum()),'guard_range_m':guard_m,'guard_angle_deg':guard_deg,'metric_background_valid':bg_fraction>=.05,'metric_domain':'fixed declared physical range ROI, +/-50deg','fixed_background_pixels':int(fixed.sum()),'fixed_background_truth_pixels':int((fixed&support).sum()),'fixed_background_energy_fraction':float(np.sum(p[fixed]**2)/max(energy,1e-30)) if fixed.any() else None,'fixed_background_mean_power':float(np.mean(p[fixed]**2)) if fixed.any() else None}

def image_panel(cfg,images,path,truth=None,common=False,range_slice=None):
    keys=list(images);fig,axs=plt.subplots(1,len(keys),figsize=(4*len(keys),3.6),squeeze=False)
    valid=(cfg.ranges>=cfg.range_origin)&(cfg.ranges<=cfg.range_origin+cfg.range_span)
    ref=max(np.max(np.abs(x[valid])) for x in images.values()) if common else None
    for ax,k in zip(axs[0],keys):
        z=images[k][valid];m=ax.pcolormesh(cfg.angles,cfg.ranges[valid],db(z,ref,floor=-50),cmap='viridis',vmin=-50,vmax=0,shading='auto',rasterized=True)
        ax.set(xlabel='Angle (deg)',ylabel='Range (m)',title=LABELS.get(k,k));fig.colorbar(m,ax=ax,label='Amplitude (dB)')
    savefig(fig,path)

def peak_width(axis,power):
    axis=np.asarray(axis);p=np.asarray(power);i=int(np.argmax(p));half=p[i]/2
    l=i
    while l>0 and p[l]>=half:l-=1
    r=i
    while r<len(p)-1 and p[r]>=half:r+=1
    if l==0 and p[l]>=half:return None
    if r==len(p)-1 and p[r]>=half:return None
    left=axis[l]+(half-p[l])/(p[l+1]-p[l])*(axis[l+1]-axis[l])
    right=axis[r-1]+(half-p[r-1])/(p[r]-p[r-1])*(axis[r]-axis[r-1])
    return float(right-left)
