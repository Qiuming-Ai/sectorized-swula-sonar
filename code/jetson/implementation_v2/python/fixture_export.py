"""Generate frozen common-input Jetson correctness fixtures, without hardware.

Usage: python fixture_export.py --out ABSOLUTE_DIR
An existing completed output directory is never overwritten. Use --noise-count
only for an explicitly labelled development subset; formal defaults remain20.
"""
from __future__ import annotations
import argparse
import importlib.util
import json
import time
from pathlib import Path
import numpy as np
from fixture_lib import (MODEL_SHA256, METHODS, FIXED_GATES, comparison, dump,
                         independent_image, matched_time, read_method, receive_fft, sha)

HERE = Path(__file__).resolve().parent
MODEL = HERE.parents[1] / "reference_snapshot" / "numerical_model.py"

def load_model():
    import sys
    if sha(MODEL) != MODEL_SHA256:
        raise RuntimeError("Frozen reference model hash mismatch")
    spec = importlib.util.spec_from_file_location("jetson_frozen_model_v2", MODEL)
    nm = importlib.util.module_from_spec(spec); sys.modules[spec.name] = nm
    spec.loader.exec_module(nm)
    return nm

def save_array(path, values, dtype, manifest):
    values = np.ascontiguousarray(values, dtype=np.dtype(dtype))
    values.tofile(path)
    manifest.append(dict(path=str(path), dtype=values.dtype.str, shape=list(values.shape),
                         bytes=path.stat().st_size, sha256=sha(path)))

def bank_files(nm, cfg, method, root, manifest):
    path = root / method / "static"; path.mkdir(parents=True)
    _, q, _, metadata = nm.spectra(cfg, method, 16384)
    g = nm.frontend_response(cfg, 16384)
    positions = nm.positions(cfg, method)
    ranges = np.linspace(30., 30.8, 321)
    angles = np.deg2rad(np.linspace(-50., 50., 401))
    sector = np.zeros(401, dtype=np.int32) if method=="dense" else np.clip(np.digitize(np.rad2deg(angles),cfg.edges[1:-1]),0,5).astype(np.int32)
    meta = dict(schema_version=2, channels=len(positions), record_samples=8192,
                fft_samples=16384, oversampling=8, range_rows=321, angle_columns=401,
                waveforms=q.shape[1], fs=400000., fc=450000., c=1500., gate_range=28.5)
    (path / "meta.txt").write_text("".join(f"{k}={v}\n" for k,v in meta.items()), encoding="utf-8")
    for name, arr, dtype in (("G.bin",g,"<c8"),("G64.bin",g,"<c16"),
            ("Q.bin",q.T,"<c8"),("Q64.bin",q.T,"<c16"),
            ("positions.bin",positions,"<f8"),("angles.bin",angles,"<f8"),
            ("ranges.bin",ranges,"<f8"),("sector_map.bin",sector,"<i4")):
        save_array(path/name, arr, dtype, manifest)
    dump(path/"waveform_metadata.json", metadata)
    return path, read_method(path)

def emitter(nm, cfg, method, ranges, angles, coefficients):
    raw_fft, _, _ = nm.point_raw(cfg, method, ranges, angles, coefficients,
                                geometry="spherical", nfft=8192)
    # Quantization precedes every independent CPU/GPU receive reconstruction.
    return np.fft.ifft(raw_fft, axis=0).T.astype(np.complex64)

def write_fixture(nm, cfg, method, bank, path, name, raw, manifest, catalog=None):
    folder = path.parent/name; folder.mkdir(parents=True,exist_ok=True)
    t0 = time.perf_counter()
    save_array(folder/"raw.bin", raw, "<c8", manifest)
    reference = independent_image(raw, bank)
    save_array(folder/"reference.bin", reference, "<c16", manifest)
    # Secondary consistency check only: not the source of reference.bin.
    raw_fft = np.fft.fft(np.asarray(raw,np.complex128).T, n=16384, axis=0)
    frozen = nm.reconstruct_points(cfg, method, raw_fft, bank["q"].T,
                                   bank["ranges"], np.rad2deg(bank["angles"]))
    audit = comparison(reference, frozen, bank)
    audit.update(independent_reference_seconds=time.perf_counter()-t0,
                 secondary_model_used_to_generate_reference=False,
                 raw_quantized_before_CPU_reference=True)
    dump(folder/"cpu_consistency.json", audit)
    if catalog is not None:
        dump(folder/"catalog.json", catalog)
    print(f"{method}/{name}: CPU {audit['relative_l2']} passed={audit['passed']} {audit['independent_reference_seconds']:.2f}s", flush=True)
    return dict(method=method, name=name, directory=str(folder),
                bank_directory=str(path), cpu_consistency=audit,
                benchmark_workload=name in ("paper_clean",) or name.startswith("noise"))

def write_override_unit(nm,cfg,method,base_bank,path,name,template,manifest):
    """Fixture-local override verifies convolution orientation, not main physics."""
    from dataclasses import replace
    folder=path.parent/name;folder.mkdir(parents=True)
    meta=dict(base_bank["meta"]);meta["unit_template_override"]=1
    (folder/"meta.txt").write_text("".join(f"{k}={v}\n" for k,v in meta.items()),encoding="utf-8")
    n=int(meta["fft_samples"]); k=int(meta["waveforms"])
    q=np.fft.fft(template/np.linalg.norm(template),n)
    # Each waveform has the same deliberately asymmetric template; sector
    # coverage remains exercised, but no physical ZC/LFM interpretation applies.
    q=np.tile(q[None,:],(k,1));g=np.ones(n,dtype=np.complex128)
    for file,arr,dtype in (("G.bin",g,"<c8"),("G64.bin",g,"<c16"),
            ("Q.bin",q,"<c8"),("Q64.bin",q,"<c16"),
            ("positions.bin",base_bank["positions"],"<f8"),
            ("angles.bin",base_bank["angles"],"<f8"),
            ("ranges.bin",base_bank["ranges"],"<f8"),
            ("sector_map.bin",base_bank["sector_map"],"<i4")):
        save_array(folder/file,arr,dtype,manifest)
    bank=read_method(folder)
    raw=np.zeros((int(meta["channels"]),8192),dtype=np.complex64)
    for channel in range(len(raw)):
        raw[channel,1090+channel%7]=(.45+.005*channel)*np.exp(1j*(.17+.13*channel))
        raw[channel,1140+channel%5]=(.25-.003j*channel)*np.exp(.37j)
    save_array(folder/"raw.bin",raw,"<c8",manifest)
    reference=independent_image(raw,bank)
    save_array(folder/"reference.bin",reference,"<c16",manifest)
    fe=receive_fft(raw,bank);mf=matched_time(fe,bank,0)
    save_array(folder/"fe_reference.bin",fe,"<c16",manifest)
    save_array(folder/"mf0_reference.bin",mf,"<c16",manifest)
    # Same discrete definition cross-check with FE disabled in frozen routine.
    frozen=nm.reconstruct_points(replace(cfg,rx_frontend=False),method,
        np.fft.fft(raw.astype(np.complex128).T,n=n,axis=0),q.T,
        bank["ranges"],np.rad2deg(bank["angles"]))
    audit=comparison(reference,frozen,bank)
    dump(folder/"cpu_consistency.json",audit)
    dump(folder/"catalog.json",dict(unit_only=True,benchmark_workload=False,
        physical_waveform=False,frontend_override="identity",template_samples=len(template),
        template_real=template.real.tolist(),template_imag=template.imag.tolist(),
        matched_definition="IFFT(FFT(raw)*conjugate FFT(template)); unshifted lag axis",
        mf_reference_waveform=0))
    print(f"{method}/{name}: override CPU {audit['relative_l2']} passed={audit['passed']}",flush=True)
    return dict(method=method,name=name,directory=str(folder),bank_directory=str(folder),
                cpu_consistency=audit,benchmark_workload=False,unit_template_override=True)

def asymmetric_unit_tests(out):
    """Independent orientation/length tests, no override of physical fixtures."""
    cases = []
    arrays = {}
    for name, h in (("odd",np.array([1+2j,-.4+.2j,.7-1j,.1+.4j,-.3-.7j])),
                    ("even",np.array([.7-.5j,-.8+.9j,.2+.3j,1.1-.2j]))):
        x = np.zeros(37,complex); x[13]=.5+.8j; x[20]=-.3+.4j
        # q is a waveform, not a pre-reversed/pre-conjugated FIR kernel.
        q = h/np.linalg.norm(h)
        full = np.convolve(x, np.conj(q[::-1]), mode="full")
        spectrum = np.fft.ifft(np.fft.fft(x,64)*np.conj(np.fft.fft(q,64)))
        expected = np.zeros(64,complex)
        for lag in range(-len(q)+1,len(x)):
            expected[lag % 64] = full[lag+len(q)-1]
        error=float(np.linalg.norm(spectrum-expected)/np.linalg.norm(expected))
        arrays.update({f"{name}_raw":x,f"{name}_template":q,
                       f"{name}_expected_correlation":expected})
        cases.append(dict(name=name,template_samples=len(q),relative_l2=error,
                          passed=error<1e-12,scope="independent unit test, not a benchmark workload"))
    np.savez(out/"asymmetric_convolution_units.npz",**arrays)
    dump(out/"asymmetric_convolution_units.json",cases)
    return cases

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--out",type=Path,default=HERE/"generated")
    parser.add_argument("--noise-count",type=int,default=20)
    parser.add_argument("--methods",nargs="+",choices=METHODS,default=list(METHODS))
    args=parser.parse_args()
    if not 0<=args.noise_count<=20:
        raise ValueError("noise-count must be0..20")
    out=args.out.resolve()
    if out.exists() and any(out.iterdir()):
        raise FileExistsError(f"Refusing to overwrite existing fixture output: {out}")
    out.mkdir(parents=True,exist_ok=True)
    nm=load_model(); cfg=nm.Config()
    files, fixtures=[],[]; began=time.perf_counter()
    paper_r=np.repeat([30.2,30.4,30.6],6)
    paper_a=np.tile([-40.,-24.,-8.,8.,24.,40.],3)
    coefficient=np.exp(1j*np.random.default_rng(2026101001).uniform(-np.pi,np.pi,18))/np.sqrt(18)
    # Derive common per-element sigma once from Dense clean whole-record power.
    dense_clean=emitter(nm,cfg,"dense",paper_r,paper_a,coefficient)
    dense_power=float(np.mean(np.abs(dense_clean.astype(np.complex128))**2))
    sigma=float(np.sqrt(dense_power/10**(-5/10)))
    dump(out/"protocol.json",dict(schema_version=2,frozen_model_path=str(MODEL),
        frozen_model_sha256=MODEL_SHA256,correctness_gates=FIXED_GATES,
        formal_noise_count=20,generated_noise_count=args.noise_count,
        formal_complete=args.noise_count==20 and tuple(args.methods)==METHODS,
        noise_sigma_complex=float(sigma),dense_clean_whole_record_power=dense_power,
        noise_definition="Each complex pre-FE sample has E|eta|^2=sigma^2; identical128-channel prefix, seed2026101001+i; Dense clean raw whole block SNR=-5dB",
        processing="raw complex64->complex128 FFT16k->G64 once->conjugate Q64->IFFT->complex Fourier resample8->linear delay sample->carrier phase->normalized Hamming DAS",
        no_fitted_gain=True,matched_Q_contains_common_Tx_and_Rx_response=True,
        FE_stage="channel-major complex spectrum FFT(raw padded16k)*G",
        MF_stage="channel-major complex time record IFFT(FE*conjQ[k])",
        nonzero_Nyquist_interpolation="scipy complex resample: split even Nyquist into positive/negative half",
        paper_coordinates=dict(ranges_m=paper_r.tolist(),angles_deg=paper_a.tolist(),
            coefficients_real=coefficient.real.tolist(),coefficients_imag=coefficient.imag.tolist())))
    for method in args.methods:
        path,bank=bank_files(nm,cfg,method,out,files)
        clean=dense_clean if method=="dense" else emitter(nm,cfg,method,paper_r,paper_a,coefficient)
        fixtures.append(write_fixture(nm,cfg,method,bank,path,"paper_clean",clean,files))
        for i in range(args.noise_count):
            rng=np.random.default_rng(2026101001+i)
            shared=(rng.standard_normal((128,8192))+1j*rng.standard_normal((128,8192)))/np.sqrt(2)
            noisy=(clean.astype(np.complex128)+sigma*shared[:len(clean)]).astype(np.complex64)
            fixtures.append(write_fixture(nm,cfg,method,bank,path,f"noise{i:02d}",noisy,files,
                dict(seed=2026101001+i,sigma_complex=sigma,shared_128_channel_prefix=True)))
        fixtures.append(write_fixture(nm,cfg,method,bank,path,"zero",np.zeros_like(clean),files))
        edges=np.linspace(-50.,50.,7)
        bc=np.exp(1j*np.linspace(.2,2.3,7))/np.sqrt(7)
        boundary=emitter(nm,cfg,method,np.full(7,30.4),edges,bc)
        fixtures.append(write_fixture(nm,cfg,method,bank,path,"boundary_points",boundary,files,
            dict(ranges_m=[30.4]*7,angles_deg=edges.tolist(),note="seven sector edges including both FOV endpoints")))
        offgrid=emitter(nm,cfg,method,[30.5037],[24.13],[np.exp(.7j)])
        fixtures.append(write_fixture(nm,cfg,method,bank,path,"offgrid",offgrid,files,
            dict(ranges_m=[30.5037],angles_deg=[24.13],coefficient_phase_rad=.7)))
        impulse=np.zeros_like(clean)
        for channel in range(len(impulse)):
            impulse[channel,1100+channel%17]=(.25+channel/len(impulse))*np.exp(1j*(.3+.11*channel))
        fixtures.append(write_fixture(nm,cfg,method,bank,path,"complex_impulse",impulse,files,
            dict(note="channel-dependent complex impulse; interpolation, conjugation and carrier-phase unit stress test")))
        for name,template in (("asymmetric_odd",np.array([1+2j,-.4+.2j,.7-1j,.1+.4j,-.3-.7j])),
                              ("asymmetric_even",np.array([.7-.5j,-.8+.9j,.2+.3j,1.1-.2j]))):
            fixtures.append(write_override_unit(nm,cfg,method,bank,path,name,template,files))
    units=asymmetric_unit_tests(out)
    # Register all text/JSON and unit assets after their final writes.
    registered={Path(x["path"]) for x in files}
    for path in sorted(out.rglob("*")):
        if path.is_file() and path not in registered:
            files.append(dict(path=str(path),dtype=None,shape=None,bytes=path.stat().st_size,sha256=sha(path)))
    manifest=dict(schema_version=2,created_elapsed_seconds=time.perf_counter()-began,
                  source=dict(path=str(MODEL),sha256=MODEL_SHA256),files=files,fixtures=fixtures,
                  total_bytes=sum(x["bytes"] for x in files),fixture_count=len(fixtures),
                  independent_unit_tests=units,all_cpu_model_checks_passed=all(x["cpu_consistency"]["passed"] for x in fixtures))
    dump(out/"manifest.json",manifest)
    print(json.dumps(dict(output=str(out),fixture_count=len(fixtures),total_bytes=manifest["total_bytes"],
                         all_cpu_model_checks_passed=manifest["all_cpu_model_checks_passed"],elapsed_s=manifest["created_elapsed_seconds"])),flush=True)

if __name__=="__main__":
    main()
