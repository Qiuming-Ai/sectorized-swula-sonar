# Jetson v2 correctness fixtures

This directory contains local fixture generation and an independently coded FP64 receive reference. It neither connects to Jetson nor changes the original CUDA project. The frozen reference source is hash-checked before generation.

## Frozen primary dimensions

The baseband sample rate is 400 kHz, reference carrier is 450 kHz, sound speed is1500 m/s, and gate range is28.5 m. Every input contains8192 complex samples/channel and is zero padded to16384 before receive filtering. The output grid is321 ranges from30 to30.8 m in0.0025 m increments and401 angles from-50 to50 degrees in0.25-degree increments. Dense has128 receive channels; CDM and FDM have32. Fourier interpolation is8x with subsequent complex linear interpolation.

All methods share the physical FIR policy. The effective unit-energy `Q` is the actual common-Tx-FIR pulse passed through the common receive frontend. `Q` is **not pre-conjugated**. The receive data undergo `G` exactly once; correlation uses `conj(Q)`.

## Layout

```
generated/
  protocol.json
  manifest.json
  dense|cdm|fdm/
    static/
      meta.txt
      G.bin, Q.bin             complex64 [N], [K,N]
      G64.bin, Q64.bin         complex128 [N], [K,N]
      positions.bin           float64 [M]
      angles.bin              float64 [A], radians
      ranges.bin              float64 [R], metres
      sector_map.bin          int32 [A]
    paper_clean|noise00..19|zero|boundary_points|offgrid|complex_impulse/
      raw.bin                 complex64 [M,8192]
      reference.bin           complex128 [R,A]
      cpu_consistency.json
    asymmetric_odd|asymmetric_even/
      same raw/reference files
      same-named complete static override bank
      fe_reference.bin        complex128 [M,N], spectral
      mf0_reference.bin       complex128 [M,N], time-domain k=0
```

All binary storage is C-major, little-endian. The two asymmetric template cases override the frontend with identity and replace the physical waveform with a deliberately non-symmetric4- or5-sample complex unit-energy sequence. They only test correlation orientation, even/odd indexing, interpolation and beamforming. They are not paper waveforms or performance workloads.

## Commands

Generate a new, non-existing output directory:

```
python fixture_export.py --out NEW_ABSOLUTE_DIRECTORY
```

The default creates27 cases/method (81 total), including20 fixed-seed noisy inputs. Common complex noise has `E|eta|²=sigma²` before the frontend; `sigma` is fixed from the Dense clean whole-record power and-5 dB SNR. CDM and FDM use the first32 channels of the same128-channel noise realization. Quantization to complex64 occurs before every CPU or GPU reconstruction.

Generate an on-demand stage reference:

```
python stage_reference.py --bank METHOD/static --raw METHOD/paper_clean/raw.bin --stage fe --out NEW_FE.bin
python stage_reference.py --bank METHOD/static --raw METHOD/paper_clean/raw.bin --stage mf --waveform 0 --out NEW_MF.bin
```

Compare GPU output without gain fitting:

```
python compare_outputs.py --bank METHOD/static --reference METHOD/paper_clean/reference.bin --candidate GPU_IMAGE.bin --kind image --out CHECK.json
```

`FE` denotes the spectral array `FFT(raw zero-padded)*G`; `MF[k]` denotes `IFFT(FE*conj(Q[k]))`. Stage arrays are channel-major. For asymmetric cases, use the case directory itself as `--bank`.

## Fixed correctness gates

- FE/MF complex relative L2 error≤1e-4.
- Final image complex relative L2 error≤1e-3, with no fitted complex gain or post-hoc rescaling.
- Nonconstant image magnitude Pearson correlation≥0.999.
- Zero input output maximum absolute value≤1e-7.
- Global peak differences≤one declared image bin (0.0025 m,0.25 degrees).

The image reference is independently coded NumPy FFT + scipy complex Fourier resampling + channel-wise fractional-delay sampling and carrier compensation. The frozen reconstruction is only a secondary consistency check and does not generate `reference.bin`. Complex Fourier resampling explicitly splits the even-length Nyquist bin into positive and negative half contributions, matching scipy's convention.
