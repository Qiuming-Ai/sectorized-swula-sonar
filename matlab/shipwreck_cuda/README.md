# Optional CUDA build: STL to A11 arrivals

This directory contains the original `stl_raytrace_arrivals_cuda3_5_1.cu` source and a local MATLAB build helper, **not a compiled MEX**. The source does not replace the released frozen arrival CSV. Neither the archived A11/A16 route nor the final manuscript image replot requires CUDA.

## Build and test

Use a MATLAB installation with `mexcuda`, an NVIDIA GPU/driver, and a C++ compiler and CUDA setup supported by that MATLAB release. Configure a supported C++ compiler with `mexcuda -setup C++` if none is selected. Do not copy MATLAB or CUDA runtime libraries into the release. The build uses MATLAB's `-R2018a` MEX interface and ordinary toolchain selection; it does not bypass compiler-compatibility checks.

From the repository root in MATLAB:

```matlab
addpath(fullfile(pwd,'matlab','shipwreck_cuda'));
binary = build_shipwreck_mex();
result = smoke_shipwreck_mex(binary);
```

The default binary is built into `outputs/matlab/shipwreck_cuda/`. You can instead pass a writable local build directory to `build_shipwreck_mex(out_dir)`. On Windows, if the compiler rejects a UNC path, use a local checkout/build directory. The smoke test traces a single triangle and checks the one-bounce output schema; it is not a full-geometry regeneration test.

## Regenerate the legacy arrival catalogue

Extract the geometry overlay at the repository root (or give its directory explicitly). Then add A11 to the MATLAB path and select the re-trace branch:

```matlab
addpath(fullfile(pwd,'matlab','paperA_revise','11_A11_shipwreck_common_scene'));
metrics = run_A11_shipwreck_common_scene("paper", [], "retrace", [], fileparts(binary));
```

For just the raw one-bounce arrival cloud, use `audit_A11_legacy_raytrace([], [], fileparts(binary))`. It retains the historical transform, procedural seabed, fixed random seeds and one-bounce CUDA model. A full A11 run projects the cloud using `sin(theta)=x/r`, applies the physical FOV/range filter and exports the new CSV beside the run outputs. `outputs/matlab/A11_paper_retrace/` is separate from the frozen archive. It is normal for GPU-generated arrival ordering to differ; use `"archive"` when comparing the frozen manuscript geometry. No measured sonar or calibrated scattering response is implied.

## Port verification

In MATLAB R2025b, the release source syntax check, A9 source audit, archived A11 scene construction, A16 configuration, and a small-grid Dense/CDM/FDM forward/reconstruction smoke all passed. The full archive scene used 207,107 arrivals on the 21,334 × 201 legacy grid; its visual reference agreed with the archive to a maximum absolute difference of `1.23e-15`. The full A11/A16 imaging runs were not repeated in this port check.

The local `mexcuda` attempt selected the NVIDIA CUDA compiler, but the bundled CUDA header rejected the installed MSVC host version (`C1189`). The compiler check was not bypassed. Build the MEX with a supported compiler/CUDA combination and run `smoke_shipwreck_mex` before a full re-trace. No compiled binary is included.
