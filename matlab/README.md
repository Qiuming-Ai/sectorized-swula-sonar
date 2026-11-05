# MATLAB source guide

The repository packages the study-related MATLAB source tree in separate modules. Add only the module or experiment you are running to the MATLAB path: some historical functions have names shared with other implementations.

| Module | Contents and starting point |
|---|---|
| [Numbered revision experiments](paperA_revise/README.md) | A1–A19 archive, shared numerical helpers, and experiment-specific commands. A18 is a Python-only schematic. Start with the A1 quick run below. |
| [Legacy SonarSim](legacy_sonarsim/README.md) | Original Paper A experiment scripts, baseband processing functions, and desktop MATLAB timing. Use the optional [archived timing inputs](../data/desktop-timing/README.md) for an exact-input replay. |
| [Shipwreck chain](paperA_revise/11_A11_shipwreck_common_scene/README.md) | A11 synthetic-scene builder and reconstruction, with its A9 audit and A16 phase experiment. Extract the [shipwreck data](../data/shipwreck/README.md) asset. See the A11 guide for using archived arrivals or compiling the MATLAB/CUDA MEX to regenerate them from STL. |
| [Figure prototypes](figure_prototypes/README.md) | MATLAB curve and schematic examples used during figure development, including the Figure 2 slice replot. |

From the repository root, with a licensed MATLAB installation:

```sh
matlab -batch "addpath('matlab/paperA_revise'); run_paperA_revision_experiments('quick',1)"
matlab -batch "addpath('matlab/figure_prototypes'); PSF"
```

Use the module notes for toolbox and CUDA prerequisites. MATLAB products and prebuilt MEX binaries are not bundled. Generated results and the companion data overlays are not part of the Git tree.

The historical MATLAB and final Python implementations differ in some scene, waveform, and noise settings. Use the [evidence map](../docs/evidence-map.md) to reproduce a specific published result, and use the C++/CUDA module to reproduce the Jetson benchmark. Regenerating arrivals with the historical CUDA ray tracer does not promise identical record order or bitwise reproduction of the archived catalogue.
