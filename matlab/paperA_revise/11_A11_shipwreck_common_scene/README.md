# A11 — Common-scene ray-traced shipwreck experiment in the legacy geometry

## Released MATLAB route

This is the **historical MATLAB control**, not the final Python common-FIR, equal total transmit-element-energy, common per-element-noise experiment. The archived geometry is shared, but the channel, sampling and metric definitions differ. For the exact manuscript Figure 13, use `code/scenes/plot_shipwreck.py` instead.

Extract the optional shipwreck-data archive at the package root so that `data/shipwreck/geometry/` contains `ship.stl` and `A11_raytrace_arrivals.csv`, or pass the absolute path of an equivalent extracted geometry directory. The default `archive` branch consumes the frozen, already ROI-filtered catalogue. It does **not** run or require CUDA, a saved MATLAB workspace, or a private SonarSim checkout. The STL is used to reproduce the diagnostic facet catalogue; it is not re-traced in this branch.

From the package root in MATLAB:

```matlab
addpath(fullfile(pwd,'matlab','paperA_revise','11_A11_shipwreck_common_scene'));
cfg = config_A11("paper");
[scene, visual, facets, seabed, arrivals, provenance] = build_shipwreck_scene(cfg);
metrics = run_A11_shipwreck_common_scene("paper", [], "archive");
```

Use `"quick"` for a lighter historical run: its archive branch thins the released, prefiltered rows and uses a coarser angle grid. It is not an exact replay of the older quick path, which thinned the unfiltered CUDA cache. Both modes retain the declared 42–122 m view. The default run writes to `outputs/matlab/A11_paper_archive/` or `outputs/matlab/A11_quick_archive/`. Optional arguments after `arrival_source` are an output directory and a locally built MEX directory. A run can create large workspaces, so place outputs on a writable disk with adequate space. The diagnostic plotter accepts `--results PATH_TO_RUN`.

For a fresh one-bounce STL-to-arrivals run, follow [`../../shipwreck_cuda/README.md`](../../shipwreck_cuda/README.md), then use `"retrace"` instead of `"archive"`. The regenerated catalogue and cache go to the run output folder, not into the released geometry inputs. CUDA append order can depend on execution scheduling, so the re-trace is not asserted to be byte-identical to the frozen archive.

## Historical model and evidence boundary

The three legacy `PaperA_exp8_*_real_scene.m` scripts use different active
models, arrays, waveforms, noise calls, and beamformers. Their
outputs cannot be compared directly. A11 retains the geometry and one-bounce
CUDA visibility/BRDF calculation from `PaperA_exp8_baseline_real_scene.m`
while replacing the method-dependent comparison chain with a shared
element-level baseband model.

The rebuilt scene uses the same `ship.stl`, scale 1, rotation
`[roll,pitch,yaw]=[90,0,-45] deg`, translation `[0,40,40] m`, and coordinate
permutation `[1,3,2]` as the legacy baseline. The sonar is at the origin. For
the one-dimensional receive array along x, the far-field phase is governed by
the direction cosine `x/r`, where `r=sqrt(x^2+y^2+z^2)`. The effective image
angle is therefore `asin(x/r)`, and the delay coordinate is the same full 3-D
slant range. Using the horizontal-plane angle `atan2(x,z)` would ignore the
unresolved y coordinate and artificially expand the hull laterally. The
45-deg direction stored by the preview helper is not the imaging axis because
the legacy `tx_cos_power` and `rx_cos_power` values are both -1 (isotropic gain).
The original STL is approximately 155 m long; after transformation it spans
roughly 46--117 m in slant range. The new display therefore uses 42--122 m
instead of copying the legacy script's inconsistent 4.5--6 m range entry.

The legacy seabed geometry is rebuilt as a 200 m by 100 m x-z grid at y=48 m
with 0.1 m vertical roughness. The complete ship and seabed meshes are passed
through the original one-bounce CUDA ray tracer. It produces 307,198 arrivals;
207,107 positive-amplitude arrivals remain within the paper window after the
direction-cosine projection. The
ray-tracer amplitudes are used directly, deterministic unresolved-surface
phases are assigned, and this identical complex scene is sent through the
dense, CDM, and FDM channel models. One paired raw-channel noise realization
is added at 10 dB before matched filtering and DAS imaging.

Panel (a) is a visual reference only. A 1.0-range-bin by 0.35-angle-bin
Gaussian smoothing and gamma=0.45 display compression make the geometry
legible at publication size. These display operations do not enter the
forward channel, reconstructions, or metrics.

This is a one-bounce ray-traced numerical experiment. It is not a water-tank
experiment or surveyed shipwreck data set. Elevation is not reconstructed,
and measured directivity, multipath, and reverberation remain absent. Because
the one-dimensional array cannot identify y and z separately, the remaining
curvature in an angle--slant-range display is the expected elevation-collapse
geometry, not deformation of the rigid STL mesh. No ground-truth-aided
dewarping is applied to the compared images.

Historical output types (run outputs now live below `outputs/matlab/`):

- `figures/A11_shipwreck_common_scene.pdf`: native direction-cosine-angle--slant-range comparison.
- `figures/A11_shipwreck_common_scene_sector.pdf`: historical receive-array sector diagnostic that informed the publication layout. The final manuscript panel is produced separately by the Python replotter; its second Cartesian axis is likewise an unresolved-plane radial coordinate, not physical forward z.
- `figures/A11_geometry_audit.pdf`: legacy 3-D pose and the exact 2-D projection.
- `results/A11_legacy_raytrace_cache.mat`: saved common CUDA arrival cloud.
- `results/A11_raytrace_arrivals.csv`: arrivals inside the paper window.
- `results/A11_ship_facets.csv`: transformed ship-facet geometry catalogue.
- `results/A11_seabed_scatterers.csv`: legacy-seabed geometry catalogue.
- `results/A11_shipwreck_metrics.csv`: common quantitative metrics.
