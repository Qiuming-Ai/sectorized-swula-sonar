# Desktop MATLAB timing inputs

The optional `desktop-timing-data_v10.zip` asset contains the three exact archived synthetic-baseband HDF5 inputs used in the older desktop MATLAB timing experiment. Extract it at the repository root, then verify it with `python scripts/check_release.py --desktop-timing`. Do not commit the extracted `overlay/` folder.

Each file contains 100 synthetic pings. The associated 100-row recorded timing CSV and the MATLAB processing entry point are in [legacy SonarSim](../../matlab/legacy_sonarsim/README.md). This is a historical desktop processing measurement, not the Jetson replay workload. The three historical method inputs have different sampling/processing settings, SNR, and velocity; a rerun measures the new host and does not validate a controlled cross-platform speedup.
