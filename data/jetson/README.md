# Recorded Jetson benchmark

These are the formal measurement records used by the v08 manuscript's Jetson resource table and Figure 15, not the earlier v2/v3 pilot timings.

| Directory | Contents |
|---|---|
| `latency/session{0,1,2}/{dense,cdm,fdm}/runner/` | Nine raw latency CSVs and lifecycle records. Each CSV contains the 50 predefined warm-ups and 200 measured calls. |
| `power/session{0,1,2}/{dense,cdm,fdm}/` | Nine raw current/average VDD_IN streams with monotonic timestamps, parsed samples, integral receipts, runner windows, block-level status/endpoint CSVs, and sanitized phase events. Use the current reading, not the tool's cumulative average. |
| `memory/session{0,1,2}/{dense,cdm,fdm}/` | Nine explicit device-allocation ledgers, lifecycle snapshots, and auxiliary sampled RAM records. The inherited filename `power_samples.json` also holds the system RAM samples for this separate run; it is not used as a primary power window. |
| `summary/` | Timing and power summaries, allocation/FFT ledgers, and compact figure input CSVs. |
| `provenance/` | Original pre-measurement source/input/binary manifest and 75-case own-reference correctness receipt. These are archival, not live replay permits. |

The original source-relative record root was `jetson_audit_20261010/formal_collection_20261010/results/main_20261010/`. Numeric files are copied byte-for-byte. Each public `invocation.json` retains status and the runner's monotonic phase timestamps and counts; it omits original process IDs, commands, machine-local paths, and duplicated stdout.

From the repository root, `python code/jetson/verify_records.py` independently recalculates the table, including all 1,800 measured calls, 9 power windows, and 85,439 completed power-window blocks. Power is in W at the Jetson module's VDD_IN rail, energy per complete block is in J, and explicit allocation is in MiB.

No board executable, last-frame complex array, original GUI/DAQ record, or synthetic fixture bank is bundled. See `code/jetson/README.md` for the regenerators and the fresh-build procedure.
