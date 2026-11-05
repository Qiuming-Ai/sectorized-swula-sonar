"""Timestamp raw Jetson telemetry; integrate only bracketed, valid windows."""
from __future__ import annotations
import csv
import json
import re
import subprocess
import threading
import time
from pathlib import Path

POWER = re.compile(r"\bVDD_IN\s+(\d+)mW/(\d+)mW")
TEMP = re.compile(r"([A-Za-z0-9_]+)@([0-9.]+)C")

def parse_line(line, timestamp_ns):
    p = POWER.search(line)
    temps = {k: float(v) for k, v in TEMP.findall(line)}
    return {"timestamp_ns": timestamp_ns, "power_W": int(p[1]) / 1000 if p else None,
            "tool_average_W": int(p[2]) / 1000 if p else None,
            "max_temperature_C": max(temps.values()) if temps else None,
            "temperatures_C": temps}

def system_ram():
    values = {}
    for line in Path('/proc/meminfo').read_text().splitlines():
        k, v = line.split(':', 1)
        values[k] = int(v.strip().split()[0]) * 1024
    return values['MemTotal'] - values['MemAvailable']

class Collector:
    def __init__(self, directory, stop_file):
        self.directory, self.stop_file = Path(directory), Path(stop_file)
        self.rows, self.errors, self.raw = [], [], []
        self.process = self.thread = None
        self.on_safety = None

    def start(self):
        self.process = subprocess.Popen(['tegrastats', '--interval', '100'],
                                        stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                        text=True, bufsize=1)
        self.thread = threading.Thread(target=self._read, daemon=True)
        self.thread.start()

    def _read(self):
        for line in self.process.stdout:
            now = time.monotonic_ns()
            line = line.rstrip('\n')
            self.raw.append({"timestamp_ns": now, "line": line})
            row = parse_line(line, now)
            try:
                row['system_RAM_used_bytes'] = system_ram()
            except (OSError, KeyError, ValueError) as e:
                row['system_RAM_used_bytes'] = None
                self.errors.append(str(e))
            self.rows.append(row)
            if row['max_temperature_C'] is not None and row['max_temperature_C'] >= 80:
                self.stop_file.write_text('temperature >= 80 C\n')
                if self.on_safety is not None:
                    self.on_safety()

    def stop(self):
        if self.process is not None and self.process.poll() is None:
            self.process.terminate()  # Only our own telemetry child.
            try:
                self.process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                self.process.kill()
                self.process.wait()
        if self.thread is not None:
            self.thread.join(timeout=5)
        self.directory.mkdir(parents=True, exist_ok=True)
        with (self.directory/'power_raw.jsonl').open('x', encoding='utf-8') as f:
            for row in self.raw:
                f.write(json.dumps(row, ensure_ascii=False)+'\n')
        (self.directory/'power_samples.json').write_text(json.dumps(
            {"rows": self.rows, "reader_errors": self.errors,
             "requested_interval_ms": 100, "clock": "CLOCK_MONOTONIC arrival timestamps",
             "power_scope": "module VDD_IN current; tool average retained only as auxiliary"},
            indent=2, allow_nan=False)+'\n')

def integrate(rows, t0_ns, t1_ns, completed_blocks, workload_valid=True):
    reasons = []
    if not workload_valid: reasons.append('workload_failed')
    if t1_ns <= t0_ns or completed_blocks <= 0: reasons.append('invalid_window_or_count')
    if any(b['timestamp_ns'] <= a['timestamp_ns'] for a, b in zip(rows, rows[1:])):
        reasons.append('nonincreasing_timestamps')
    valid = [r for r in rows if r.get('power_W') is not None]
    before = [r for r in valid if r['timestamp_ns'] <= t0_ns]
    after = [r for r in valid if r['timestamp_ns'] >= t1_ns]
    if not before or not after: reasons.append('unbracketed_boundaries')
    result = {"valid": False, "reasons": reasons, "t0_ns": t0_ns, "t1_ns": t1_ns,
              "completed_blocks": completed_blocks, "idle_subtracted": False,
              "scope": "Jetson module VDD_IN, not carrier/system wall input"}
    if reasons: return result
    left, right = before[-1]['timestamp_ns'], after[0]['timestamp_ns']
    region = [r for r in rows if left <= r['timestamp_ns'] <= right]
    if any(r.get('power_W') is None for r in region): reasons.append('missing_power_field')
    if any(b['timestamp_ns']-a['timestamp_ns'] > 500_000_000 for a, b in zip(region, region[1:])):
        reasons.append('sample_gap_over_500_ms')
    if any(r.get('max_temperature_C') is not None and r['max_temperature_C'] >= 80 for r in region):
        reasons.append('temperature_stop')
    if reasons: return result
    def boundary(t):
        for a, b in zip(region, region[1:]):
            if a['timestamp_ns'] <= t <= b['timestamp_ns']:
                q = (t-a['timestamp_ns'])/(b['timestamp_ns']-a['timestamp_ns'])
                return {"timestamp_ns": t, "power_W": a['power_W']+q*(b['power_W']-a['power_W'])}
        if region[0]['timestamp_ns'] == t: return region[0]
        if region[-1]['timestamp_ns'] == t: return region[-1]
        raise ValueError('No bracket for integration boundary')
    points = [boundary(t0_ns)] + [r for r in region if t0_ns < r['timestamp_ns'] < t1_ns] + [boundary(t1_ns)]
    energy = sum((a['power_W']+b['power_W'])/2*(b['timestamp_ns']-a['timestamp_ns'])/1e9
                 for a, b in zip(points, points[1:]))
    duration = (t1_ns-t0_ns)/1e9
    result.update(valid=True, energy_J=energy, average_module_power_W=energy/duration,
                  energy_J_per_block=energy/completed_blocks, actual_window_seconds=duration,
                  valid_samples=len(region), maximum_interval_ms=max(
                      (b['timestamp_ns']-a['timestamp_ns'])/1e6 for a, b in zip(region, region[1:])))
    return result

def selftest():
    rows = [{"timestamp_ns": i*100_000_000, "power_W": 10+2*i/10,
             "max_temperature_C": 40} for i in range(11)]
    r = integrate(rows, 50_000_000, 950_000_000, 100)
    assert r['valid'] and abs(r['energy_J']-9.9) < 1e-12
    assert abs(r['average_module_power_W']-11) < 1e-12
    assert not integrate(rows, -1, 950_000_000, 100)['valid']
    assert not integrate(rows, 50_000_000, 950_000_000, 0)['valid']
    broken = [dict(r) for r in rows]; broken[5]['power_W'] = None
    assert 'missing_power_field' in integrate(broken, 50_000_000, 950_000_000, 100)['reasons']
    assert not integrate(rows[::-1], 50_000_000, 950_000_000, 100)['valid']
    assert not integrate(rows, 50_000_000, 950_000_000, 100, False)['valid']
    assert not integrate(rows[::6], 0, 600_000_000, 100)['valid']
    print('TELEMETRY_SELFTEST_PASS')

if __name__ == '__main__': selftest()
