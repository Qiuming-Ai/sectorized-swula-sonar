"""Portable release controller: fresh local build, freeze, and verification.

The historical paper records are read-only under data/jetson. A replay creates
its own binary identity, fixture manifest, implementation receipt, and logs.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
import statistics
import subprocess
import sys
import threading
import time
from pathlib import Path
from telemetry import Collector, integrate

ROOT = Path(__file__).resolve().parent
AUDIT = ROOT.parent
RELEASE = ROOT.parents[2]
HISTORICAL_FREEZE = RELEASE/'data'/'jetson'/'provenance'/'FROZEN_FORMAL_FINAL.json'
V2 = AUDIT/'implementation_v2'
V3 = AUDIT/'implementation_v3'
DATA = V2/'python'/'generated'
BANK = V3/'python'/'generated'/'d4_sinc16'
EXE = ROOT/'build'/'formal_runner'
FREEZE = ROOT/'FROZEN_FORMAL_FINAL.json'
CASES = sorted(['paper_clean','zero','boundary_points','offgrid','complex_impulse']+
               [f'noise{i:02d}' for i in range(20)])
METHODS = ('dense','cdm','fdm')

def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1_048_576), b''): h.update(block)
    return h.hexdigest()

def save(path, value, exclusive=False):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x' if exclusive else 'w', encoding='utf-8') as f:
        json.dump(value, f, indent=2, ensure_ascii=False, allow_nan=False); f.write('\n')

def fixture_reference(method, case):
    return BANK/'cases'/case/'reference.bin' if method=='fdm' else DATA/method/case/'reference.bin'

def freeze():
    # Must be run on the Jetson after compiling, before formal measurements.
    if not EXE.is_file(): raise RuntimeError('Compile the new runner before freezing')
    files = list((ROOT/'core').glob('*')) + list(ROOT.glob('*.py'))
    files += [ROOT/'contract.json', AUDIT/'README.md', AUDIT/'reference_snapshot'/'numerical_model.py']
    files += [V2/'core'/'engine.h', V2/'core'/'engine.cu', V3/'core'/'engine_v3.h', V3/'core'/'engine_v3.cu']
    files += [V2/'python'/'fixture_export.py', V2/'python'/'fixture_lib.py',
              V3/'python'/'export_references.py', V3/'python'/'fdm_reference.py']
    for method in METHODS:
        files += [p for p in (DATA/method/'static').iterdir() if p.is_file()]
        for case in CASES: files += [DATA/method/case/'raw.bin', fixture_reference(method,case)]
    files += [p for p in BANK.iterdir() if p.is_file() and p.suffix in ('.bin','.txt')]
    missing=[str(p) for p in files if not p.is_file()]
    if missing: raise RuntimeError('Missing required frozen files: '+str(missing))
    entries = [{'path':p.relative_to(AUDIT).as_posix(),'bytes':p.stat().st_size,'sha256':sha(p)}
               for p in sorted(set(files))]
    manifest = dict(files=entries, binary_sha256=sha(EXE), methods=list(METHODS), cases=CASES,
                    equivalence_claimed=False, adopted_engineering_candidate='d4_sinc16')
    manifest['record_kind']='new replay; not the original measured campaign'
    if HISTORICAL_FREEZE.is_file():
        manifest['historical_paper_manifest_sha256']=sha(HISTORICAL_FREEZE)
    save(FREEZE, manifest, True)
    write_permit(False)
    print(json.dumps(dict(event='frozen',files=len(entries),binary_sha256=manifest['binary_sha256'])), flush=True)

def verify_freeze():
    manifest=json.loads(FREEZE.read_text())
    if manifest['binary_sha256'] != sha(EXE): raise RuntimeError('Frozen binary changed')
    for row in manifest['files']:
        p = AUDIT/row['path']
        if p.stat().st_size!=row['bytes'] or sha(p)!=row['sha256']:
            raise RuntimeError('Frozen source/input changed: '+row['path'])
    return manifest

def write_permit(verified):
    # Separate licence for the selected engineering chain. Old image_valid stays false.
    values = dict(selected_engineering=1, equivalence_claimed=0, implementation_verified=int(verified),
                  candidate='d4_sinc16', frozen_manifest_sha256=sha(FREEZE),
                  binary_sha256=sha(EXE))
    if verified: values['formal_implementation_receipt_sha256']=sha(ROOT/'IMPLEMENTATION_RECEIPT_FORMAL.json')
    (ROOT/'permit.txt').write_text(''.join(f'{k}={v}\n' for k,v in values.items()))

def permission():
    verify_freeze()
    receipt=json.loads((ROOT/'IMPLEMENTATION_RECEIPT_FORMAL.json').read_text())
    if (not receipt['all_passed'] or receipt['comparison_count']!=75 or
        receipt['binary_sha256']!=sha(EXE) or receipt.get('frozen_manifest_sha256')!=sha(FREEZE)):
        raise RuntimeError('Current complete formal implementation receipt required')
    for row in receipt['output_files']:
        if sha(ROOT/row['path'])!=row['sha256']: raise RuntimeError('Changed formal validation output')
    write_permit(True)

def invoke(mode, method, session, out):
    out=Path(out); out.mkdir(parents=True, exist_ok=False)
    runner_out=out/'runner'; stop=out/'safety_stop.txt'
    cmd=[str(EXE),'--mode',mode,'--method',method,'--data',str(DATA),'--fdm-bank',str(BANK),
         '--permit',str(ROOT/'permit.txt'),'--session',str(session),'--out',str(runner_out),
         '--stop-file',str(stop)]
    collector=Collector(out,stop) if mode in ('power','memory') else None
    events=[]; lines=[]; code=None; watchdog=None; watchdog_expired=[]
    if collector: collector.start()
    try:
        p=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,bufsize=1)
        if collector:
            collector.on_safety=lambda: p.terminate() if p.poll() is None else None
        def timeout():
            if p.poll() is None:
                watchdog_expired.append(True)
                stop.write_text('whole runner safety timeout; per-call blocking duration not observed\n')
                p.terminate()  # Only our own workload child; raw telemetry and events remain.
                try:p.wait(timeout=5)
                except subprocess.TimeoutExpired:p.kill()
        watchdog=threading.Timer(300 if mode=='power' else 120,timeout)
        watchdog.daemon=True;watchdog.start()
        for line in p.stdout:
            lines.append(dict(arrival_ns=time.monotonic_ns(),line=line.rstrip()))
            print(line,end='',flush=True)
            try: events.append(json.loads(line))
            except json.JSONDecodeError: pass
        code=p.wait()
        if collector: time.sleep(0.65)  # Bracket t1; this is outside the measured workload.
    finally:
        if watchdog:watchdog.cancel()
        if collector: collector.stop()
        save(out/'invocation.json',dict(command=cmd,returncode=code,events=events,stdout=lines,
             whole_runner_timeout=bool(watchdog_expired)),True)
    if code!=0: raise RuntimeError(f'{mode}/{method}/{session} failed; all partial output retained')
    if mode=='power':
        window=json.loads((runner_out/'window.json').read_text())
        result=integrate(collector.rows,window['t0_ns'],window['t1_ns'],window['completed_blocks'],
                         window.get('workload_valid',False) and code==0 and not stop.exists())
        save(out/'POWER_INTEGRAL.json',result,True)
        if not result['valid']: raise RuntimeError('Invalid power window: '+str(result['reasons']))
    if mode=='memory':
        used=[r['system_RAM_used_bytes'] for r in collector.rows if r.get('system_RAM_used_bytes') is not None]
        save(out/'SYSTEM_RAM_SAMPLED.json',dict(sample_count=len(used),requested_interval_ms=100,
              sampled_max_used_bytes=max(used) if used else None,
              definition='MemTotal minus MemAvailable; sampled auxiliary, not continuous peak'),True)
    return runner_out

def checks():
    verify_freeze(); outputs=[]; all_rows=[]
    catalog=ROOT/'results'/('implementation_verification_'+time.strftime('%Y%m%d_%H%M%S'))
    for method in METHODS:
        out=invoke('verify',method,0,catalog/method)
        with (out/'verification.csv').open() as f: rows=list(csv.DictReader(f))
        if len(rows)!=25 or {r['fixture'] for r in rows}!=set(CASES): raise RuntimeError('Wrong verification catalog')
        if any(r['overall_pass']!='1' for r in rows): raise RuntimeError('Formal runner own-reference validation failed')
        all_rows.extend({**r,'method':method} for r in rows)
        outputs += [dict(path=p.relative_to(ROOT).as_posix(),sha256=sha(p)) for p in out.iterdir() if p.is_file()]
    receipt=dict(all_passed=True,comparison_count=75,rows=all_rows,output_files=outputs,
                 binary_sha256=sha(EXE),frozen_manifest_sha256=sha(FREEZE),
                 equivalence_claimed=False,record_kind='new replay implementation receipt')
    save(ROOT/'IMPLEMENTATION_RECEIPT_FORMAL.json',receipt,True);write_permit(True)
    print('FORMAL_IMPLEMENTATION_GATE_PASS 75/75; NOT fullrate FDM equivalence',flush=True)

def campaign(mode, run_id, only_session=None, only_method=None):
    permission()
    if not run_id or '/' in run_id or '\\' in run_id or run_id in ('.','..'):
        raise ValueError('Explicit safe run id required')
    contract=json.loads((ROOT/'contract.json').read_text())
    for session,order in enumerate(contract['latency']['method_orders']):
        if only_session is not None and session!=only_session: continue
        for method in order:
            if only_method is not None and method!=only_method: continue
            print(json.dumps(dict(event='collection_begin',mode=mode,method=method,session=session)),flush=True)
            invoke(mode,method,session,ROOT/'results'/run_id/mode/f'session{session}'/method)
            print(json.dumps(dict(event='collection_complete',mode=mode,method=method,session=session)),flush=True)
    verify_freeze()

def percentile(values,q):
    values=sorted(values);h=(len(values)-1)*q;i=math.floor(h);j=math.ceil(h)
    return values[i]+(h-i)*(values[j]-values[i])

def summarize(run_id):
    directory=ROOT/'results'/run_id; result=dict(run_id=run_id,formal=True,equivalence_claimed=False,latency={})
    for method in METHODS:
        values=[];sessions=[]
        for session in range(3):
            path=directory/'latency'/f'session{session}'/method/'runner'/'latency_raw.csv'
            if not path.exists(): sessions.append(dict(session=session,status='missing')); continue
            with path.open() as f: rows=list(csv.DictReader(f))
            measured=[r for r in rows if r['phase']=='measure']
            x=[float(r['wall_ms']) for r in measured if r['status']=='ok']
            if not all(math.isfinite(v) and v>=0 for v in x): raise RuntimeError('Invalid successful formal timing values')
            values+=x;sessions.append(dict(session=session,N=len(x),attempted=len(measured),
                failures=sum(r['status']!='ok' for r in measured),
                complete=len(x)==200 and len(measured)==200,
                mean_ms=statistics.mean(x) if x else None))
        summary=dict(N=len(values),complete=len(values)==600 and all(s.get('complete',False) for s in sessions),
                     failures=sum(s.get('failures',0) for s in sessions),sessions=sessions)
        if values: summary.update(mean_ms=statistics.mean(values),sample_SD_ms=statistics.stdev(values),
                                 median_ms=statistics.median(values),P95_ms=percentile(values,.95),max_ms=max(values))
        result['latency'][method]=summary
    result['power_windows']=[]
    for p in sorted((directory/'power').glob('session*/*/POWER_INTEGRAL.json')):
        result['power_windows'].append(dict(path=p.relative_to(ROOT).as_posix(),**json.loads(p.read_text())))
    result['power_complete']=len(result['power_windows'])==9 and all(r['valid'] for r in result['power_windows'])
    save(directory/'SUMMARY_FORMAL.json',result)
    print(json.dumps(result,indent=2,ensure_ascii=False),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['freeze','checks','latency','memory','power','summarize'])
    p.add_argument('--run-id');p.add_argument('--session',type=int,choices=range(3));p.add_argument('--method',choices=METHODS)
    a=p.parse_args()
    if a.mode=='freeze':freeze()
    elif a.mode=='checks':checks()
    elif a.mode=='summarize':summarize(a.run_id)
    else:campaign(a.mode,a.run_id,a.session,a.method)
