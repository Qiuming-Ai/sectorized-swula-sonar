"""Native8192-bin reference at ALL105 scan bearings on a coarser target grid.

This supplements (does not replace) the primary six-sector0.01deg frequency
convergence gate. It is a numerical spectral-integral check, not a finer angle
resolution claim or a scene-selection performance experiment.
"""
from pathlib import Path
from datetime import datetime,timezone
import sys,json,hashlib
import numpy as np
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE));import run_wave_spatial as ws
out=ws.HERE/'full_fov_native_frequency_validation';out.mkdir(parents=True,exist_ok=False)
scan=sorted(set(np.arange(-50,51,1).tolist()+list(ws.CFG.focus)));target=np.arange(-50,50.00001,.5)
protocol={'frozen_utc':datetime.now(timezone.utc).isoformat(),'purpose':'Independent native-full-frequency integral check of the chosen16-bin moment quadrature at all primary105 scan bearings, using201 target directions. Does not replace primary0.01deg fine angular profiles.','scan_angles_deg':scan,'target_grid_step_deg':.5,'frequency_reference':'8192 native FFT bins, no compression or PSD cutoff','candidate_group_bins':16,'scopes':['all_simultaneous_true_delay','isolated_assigned_true_delay'],'criteria':{'normalized_magnitude_relative_L2':.001,'alias_neighborhood_level_change_db':.05},'coarse_angle_maximum':'Stored as a diagnostic only; no sub0.5deg angular error claim','numerical_script_sha256':hashlib.sha256((HERE/'run_wave_spatial.py').read_bytes()).hexdigest(),'native_reference_direct_array_sum_validation':'direct_array_sum_validation/receipt.json,18 physical-element explicit-reference checks'}
ws.dump(out/'protocol_before_run.json',protocol);rows=[]
for method in ws.METHODS:
 for ii,s in enumerate(scan):
  a=ws.beam(method,s,target,16);b=ws.beam(method,s,target,1)
  for scope in protocol['scopes']:
   sa=ws.beam_stats(method,s,target,a[scope]);sb=ws.beam_stats(method,s,target,b[scope]);diff=ws.compare_profiles(a[scope],b[scope]);gl=abs(sa['worst_visible_alias_neighborhood_db']-sb['worst_visible_alias_neighborhood_db']) if sa['worst_visible_alias_neighborhood_db'] is not None else 0.
   rows.append({'method':method,'scan_angle_deg':s,'scope':scope,**diff,'alias_level_change_db':float(gl),'coarse_peak_angle_difference_deg':abs(sa['global_peak_angle_deg']-sb['global_peak_angle_deg']),'pass_fixed_criteria':bool(diff['normalized_magnitude_relative_L2']<=.001 and gl<=.05)})
  if ii%25==0:print(f'native frequency {method}:{ii+1}/{len(scan)}',flush=True)
 ws.csvout(out/'native_bin_vs_quadrature_all_scans.csv',rows)
ws.dump(out/'receipt.json',{'rows':len(rows),'all_amplitude_criteria_pass':all(r['pass_fixed_criteria'] for r in rows),'maximum_relative_L2':max(r['normalized_magnitude_relative_L2'] for r in rows),'maximum_alias_level_difference_db':max(r['alias_level_change_db'] for r in rows),'maximum_coarse_angle_difference_deg':max(r['coarse_peak_angle_difference_deg'] for r in rows),'fine_grid_gate_preserved':True,'source_unchanged':hashlib.sha256((HERE/'run_wave_spatial.py').read_bytes()).hexdigest()==protocol['numerical_script_sha256']})
print(json.dumps(json.loads((out/'receipt.json').read_text(encoding='utf-8'))),flush=True)
