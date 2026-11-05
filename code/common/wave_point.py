"""Frozen uniform-frontend waveform and continuous-point experiments.

Only writes its own results directory. The common numerical model is imported,
never changed. All adverse first outcomes and interpolation checks are retained.
"""
from __future__ import annotations
from dataclasses import asdict, replace
from pathlib import Path
from datetime import datetime, timezone, timedelta
import argparse, contextlib, hashlib, json, shutil, sys, time, traceback
import numpy as np
from scipy import signal, optimize
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numerical_model as nm

METHODS = ('dense', 'cdm', 'fdm')
LABELS = {'dense': 'Dense ULA', 'cdm': 'CDM SW-ULA', 'fdm': 'packed FDM SW-ULA'}
COLORS = {'dense': '#0072B2', 'cdm': '#D55E00', 'fdm': '#009E73'}
STYLES = {'dense': '-', 'cdm': '--', 'fdm': ':'}
MARKERS = {'dense': 'o', 'cdm': 's', 'fdm': '^'}
ROOT = Path(__file__).resolve().parents[1]
OUT_ROOT = ROOT / 'results' / 'wave_point'


def _power_db(x, ref=None, floor=-100):
    p = np.abs(x) ** 2
    refp = float(np.max(p)) if ref is None else float(ref)
    return np.maximum(10 * np.log10(np.maximum(p, 1e-30) / max(refp, 1e-30)), floor)


def _field_width(axis, amplitude):
    return nm.peak_width(np.asarray(axis), np.abs(amplitude) ** 2)


def _save_figure(fig, path):
    fig.tight_layout()
    fig.savefig(str(path) + '.png', dpi=300, bbox_inches='tight')
    fig.savefig(str(path) + '.svg', bbox_inches='tight')
    plt.close(fig)


def _series(ax, x, y, method, label=None, marker=True):
    ax.plot(x, y, color=COLORS[method], ls=STYLES[method],
            marker=MARKERS[method] if marker else None,
            markevery=max(1, len(x) // 8), linewidth=1.6,
            markersize=4, label=LABELS[method] if label is None else label)
    ax.grid(alpha=.17)


def _template_time(ft):
    """Rotate the padded causal/acausal record before taking its finite support.

    The Rx delay-compensated filter has samples at negative time indices. A
    direct [:N] cut would lose these samples. Discard only numerical roundoff
    below 1e-13 of peak amplitude and explicitly export its discarded energy.
    """
    x = np.fft.fftshift(np.fft.ifft(ft))
    valid = np.flatnonzero(np.abs(x) > 1e-13 * np.max(np.abs(x)))
    lo, hi = int(valid[0]), int(valid[-1]) + 1
    kept = x[lo:hi]
    dropped = float((np.sum(abs(x) ** 2) - np.sum(abs(kept) ** 2))
                    / max(np.sum(abs(x) ** 2), 1e-30))
    return kept, {'support_start_fftshift_index': lo, 'support_stop_fftshift_index': hi,
                  'discarded_roundoff_energy_fraction': max(0., dropped)}


def _autocorrelation(x, cfg, oversampling=16):
    ac = signal.correlate(x, x, mode='full', method='fft')
    fine = signal.resample(ac, len(ac) * oversampling)
    lag = np.arange(len(fine)) / (cfg.fs * oversampling) - (len(x) - 1) / cfg.fs
    return lag, fine


def _side_metrics(lag, ac):
    power = abs(ac) ** 2
    peak = int(np.argmax(power))
    minima, _ = signal.find_peaks(-power)
    left, right = minima[minima < peak], minima[minima > peak]
    if not len(left) or not len(right):
        return {'aperiodic_PSLR_db': None, 'aperiodic_ISLR_db': None,
                'first_minimum_mainlobe_available': False}
    a, b = int(left[-1]), int(right[0])
    side = np.r_[power[:a], power[b + 1:]]
    return {'aperiodic_PSLR_db': float(10 * np.log10(max(side.max(), 1e-30) / power[peak])),
            'aperiodic_ISLR_db': float(10 * np.log10(max(side.sum(), 1e-30) / power[a:b + 1].sum())),
            'first_minimum_mainlobe_available': True}


def _spectrum_statistics(x, cfg):
    nfft = 65536
    freq = np.fft.fftshift(np.fft.fftfreq(nfft, 1 / cfg.fs))
    power = abs(np.fft.fftshift(np.fft.fft(x, nfft))) ** 2
    total = float(power.sum())
    cumulative = np.cumsum(power) / total
    low = float(freq[np.searchsorted(cumulative, .005)])
    high = float(freq[min(np.searchsorted(cumulative, .995), len(freq) - 1)])
    return {'outside_nominal_100khz_fraction': float(power[abs(freq) > cfg.bandwidth / 2].sum() / total),
            'equal_tail_99pct_lower_hz': low, 'equal_tail_99pct_upper_hz': high,
            'equal_tail_99pct_bandwidth_hz': high - low}, freq, power


def _tx_energy_and_papr(cfg, method, send, nu):
    nfft = len(nu)
    nw = len(nm.waveforms(cfg, method)[0])
    if method == 'dense':
        excitation = send[:, :1]
        delay_min, delay_max = 0., 0.
    else:
        tx = nm.positions(cfg, method, True)
        wt = nm.weights(cfg.n_tx)
        delay = tx[:, None] * np.sin(np.deg2rad(cfg.focus))[None, :] / cfg.c
        delay_min, delay_max = float(delay.min()), float(delay.max())
        excitation = np.sum(send[:, None, :] * np.exp(
            -2j * np.pi * (cfg.fc + nu)[:, None, None] * delay[None, :, :]), axis=2) * wt[None, :]
    # Shift every element by the same delay so that intended TTD support is causal.
    excitation *= np.exp(2j * np.pi * nu[:, None] * delay_min)
    sequence = np.fft.ifft(excitation, axis=0)
    active_n = nw + int(np.ceil((delay_max - delay_min) * cfg.fs)) + 1
    active = abs(sequence[:active_n]) ** 2
    total_energy = float(np.sum(abs(sequence) ** 2))
    records = []
    for j in range(excitation.shape[1]):
        mean = float(active[:, j].mean())
        records.append({'method': method, 'tx_mode': 'full' if cfg.tx_full else 'primary',
                        'element': j + 1,
                        'element_energy_full_record': float(np.sum(abs(sequence[:, j]) ** 2)),
                        'active_window_samples': active_n,
                        'active_window_duration_ms': active_n / cfg.fs * 1000,
                        'active_mean_sample_power': mean,
                        'active_peak_sample_power': float(active[:, j].max()),
                        'active_window_PAPR_db': float(10 * np.log10(active[:, j].max() / max(mean, 1e-30))),
                        'energy_outside_declared_active_window_fraction': float(
                            np.sum(abs(sequence[active_n:, j]) ** 2) / max(np.sum(abs(sequence[:, j]) ** 2), 1e-30))})
    return total_energy, records, sequence, active_n


def run_waveforms(out):
    out.mkdir(parents=True, exist_ok=True)
    rows, pair_rows, energy_rows, bank_rows, arrays = [], [], [], [], {}
    plot_data = {}
    for mode in ('primary', 'full'):
        cfg = replace(nm.Config(), tx_full=(mode == 'full'))
        for method in METHODS:
            tx, meta = nm.waveforms(cfg, method)
            send, matched, nu, smeta = nm.spectra(cfg, method, 8192)
            actual_e, erec, excitation, active_n = _tx_energy_and_papr(cfg, method, send, nu)
            energy_rows += erec
            arrays[f'{mode}_{method}_tx_prototypes'] = tx
            arrays[f'{mode}_{method}_actual_Tx_elements'] = excitation
            bank_rows.append({'method': method, 'tx_mode': mode, 'Tx_support_samples': len(tx),
                              'Tx_support_ms': len(tx) / cfg.fs * 1000,
                              'actual_Tx_total_energy': actual_e,
                              'unscaled_Tx_total_energy': smeta['unscaled_tx_element_excitation_energy'],
                              'applied_send_gain': smeta['applied_send_gain'],
                              'crop_retained_energy_fraction_before_normalization': meta['crop_retained_energy_fraction_before_normalization'],
                              'actual_Tx_max_element_PAPR_db': max(r['active_window_PAPR_db'] for r in erec),
                              'PAPR_window_samples': active_n})
            templates = []
            for k in range(tx.shape[1]):
                q, qmeta = _template_time(matched[:, k])
                templates.append(q)
                arrays[f'{mode}_{method}_effective_Rx_template_{k + 1}'] = q
                for component, w in (('Tx_prototype', tx[:, k]), ('Tx_through_common_Rx', q)):
                    lag, ac = _autocorrelation(w, cfg)
                    spec, freq, power = _spectrum_statistics(w, cfg)
                    record = {'method': method, 'tx_mode': mode, 'channel': k + 1,
                              'component': component, 'support_samples': len(w),
                              'unit_template_energy': float(np.sum(abs(w) ** 2)),
                              'power_FWHM_mm': _field_width(lag * cfg.c / 2 * 1000, ac),
                              **spec, **_side_metrics(lag, ac),
                              'temporal_taper': meta['taper'],
                              'template_roundoff_discarded_energy_fraction': qmeta['discarded_roundoff_energy_fraction'] if component == 'Tx_through_common_Rx' else 0.}
                    rows.append(record)
                    arrays[f'{mode}_{method}_{component}_{k + 1}_autocorrelation'] = ac
                    arrays[f'{mode}_{method}_{component}_{k + 1}_lag_s'] = lag
                    if k == 0:
                        plot_data[(mode, method, component)] = (freq, power, lag, ac)
            # Both unfiltered input and filtered pulses are exported; no unfiltered
            # imaging trial is performed or treated as a final transmit design.
            native, _ = nm.waveforms(cfg, method, shaping=False)
            for k in range(native.shape[1]):
                stats, _, _ = _spectrum_statistics(native[:, k], cfg)
                rows.append({'method': method, 'tx_mode': mode, 'channel': k + 1,
                             'component': 'before_common_Tx_FIR_diagnostic',
                             'support_samples': len(native),
                             'unit_template_energy': float(np.sum(abs(native[:, k]) ** 2)),
                             'power_FWHM_mm': None, **stats,
                             'aperiodic_PSLR_db': None, 'aperiodic_ISLR_db': None,
                             'first_minimum_mainlobe_available': False,
                             'temporal_taper': meta['taper'],
                             'template_roundoff_discarded_energy_fraction': 0.})
            for component, bank in (('Tx_prototype', [tx[:, k] for k in range(tx.shape[1])]),
                                    ('Tx_through_common_Rx', templates)):
                for i in range(len(bank)):
                    for j in range(i + 1, len(bank)):
                        cross = signal.correlate(bank[i], bank[j], mode='full', method='fft')
                        lag = (np.arange(len(cross)) - (len(bank[j]) - 1)) / cfg.fs
                        norm = np.sqrt(np.sum(abs(bank[i]) ** 2) * np.sum(abs(bank[j]) ** 2))
                        amplitude = abs(cross) / norm
                        iglobal = int(np.argmax(amplitude))
                        local = abs(lag) <= .0006 + 1e-15
                        ilocal = np.where(local)[0][int(np.argmax(amplitude[local]))]
                        pair_rows.append({'method': method, 'tx_mode': mode, 'component': component,
                                          'channel_i': i + 1, 'channel_j': j + 1,
                                          'normalized_full_delay_cross_peak_db': float(20 * np.log10(amplitude[iglobal])),
                                          'global_peak_lag_ms': float(lag[iglobal] * 1000),
                                          'normalized_restricted_0p6ms_cross_peak_db': float(20 * np.log10(amplitude[ilocal])),
                                          'restricted_peak_lag_ms': float(lag[ilocal] * 1000)})
            nm.dump(out / f'{mode}_{method}_configuration.json', {'configuration': asdict(cfg), 'waveform': meta, 'spectra': smeta})
    nm.write_csv(out / 'all_channel_waveform_metrics.csv', rows)
    nm.write_csv(out / 'all_pair_cross_correlations.csv', pair_rows)
    nm.write_csv(out / 'Tx_bank_energy_and_support.csv', bank_rows)
    nm.write_csv(out / 'Tx_element_energy_and_PAPR.csv', energy_rows)
    np.savez_compressed(out / 'waveform_banks_and_correlations.npz', **arrays)
    fig, axs = plt.subplots(2, 3, figsize=(12, 6.3))
    for col, method in enumerate(METHODS):
        for mode, mode_style in (('primary', '-'), ('full', '--')):
            freq, power, lag, ac = plot_data[(mode, method, 'Tx_through_common_Rx')]
            axs[0, col].plot(freq / 1000, 10 * np.log10(np.maximum(power / power.max(), 1e-12)),
                             ls=mode_style, color=COLORS[method], linewidth=1.5,
                             label=f'{mode}: actual Tx + common Rx')
            axs[1, col].plot(lag * nm.Config().c / 2 * 1000, _power_db(ac),
                             ls=mode_style, color=COLORS[method], linewidth=1.5, label=mode)
        axs[0, col].axvline(-50, color='.4', ls=':', lw=.8)
        axs[0, col].axvline(50, color='.4', ls=':', lw=.8)
        axs[0, col].set(xlabel='Baseband frequency (kHz)', ylabel='Normalized power (dB)',
                        title=LABELS[method] + ', channel 1', ylim=(-100, 3), xlim=(-200, 200))
        axs[1, col].set(xlabel='Equivalent range lag (mm)', ylabel='Normalized correlation power (dB)',
                        xlim=(-160, 160), ylim=(-80, 3))
        for row in range(2):
            axs[row, col].grid(alpha=.17)
            axs[row, col].legend(fontsize=8)
    _save_figure(fig, out / 'uniform_waveforms_primary_vs_full')
    nm.dump(out / 'scope.json', {'energy': 'Actual simultaneous TTD transmit-element excitation energy=1; not unit waveform-bank energy.',
                                'CDM': 'frequency-domain ZC -> IFFT -> Kaiser5 polyphase resample, no temporal taper',
                                'FIR': 'all methods: common129 tap Kaiser8.6 cutoff50kHz; primary centred2000 sample crop; full2128 sample sensitivity control; common Rx frontend on guarded record',
                                'FWHM': 'half peak POWER of finite aperiodic autocorrelation, 16x Fourier interpolation; Tx-only and effective Rx templates separate',
                                'cross': 'All finite aperiodic pairs/lags; abs(lag)<=0.6ms separately, unit-energy prototypes.',
                                'PAPR': 'Peak/mean power within pulse support plus global TTD span; numerical fractional-delay tails outside this window explicitly counted; padded8192 record not used as PAPR mean.',
                                'pulse_duration': 'len/fs storage convention; last sample occurs one sample period earlier. full extends duration rather than pretending it is5ms.',
                                'raw_input': 'before_common_Tx_FIR is PSD diagnostic only, no unfiltered performance experiment',
                                'precision': 'Full float64/complex128 calculations and exports.'})
    return {'waveform_rows': len(rows), 'pair_rows': len(pair_rows), 'Tx_bank_rows': len(bank_rows), 'Tx_element_rows': len(energy_rows)}


def _blind_profile_resolution(r, profile, true_ranges):
    power = abs(profile) ** 2
    peaks, _ = signal.find_peaks(power, height=.25 * power.max())
    separation = float(true_ranges[1] - true_ranges[0])
    gate = max(.002, separation / 4)
    candidates = np.asarray(r)[peaks]
    matches = []
    used = set()
    for target in true_ranges:
        choices = [(abs(float(x) - target), k) for k, x in enumerate(candidates) if k not in used]
        if choices and min(choices)[0] <= gate:
            err, k = min(choices)
            used.add(k)
            matches.append((k, err))
    midpoint = (np.asarray(r) >= true_ranges[0]) & (np.asarray(r) <= true_ranges[1])
    valley = float(np.min(power[midpoint]) / power.max())
    resolved = len(matches) == 2 and valley < .8
    return {'blind_peak_count_above_minus6db_power': len(peaks),
            'matched_truth_count': len(matches),
            'false_candidates_above_threshold': len(peaks) - len(matches),
            'central_valley_to_global_peak_power': valley,
            'two_targets_resolved_under_declared_rule': resolved,
            'candidate_range_m': ';'.join(f'{v:.7f}' for v in candidates)}


def run_range(out):
    out.mkdir(parents=True, exist_ok=True)
    rows, arrays = [], {}
    rg = np.arange(30.32, 30.680001, .0005)
    for mode in ('primary', 'full'):
        cfg = replace(nm.Config(), tx_full=(mode == 'full'))
        for method in METHODS:
            for separation_mm in (5, 10, 20, 40, 80, 120):
                truth = np.asarray([30.5 - separation_mm / 2000, 30.5 + separation_mm / 2000])
                coefficients = np.asarray([1., np.exp(.7j)]) / np.sqrt(2)
                raw, matched, meta = nm.point_raw(cfg, method, truth, [8., 8.], coefficients, geometry='farfield')
                im = nm.reconstruct_points(cfg, method, raw, matched, rg, [8.])[:, 0]
                result = _blind_profile_resolution(rg, im, truth)
                rows.append({'method': method, 'tx_mode': mode, 'target_separation_mm': separation_mm,
                             'scatterer_phase_difference_rad': .7, 'target_angle_deg': 8.,
                             'true_range_1_m': float(truth[0]), 'true_range_2_m': float(truth[1]),
                             'range_grid_step_mm': .5, **result})
                arrays[f'{mode}_{method}_{separation_mm}mm'] = im
    arrays['range_m'] = rg
    np.savez_compressed(out / 'all_two_target_profiles.npz', **arrays)
    nm.write_csv(out / 'metrics.csv', rows)
    fig, axs = plt.subplots(1, 3, figsize=(12, 3.7))
    for ax, method in zip(axs, METHODS):
        for mode, ls in (('primary', '-'), ('full', '--')):
            profile = arrays[f'{mode}_{method}_80mm']
            ax.plot((rg - 30.5) * 1000, _power_db(profile), ls=ls,
                    color=COLORS[method], label=f'{mode}, separation 80 mm')
        for x in (-40, 40):
            ax.axvline(x, ls=':', color='.4', lw=.8)
        ax.set(xlabel='Range offset from 30.5 m (mm)', ylabel='Normalized power (dB)',
               title=LABELS[method], xlim=(-160, 160), ylim=(-45, 3))
        ax.legend(fontsize=8); ax.grid(alpha=.17)
    _save_figure(fig, out / 'physical_two_target_profiles')
    nm.dump(out / 'scope.json', {'model': 'Actual simultaneous TTD physical-frequency Tx and common Rx frontend, far-field point generation and DAS; complete waveform bank, not a single autocorrelation surrogate.',
                                'truth_phase': 'Two reflectivity coefficients differ by0.7rad. Propagation carrier phase is retained once; no truth-dependent phase compensation.',
                                'detector': 'Blind full fixed1D ROI local maxima above0.25 global maximum power, then one-to-one truth association within max(2mm,separation/4), central valley below0.8.',
                                'not_universal_resolution': 'Fixed phase, angle8deg, tested separations, fixed0.36m profile ROI; not a universal Rayleigh limit or full-FOV false-alarm experiment.',
                                'controls': 'Primary common5ms source and full5.32ms source both all3methods; no CDM temporal window.'})
    return {'rows': len(rows)}


def run_nearfield(out):
    out.mkdir(parents=True, exist_ok=True)
    rows, arrays = [], {}
    angle = np.arange(-50., 50.000001, .25)
    for distance in (5., 15., 30., 60.):
        cfg = replace(nm.Config(), range_origin=distance - .4)
        rg = np.arange(distance - .12, distance + .120001, .002)
        for target in (0., 40.):
            for method in METHODS:
                raw, matched, meta = nm.point_raw(cfg, method, [distance], [target], geometry='spherical')
                images = {focus: nm.reconstruct_points(cfg, method, raw, matched, rg, angle, focus_model=focus)
                          for focus in ('spherical', 'farfield')}
                ref = float(abs(images['spherical']).max())
                for focus, image in images.items():
                    ir, ia = np.unravel_index(np.argmax(abs(image)), image.shape)
                    profile = image[ir]
                    outside = abs(angle - target) > 8
                    peak = float(abs(image[ir, ia]))
                    rows.append({'method': method, 'distance_m': distance, 'target_angle_deg': target,
                                 'receive_focus': focus, 'angle_estimate_deg': float(angle[ia]),
                                 'range_estimate_m': float(rg[ir]),
                                 'angle_error_deg': float(angle[ia] - target),
                                 'range_error_mm': float((rg[ir] - distance) * 1000),
                                 'peak_amplitude': peak,
                                 'peak_change_to_spherical_db': float(20 * np.log10(peak / ref)),
                                 'local_power_HPBW_deg': _field_width(angle, profile),
                                 'local_range_power_FWHM_mm': _field_width(rg * 1000, image[:, ia]),
                                 'outside_8deg_peak_db': float(20 * np.log10(abs(profile[outside]).max() / peak))})
                    if distance == 5:
                        arrays[f'{method}_{target:g}deg_{focus}'] = image
    arrays['angles_deg'] = angle
    arrays['range_m_for_5m'] = np.arange(4.88, 5.120001, .002)
    np.savez_compressed(out / 'nearfield_5m_raw_images.npz', **arrays)
    nm.write_csv(out / 'metrics.csv', rows)
    fig, axs = plt.subplots(1, 2, figsize=(9, 3.6))
    for ax, target in zip(axs, (0., 40.)):
        for method in METHODS:
            group = [r for r in rows if r['method'] == method and r['receive_focus'] == 'farfield' and r['target_angle_deg'] == target]
            _series(ax, [r['distance_m'] for r in group], [r['peak_change_to_spherical_db'] for r in group], method)
        ax.set(xlabel='True range (m)', ylabel='Plane/spherical Rx-focus peak change (dB)', title=f'Target angle {target:g} deg')
        ax.legend(fontsize=8)
    _save_figure(fig, out / 'nearfield_receive_focus_loss')
    nm.dump(out / 'scope.json', {'generator': 'Exact spherical Tx-target-Rx paths, common Tx and common Rx frontend.',
                                'focusing': 'Spherical vs plane-wave Rx focus; Tx remains fixed sector far-field TTD, not nearfield Tx focusing.',
                                'angle_grid_deg': .25, 'range_grid_m': .002,
                                'source_energy': 'Final total Tx element excitation energy1 for all methods',
                                'finite_domain': 'Only ranges5/15/30/60m, angles0/40deg; grid-aligned target estimates do not demonstrate arbitrarily fine localization.'})
    return {'rows': len(rows)}


def run_offgrid(out):
    out.mkdir(parents=True, exist_ok=True)
    cfg = nm.Config()
    rg = np.arange(30.483, 30.524001, .001)
    true_range = 30.5037
    targets = (-41.37, -17.23, .37, 24.13, 47.41)
    rows = []
    for method in METHODS:
        for target in targets:
            raw, matched, meta = nm.point_raw(cfg, method, [true_range], [target], geometry='farfield')
            for step in (1., .25, .1, .05):
                angle = np.arange(-50., 50.000001, step)
                image = nm.reconstruct_points(cfg, method, raw, matched, rg, angle)
                ir, ia = np.unravel_index(np.argmax(abs(image)), image.shape)
                rows.append({'method': method, 'angular_step_deg': step,
                             'truth_angle_deg': target, 'truth_range_m': true_range,
                             'estimate_angle_deg': float(angle[ia]), 'estimate_range_m': float(rg[ir]),
                             'absolute_angle_error_deg': float(abs(angle[ia] - target)),
                             'absolute_range_error_mm': float(abs(rg[ir] - true_range) * 1000),
                             'peak_amplitude': float(abs(image[ir, ia]))})
    nm.write_csv(out / 'metrics.csv', rows)
    fig, axs = plt.subplots(1, 2, figsize=(9, 3.7))
    for method in METHODS:
        for ax, metric, label in ((axs[0], 'absolute_angle_error_deg', 'Maximum absolute angle error (deg)'),
                                   (axs[1], 'absolute_range_error_mm', 'Maximum absolute range error (mm)')):
            steps = [1., .25, .1, .05]
            maximum = [max(r[metric] for r in rows if r['method'] == method and r['angular_step_deg'] == step) for step in steps]
            _series(ax, steps, maximum, method)
            ax.set(xlabel='Angular grid step (deg)', ylabel=label, xscale='log')
            ax.legend(fontsize=8)
    _save_figure(fig, out / 'offgrid_convergence')
    nm.dump(out / 'scope.json', {'generator': 'Continuous-coordinate farfield point echoes; full physical-frequency phase/TTD and common frontend.',
                                'truth_angles_deg': targets, 'truth_range_m': true_range,
                                'reconstruction_range_step_m': .001,
                                'not_large_scene_accuracy': 'This check does not assign its point precision to201-node grid scene reconstruction.'})
    return {'rows': len(rows)}


def run_point_scene(out):
    out.mkdir(parents=True, exist_ok=True)
    cfg = nm.Config()
    truth_r = np.repeat([30.2, 30.4, 30.6], 6)
    truth_a = np.tile([-40., -24., -8., 8., 24., 40.], 3)
    rng = np.random.default_rng(2026101001)
    coefficient = np.exp(1j * rng.uniform(-np.pi, np.pi, 18)) / np.sqrt(18)
    rg = np.arange(30., 30.800001, .0025)
    angle = np.arange(-50., 50.000001, .25)
    images, records = {}, []
    arrays = {'range_m': rg, 'angles_deg': angle, 'true_range_m': truth_r,
              'true_angle_deg': truth_a, 'complex_reflectivity': coefficient}
    for method in METHODS:
        raw, matched, meta = nm.point_raw(cfg, method, truth_r, truth_a, coefficient, geometry='spherical')
        image = nm.reconstruct_points(cfg, method, raw, matched, rg, angle, focus_model='farfield')
        images[method] = image
        arrays[f'{method}_complex_image'] = image
        nm.dump(out / f'{method}_configuration.json', {'config': asdict(cfg), 'meta': meta})
        support = np.zeros(image.shape, bool)
        guard = np.zeros(image.shape, bool)
        for rr, aa in zip(truth_r, truth_a):
            ir = int(np.argmin(abs(rg - rr))); ia = int(np.argmin(abs(angle - aa)))
            support[ir, ia] = True
            guard |= (abs(rg[:, None] - rr) <= .06 + 1e-12) & (abs(angle[None, :] - aa) <= 1 + 1e-12)
        mag = abs(image)
        background = ~guard
        peaks, vals = [], []
        # Local target samples are a labelled diagnostic, not blind Pd.
        for rr, aa in zip(truth_r, truth_a):
            region = (abs(rg[:, None] - rr) <= .02 + 1e-12) & (abs(angle[None, :] - aa) <= 1 + 1e-12)
            vals.append(float(mag[region].max()))
        records.append({'method': method, 'source_total_reflectivity_power': float(np.sum(abs(coefficient) ** 2)),
                        'actual_Tx_element_energy': meta['final_tx_element_excitation_energy'],
                        'global_image_peak_amplitude': float(mag.max()),
                        'minimum_local_truth_neighborhood_peak_amplitude': min(vals),
                        'maximum_local_truth_neighborhood_peak_amplitude': max(vals),
                        'support_leakage': float(np.sum(mag[background] ** 2) / np.sum(mag ** 2)),
                        'background_fraction': float(background.mean()),
                        'CNR_linear_magnitude': float((mag[support].mean() - mag[background].mean()) / mag[background].std()),
                        'dense_image_magnitude_correlation': None})
    for record in records:
        a = abs(images[record['method']]).ravel()
        b = abs(images['dense']).ravel()
        record['dense_image_magnitude_correlation'] = float(np.corrcoef(a, b)[0, 1])
    nm.write_csv(out / 'metrics.csv', records)
    nm.write_csv(out / 'continuous_truth.csv', [{'target': i + 1, 'range_m': float(rr), 'angle_deg': float(aa),
                                              'coefficient_real': float(coefficient[i].real), 'coefficient_imag': float(coefficient[i].imag)}
                                             for i, (rr, aa) in enumerate(zip(truth_r, truth_a))])
    np.savez_compressed(out / 'raw_18point_images_and_truth.npz', **arrays)
    for common in (False, True):
        fig, axs = plt.subplots(1, 3, figsize=(12, 3.9))
        ref = max(abs(im).max() for im in images.values()) if common else None
        for ax, method in zip(axs, METHODS):
            mag = abs(images[method]); reference = mag.max() if ref is None else ref
            dB = 20 * np.log10(np.maximum(mag, 1e-30) / reference)
            mesh = ax.pcolormesh(angle, rg, np.maximum(dB, -50), shading='auto',
                                 cmap='viridis', vmin=-50, vmax=0, rasterized=True)
            ax.scatter(truth_a, truth_r, facecolors='none', edgecolors='white', s=22,
                       marker='o', linewidths=.6, label='Continuous truth')
            ax.set(xlabel='Angle (deg)', ylabel='Range (m)', title=LABELS[method], xlim=(-50, 50), ylim=(30, 30.8))
            ax.legend(fontsize=7, loc='lower left')
            fig.colorbar(mesh, ax=ax, label='Amplitude (dB)')
        _save_figure(fig, out / ('18point_common_reference' if common else '18point_per_image_reference'))
    nm.dump(out / 'scope.json', {'scene': '18 continuous points,3 exact range layers, single geometry and phase seed2026101001, noiseless.',
                                'generator': 'Exact spherical Tx-target-Rx; far-field Rx focus, physical TTD, common Tx and Rx frontend.',
                                'all_sources': 'Actual total transmit-element excitation energy1; sum individual reflectivity powers1.',
                                'evaluation': 'Fixed30–30.8m, +/-50deg; strict coordinate±0.06m/±1deg protection mask.',
                                'local_truth_diagnostic': 'Neighborhood peak amplitudes use truth only as labelled amplitude diagnostics, not blind detection evidence.',
                                'display': 'Separate per-image and common reference; no TVG/gamma/smoothing in metrics.',
                                'storage': 'Complex128 raw reconstructed images and coefficients; no float32 metric calculations.'})
    return {'rows': len(records), 'target_count': 18, 'image_shape': list(images['dense'].shape)}


def _global_and_local_correlation(y_fft, matched_fft, cfg, os_corr=8):
    corr = np.fft.fftshift(np.fft.ifft(y_fft * np.conj(matched_fft)))
    fine = signal.resample(corr, len(corr) * os_corr)
    lag_s = np.arange(len(fine)) / (cfg.fs * os_corr) - len(corr) / (2 * cfg.fs)
    global_i = int(np.argmax(abs(fine)))
    local_indices = np.flatnonzero(abs(lag_s) <= .0006 + 1e-15)
    local_i = int(local_indices[np.argmax(abs(fine[local_indices]))])
    return lag_s, fine, global_i, local_i


def _continuous_tx_echo(tx, beta, fd, high_rate, cfg, nfft, interpolated=None):
    """Explicit continuous-time query, NOT rounded pulse-length resampling.

    Zero-padded Fourier interpolation represents the same finite digital pulse.
    beta*t is queried on its high-rate interpolant with linear fractional index
    interpolation. Exact transmitted finite support[0,len(tx)/fs) is imposed.
    A fixed8192 sample receive record is shared by all velocities. Pulse samples
    are not stretched by choosing an integer output length.
    """
    padded = np.pad(tx, (512, nfft - 512 - len(tx)))
    high = signal.resample(padded, nfft * high_rate) if interpolated is None else interpolated
    receive_t = np.arange(nfft) / cfg.fs
    query = beta * np.arange(nfft)
    inside = (query >= 0) & (query < len(tx))
    index = (512 + query[inside]) * high_rate
    low = np.floor(index).astype(int); frac = index - low
    y = np.zeros(nfft, dtype=complex)
    y[inside] = high[low] * (1 - frac) + high[low + 1] * frac
    y *= np.exp(2j * np.pi * fd * receive_t)
    energy = float(np.sum(abs(y) ** 2))
    return y / np.sqrt(energy), energy


def run_doppler(out):
    out.mkdir(parents=True, exist_ok=True)
    cfg = nm.Config()
    nfft = 8192
    velocities = (-10., -5., -2., -1., 0., 1., 2., 5., 10.)
    interpolation_rates = (8, 16, 32)
    records, arrays, summary, convergence = [], {}, [], []
    g = nm.frontend_response(cfg, nfft)
    for method in METHODS:
        tx, meta = nm.waveforms(cfg, method)
        send, matched, nu, smeta = nm.spectra(cfg, method, nfft)
        for k in range(tx.shape[1]):
            template = matched[:, k]
            nominal = np.fft.fft(tx[:, k], nfft) * g
            lag, zero_corr, izero, _ = _global_and_local_correlation(nominal, template, cfg)
            ref = float(abs(zero_corr[izero]))
            for rate in interpolation_rates:
                padded = np.pad(tx[:, k], (512, nfft - 512 - len(tx)))
                interpolated = signal.resample(padded, nfft * rate)
                for velocity in velocities:
                    beta = (cfg.c - velocity) / (cfg.c + velocity)
                    fd = (beta - 1) * cfg.fc
                    echo, pre_frontend_e = _continuous_tx_echo(tx[:, k], beta, fd, rate, cfg, nfft, interpolated)
                    received_fft = np.fft.fft(echo) * g
                    lag, corr, ig, il = _global_and_local_correlation(received_fft, template, cfg)
                    record = {'method': method, 'channel': k + 1, 'radial_velocity_mps': velocity,
                              'model': 'continuous_time_query_and_carrier_shift', 'interpolation_rate': rate,
                              'time_scale_beta': beta, 'carrier_shift_hz': fd,
                              'receive_record_samples': nfft,
                              'received_pulse_support_duration_ms': cfg.duration / beta * 1000,
                              'pre_frontend_energy_before_pulse_normalization': pre_frontend_e,
                              'Rx_frontend_energy_after_pulse_normalization': float(np.sum(abs(received_fft) ** 2) / nfft),
                              'global_peak_change_db': float(20 * np.log10(abs(corr[ig]) / ref)),
                              'global_peak_lag_ms': float(lag[ig] * 1000),
                              'global_peak_equivalent_range_bias_mm': float(lag[ig] * cfg.c / 2 * 1000),
                              'local_0p6ms_peak_change_db': float(20 * np.log10(abs(corr[il]) / ref)),
                              'local_0p6ms_peak_lag_ms': float(lag[il] * 1000),
                              'local_0p6ms_equivalent_range_bias_mm': float(lag[il] * cfg.c / 2 * 1000)}
                    records.append(record)
                    if rate == 32:
                        arrays[f'{method}_ch{k + 1}_v{velocity:g}_continuous32_correlation'] = corr
                del interpolated
            for velocity in velocities:
                beta = (cfg.c - velocity) / (cfg.c + velocity)
                fd = (beta - 1) * cfg.fc
                echo = np.pad(tx[:, k], (0, nfft - len(tx))) * np.exp(2j * np.pi * fd * np.arange(nfft) / cfg.fs)
                echo /= np.sqrt(np.sum(abs(echo) ** 2))
                received_fft = np.fft.fft(echo) * g
                lag, corr, ig, il = _global_and_local_correlation(received_fft, template, cfg)
                records.append({'method': method, 'channel': k + 1, 'radial_velocity_mps': velocity,
                                'model': 'frequency_shift_only', 'interpolation_rate': 0,
                                'time_scale_beta': beta, 'carrier_shift_hz': fd,
                                'receive_record_samples': nfft,
                                'received_pulse_support_duration_ms': cfg.duration * 1000,
                                'pre_frontend_energy_before_pulse_normalization': 1.,
                                'Rx_frontend_energy_after_pulse_normalization': float(np.sum(abs(received_fft) ** 2) / nfft),
                                'global_peak_change_db': float(20 * np.log10(abs(corr[ig]) / ref)),
                                'global_peak_lag_ms': float(lag[ig] * 1000),
                                'global_peak_equivalent_range_bias_mm': float(lag[ig] * cfg.c / 2 * 1000),
                                'local_0p6ms_peak_change_db': float(20 * np.log10(abs(corr[il]) / ref)),
                                'local_0p6ms_peak_lag_ms': float(lag[il] * 1000),
                                'local_0p6ms_equivalent_range_bias_mm': float(lag[il] * cfg.c / 2 * 1000)})
                arrays[f'{method}_ch{k + 1}_v{velocity:g}_frequency_only_correlation'] = corr
    for method in METHODS:
        for velocity in velocities:
            for model in ('continuous_time_query_and_carrier_shift', 'frequency_shift_only'):
                group = [r for r in records if r['method'] == method and r['radial_velocity_mps'] == velocity
                         and r['model'] == model and r['interpolation_rate'] in (0, 32)]
                worst_loss = min(group, key=lambda x: x['global_peak_change_db'])
                worst_bias = max(group, key=lambda x: abs(x['global_peak_equivalent_range_bias_mm']))
                worst_local = min(group, key=lambda x: x['local_0p6ms_peak_change_db'])
                summary.append({'method': method, 'radial_velocity_mps': velocity, 'model': model,
                                'worst_global_peak_change_db': worst_loss['global_peak_change_db'],
                                'worst_peak_change_channel': worst_loss['channel'],
                                'maximum_absolute_global_range_bias_mm': abs(worst_bias['global_peak_equivalent_range_bias_mm']),
                                'worst_range_bias_channel': worst_bias['channel'],
                                'worst_local_0p6ms_peak_change_db': worst_local['local_0p6ms_peak_change_db']})
            for rate in (8, 16):
                errors = []
                for k in range(nm.waveforms(cfg, method)[0].shape[1]):
                    q = [r for r in records if r['method'] == method and r['channel'] == k + 1 and
                         r['radial_velocity_mps'] == velocity and r['model'] == 'continuous_time_query_and_carrier_shift']
                    a = next(r for r in q if r['interpolation_rate'] == rate)
                    b = next(r for r in q if r['interpolation_rate'] == 32)
                    errors.append((abs(a['global_peak_change_db'] - b['global_peak_change_db']),
                                   abs(a['global_peak_equivalent_range_bias_mm'] - b['global_peak_equivalent_range_bias_mm']),
                                   abs(a['local_0p6ms_peak_change_db'] - b['local_0p6ms_peak_change_db'])))
                loss = max(x[0] for x in errors); bias = max(x[1] for x in errors)
                convergence.append({'method': method, 'radial_velocity_mps': velocity,
                                    'tested_interpolation_rate': rate, 'reference_interpolation_rate': 32,
                                    'max_global_peak_change_difference_db': loss,
                                    'max_global_range_bias_difference_mm': bias,
                                    'max_local_peak_change_difference_db': max(x[2] for x in errors),
                                    'fixed_loss_tolerance_db': .002,
                                    'fixed_range_tolerance_mm': .5,
                                    'convergence_under_fixed_criteria': loss <= .002 and bias <= .5})
    arrays['lag_s'] = lag
    np.savez_compressed(out / 'all_global_correlations.npz', **arrays)
    nm.write_csv(out / 'all_channel_metrics.csv', records)
    nm.write_csv(out / 'worst_channel_summary.csv', summary)
    nm.write_csv(out / 'continuous_interpolation_convergence.csv', convergence)
    fig, axs = plt.subplots(2, 2, figsize=(10, 7))
    for method in METHODS:
        for col, model in enumerate(('frequency_shift_only', 'continuous_time_query_and_carrier_shift')):
            group = [r for r in summary if r['method'] == method and r['model'] == model]
            _series(axs[0, col], velocities, [r['worst_global_peak_change_db'] for r in group], method)
            _series(axs[1, col], velocities, [r['maximum_absolute_global_range_bias_mm'] for r in group], method)
    for col, title in enumerate(('Frequency shift only', 'Continuous time query + carrier shift')):
        axs[0, col].set(xlabel='Radial velocity (m/s; positive receding)', ylabel='Worst global peak change (dB)', title=title)
        axs[1, col].set(xlabel='Radial velocity (m/s; positive receding)', ylabel='Maximum absolute equivalent range bias (mm)')
        for row in range(2):
            axs[row, col].legend(fontsize=8)
    _save_figure(fig, out / 'continuous_doppler_global_peak_and_bias')
    nm.dump(out / 'scope.json', {'velocity': 'positive receding; beta=(c-v)/(c+v), carrier shift=(beta-1)fc. FDM baseband offsets are time-scaled automatically.',
                                'continuous_query': 'Explicit u(beta*t) on zero-padded high-rate Fourier interpolants with fractional-index linear sampling; transmitted finite support fixed; common8192 receive record, no round(n/beta) construction.',
                                'ordering': 'Time scale actual Tx pulse, apply carrier shift and unit received-pulse energy normalization, then common Rx frontend, then effective nominal unit-energy matched filter.',
                                'interpolation_controls': '8x/16x/32x; canonical32x result; fixed0.002dB/0.5mm numerical convergence checks, unfavorable/failing checks retained.',
                                'correlation': 'Full padded8192 circular FFT correlation equals finite linear lag record within its support;8x Fourier correlation interpolation. Global maximum from all lags. Fixedabs(lag)<=0.6ms is a separate local diagnostic, never substitutes for global result.',
                                'worst_channels': 'Most negative peak change and largest absolute lag are independent maxima and may belong to different channels.',
                                'limitation': 'Uncompensated waveform mismatch only; not motion-compensated3D image reconstruction or measured platform data; continuous signal reconstruction is a numerical interpolation model.'})
    return {'all_channel_rows': len(records), 'summary_rows': len(summary), 'convergence_rows': len(convergence),
            'convergence_pass_rows': sum(r['convergence_under_fixed_criteria'] for r in convergence)}


RUNNERS = {'waveforms': run_waveforms, 'range': run_range, 'nearfield': run_nearfield,
           'offgrid': run_offgrid, 'point_scene': run_point_scene, 'doppler': run_doppler}


def run_all(selected=None):
    selected = list(RUNNERS) if selected is None else list(selected)
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now(timezone(timedelta(hours=8))).strftime('%Y%m%d_%H%M%S_%f')
    out = OUT_ROOT / f'run_{stamp}'
    out.mkdir()
    snapshot = out / 'code_snapshot'; snapshot.mkdir()
    for path in (Path(__file__), Path(nm.__file__)):
        shutil.copyfile(path, snapshot / path.name)
    hashes = {path.name: nm.sha(path) for path in snapshot.iterdir()}
    nm.dump(out / 'protocol_before_run.json', {'configuration': asdict(nm.Config()), 'selected': selected,
                                              'code_sha256': hashes, 'frozen_axes': {
                                                  'Tx_control': ['primary2000samples', 'full2128samples'],
                                                  'double_target_separation_mm': [5, 10, 20, 40, 80, 120],
                                                  'double_target_phase_rad': .7, 'double_target_angle_deg': 8.,
                                                  'nearfield_range_m': [5, 15, 30, 60], 'nearfield_angle_deg': [0, 40],
                                                  'offgrid_steps_deg': [1, .25, .1, .05],
                                                  'Doppler_velocity_m_s': [-10, -5, -2, -1, 0, 1, 2, 5, 10],
                                                  'continuous_interpolation_rates': [8, 16, 32],
                                                  'numerical_convergence_db': .002, 'numerical_convergence_mm': .5},
                                              'no_favorable_result_gate': True,
                                              'author_validation': 'pending', 'originals_untouched': True})
    started = time.perf_counter()
    results = {}
    with (out / 'execution.log').open('w', encoding='utf-8') as log:
        class Tee:
            def write(self, value):
                sys.__stdout__.write(value); sys.__stdout__.flush(); log.write(value); log.flush()
            def flush(self):
                sys.__stdout__.flush(); log.flush()
        with contextlib.redirect_stdout(Tee()), contextlib.redirect_stderr(Tee()):
            for name in selected:
                sub = out / name; begin = time.perf_counter()
                print(f'START {name}', flush=True)
                try:
                    data = RUNNERS[name](sub)
                    data.update({'status': 'computed_pending_author_validation', 'elapsed_seconds': time.perf_counter() - begin})
                    results[name] = data
                    nm.dump(sub / 'receipt.json', data)
                    print(f'FINISH {name}: {data}', flush=True)
                except Exception as exc:
                    traceback.print_exc()
                    data = {'status': 'failed_diagnostic_preserved', 'error': repr(exc), 'elapsed_seconds': time.perf_counter() - begin}
                    results[name] = data
                    nm.dump(sub / 'receipt.json', data)
                    print(f'FAIL {name}: {data}', flush=True)
    artifacts = [{'path': str(p.relative_to(out)), 'size_bytes': p.stat().st_size, 'sha256': nm.sha(p)}
                 for p in out.rglob('*') if p.is_file()]
    receipt = {'run_path': str(out.relative_to(ROOT)), 'selected': selected, 'results': results,
               'elapsed_seconds': time.perf_counter() - started,
               'all_requested_computed': all(q['status'] == 'computed_pending_author_validation' for q in results.values()),
               'code_hashes': hashes, 'artifacts': artifacts,
               'source_model_unchanged_during_run': nm.sha(Path(nm.__file__)) == hashes['numerical_model.py'],
               'no_manuscript_edit': True, 'author_validation': 'pending'}
    nm.dump(out / 'receipt.json', receipt)
    nm.dump(OUT_ROOT / 'latest.json', {'path': str(out.relative_to(ROOT)), 'all_requested_computed': receipt['all_requested_computed'], 'selected': selected})
    print(json.dumps({'run': str(out), 'results': results}, ensure_ascii=False), flush=True)
    return receipt


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--only', nargs='+', choices=list(RUNNERS))
    args = parser.parse_args()
    result = run_all(args.only)
    raise SystemExit(0 if result['all_requested_computed'] else 1)
