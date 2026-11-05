"""Re-export Figure 5 from frozen plotted data; change panel (c) encoding only.

Run: python -X utf8 plot_figure5_legend.py
No simulation imports, new observations, or numeric metric calculations occur.
"""
from pathlib import Path
import csv
import hashlib
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
HERE = ROOT / 'outputs/figures/figure05'
DATA = ROOT / 'data/waveforms/figure05'
COLORS = {'cdm': '#D55E00', 'fdm': '#009E73'}
plt.rcParams.update({
    'font.family': 'serif', 'font.size': 10, 'axes.labelsize': 10,
    'axes.titlesize': 10, 'xtick.labelsize': 9, 'ytick.labelsize': 9,
    'legend.fontsize': 9, 'pdf.fonttype': 42, 'axes.spines.top': False,
    'axes.spines.right': False, 'savefig.dpi': 600,
})
with (DATA / 'finite_ambiguity_worst_pairs.csv').open(encoding='utf-8-sig', newline='') as file:
    SUMMARY = list(csv.DictReader(file))
CFG = json.loads((DATA / 'configuration.json').read_text(encoding='utf-8'))
RAW = np.load(DATA / 'figure5_actual_arrays.npz')
LAG = RAW['lag_s']
GEOMETRY = json.loads((DATA / 'original_geometry.json').read_text(encoding='utf-8'))


def sha(path):
    with Path(path).open('rb') as file:
        return hashlib.file_digest(file, 'sha256').hexdigest()


def build(repair: bool):
    fig, axes = plt.subplots(2, 2, figsize=(6, 4.9))
    curves = {}
    for ax, key in [(axes[0, 0], 'full_correlation'), (axes[0, 1], 'local_correlation')]:
        power = 20 * np.log10(np.maximum(abs(RAW[key]), 1e-30))
        ax.plot(LAG * 1000, power, color=COLORS['fdm'], lw=1.25)
        curves[key] = np.stack((LAG * 1000, power))
    axes[0, 0].axvspan(-.6, .6, color='.9')
    axes[0, 0].set(xlim=(-5.35, 5.35), ylim=(-110, 2), xlabel='Delay (ms)',
                   ylabel='Normalized cross power (dB)', title='(a) Full effective-pulse support')
    axes[0, 1].set(xlim=(-.6, .6), ylim=(-105, -65), xlabel='Delay (ms)',
                   ylabel='Normalized cross power (dB)', title='(b) Evaluated local delay interval')
    for method in ('cdm', 'fdm'):
        block = [row for row in SUMMARY if row['method'] == method
                 and row['pair_type'] == 'distinct_cross']
        for scope, linestyle in [('global', '-'), ('local_0p6ms', '--')]:
            velocity = np.array([float(row['velocity_mps_positive_receding']) for row in block])
            values = np.array([float(row['worst_' + scope + '_peak_db']) for row in block])
            axes[1, 0].plot(velocity, values, color=COLORS[method], ls=linestyle,
                            marker='s' if method == 'cdm' else '^', markersize=3, lw=1.2,
                            label=method.upper() + (' full' if scope == 'global' else ' local'))
            curves[method + '_' + scope] = np.stack((velocity, values))
    axes[1, 0].set(xlabel='Radial velocity (m/s; receding +)',
                   ylabel='Worst distinct-pair peak (dB)', ylim=(-100, 0),
                   title='(c) Continuous Doppler model')
    axes[1, 0].legend(fontsize=9, ncol=2, loc='center', bbox_to_anchor=(.5, .45),
                      framealpha=.9, handlelength=1.4, columnspacing=.8)
    transmit = RAW['actual_fdm_tx']
    frequencies = CFG['fc'] + np.fft.fftshift(np.fft.fftfreq(65536, 1 / CFG['fs']))
    styles = ['-', '--', ':', '-.', '-', '--']
    palette = ['#0072B2', '#D55E00', '#009E73', '#CC79A7', '#E69F00', '#56B4E9']
    for channel in range(6):
        power = abs(np.fft.fftshift(np.fft.fft(transmit[:, channel], 65536))) ** 2
        spectrum = 10 * np.log10(np.maximum(power / power.max(), 1e-12))
        axes[1, 1].plot(frequencies / 1000, spectrum, color=palette[channel],
                        ls=styles[channel], lw=1.15, label=f"{CFG['centres'][channel]/1000:g} kHz")
        curves[f'spectrum_{channel+1}'] = np.stack((frequencies / 1000, spectrum))
    axes[1, 1].set(xlim=(395, 505), ylim=(-60, 2), xlabel='RF frequency (kHz)',
                   ylabel='Normalized Tx power (dB)', title='(d) Actual six-band Tx spectra')
    axes[1, 1].legend(fontsize=9, ncol=2, loc='lower center', framealpha=.85)
    for ax in axes.flat:
        ax.grid(alpha=.15)
    # Freeze the same panel placement as the original before repairing the legend.
    fig.tight_layout(pad=.7, h_pad=1.1, w_pad=.9)
    # Matplotlib tight_layout depends slightly on installed renderer/font versions.
    # Restore the axes rectangles read directly from the original vector PDF,
    # rather than allowing a style-only repair to shift the other three panels.
    page_height = GEOMETRY['pdf_page_rect_pt'][3]
    for ax, rect in zip(axes.flat, GEOMETRY['panel_rectangles_pt']):
        left, top, right, bottom = rect
        ax.set_position(((left + GEOMETRY['crop_left_pt']) / 432,
                         (page_height - bottom + GEOMETRY['crop_bottom_pt']) / 352.8,
                         (right - left) / 432, (bottom - top) / 352.8))
    if repair:
        for line in axes[1, 0].lines:
            local = line.get_label().endswith(' local')
            if local:
                line.set_linestyle((0, (2.5, 1.8)))
                line.set_markerfacecolor('white')
                line.set_markeredgewidth(.9)
            else:
                line.set_markerfacecolor(line.get_color())
        # A longer sample visibly contains several on/off intervals. Local markers
        # are hollow, so the encoding remains legible when printed in grayscale.
        axes[1, 0].legend(fontsize=9, ncol=2, loc='center', bbox_to_anchor=(.5, .48),
                          framealpha=.9, handlelength=2.3, handletextpad=.3,
                          columnspacing=.4, numpoints=1)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    legend_box = axes[1, 0].get_legend().get_window_extent(renderer)
    axes_box = axes[1, 0].get_window_extent(renderer)
    return fig, axes, curves, {
        'axes_positions': [list(ax.get_position().bounds) for ax in axes.flat],
        'axes_limits': [{'xlim': list(ax.get_xlim()), 'ylim': list(ax.get_ylim())} for ax in axes.flat],
        'legend_width_px': float(legend_box.width), 'panel_c_width_px': float(axes_box.width),
        'legend_inside_panel_c': bool(legend_box.x0 >= axes_box.x0 and legend_box.x1 <= axes_box.x1),
    }


def export_main():
    HERE.mkdir(parents=True, exist_ok=True)
    original, original_axes, original_curves, before = build(False)
    original.savefig(HERE / 'reproduced_original_style.pdf', bbox_inches='tight')
    original.savefig(HERE / 'reproduced_original_style.png', bbox_inches='tight')
    plt.close(original)
    repaired, repaired_axes, repaired_curves, after = build(True)
    for key in original_curves:
        if not np.array_equal(original_curves[key], repaired_curves[key]):
            raise AssertionError(f'Curve changed: {key}')
    if before['axes_positions'] != after['axes_positions'] or before['axes_limits'] != after['axes_limits']:
        raise AssertionError('Panel positions/ranges changed')
    repaired.savefig(HERE / 'figure_fdm_ambiguity_scope.pdf', bbox_inches='tight')
    repaired.savefig(HERE / 'figure_fdm_ambiguity_scope.png', bbox_inches='tight')
    repaired.savefig(HERE / 'figure_fdm_ambiguity_scope.svg', bbox_inches='tight')
    panel_bbox = repaired_axes[1, 0].get_tightbbox(repaired.canvas.get_renderer()).transformed(repaired.dpi_scale_trans.inverted())
    repaired.savefig(HERE / 'panel_c.png', bbox_inches=panel_bbox.expanded(1.03, 1.03), dpi=600)
    plt.close(repaired)
    np.savez_compressed(HERE / 'figure5_plotted_curves.npz', **original_curves)
    receipt = {
        'scope': 'Panel (c) line/marker/legend style only; no experiments rerun.',
        'raw_curve_arrays_bitwise_unchanged': True,
        'panel_positions_unchanged': before['axes_positions'] == after['axes_positions'],
        'axis_limits_unchanged': before['axes_limits'] == after['axes_limits'],
        'before': before, 'after': after,
        'encoding': {'method': 'CDM square vs FDM triangle plus original color',
                     'delay_scope': 'full solid/filled vs local short-dashed/hollow',
                     'legend_handlelength': 2.3, 'local_dash_sequence_pt': [3.0, 2.16],
                     'source_font_minimum_pt': 9.0},
        'outputs_sha256': {p.name: sha(p) for p in sorted(HERE.glob('figure_fdm_ambiguity_scope.*'))},
        'snapshot_sha256': {p.name: sha(p) for p in sorted(DATA.iterdir())},
    }
    (HERE / 'FIGURE5_LEGEND_QA.json').write_text(json.dumps(receipt, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(receipt, ensure_ascii=False))


if __name__ == '__main__':
    export_main()
