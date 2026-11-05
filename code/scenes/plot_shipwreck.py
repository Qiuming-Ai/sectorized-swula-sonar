"""Redraw v08 manuscript Figure 13 from the released shipwreck-data overlay.

Extract the overlay at the package root before running. This uses only the
saved common-noise image arrays and the display-only geometric reference.
It does not read the STL or arrival catalogue or rerun the sonar simulation.
"""
import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat

SCRIPT = Path(__file__).resolve()
ROOT = SCRIPT.parents[2]
GEOMETRY = ROOT / 'data/shipwreck/geometry'
VISUAL = GEOMETRY / 'A11_plot_data.mat'
IMAGES = GEOMETRY / 'noisy_common_per_element_noise_10dB.npz'
STEM = 'shipwreck_common_per_element_noise_common_scale'
METHODS = ('dense', 'cdm', 'fdm')
LABELS = ('Dense ULA', 'CDM SW-ULA', 'FDM SW-ULA')


def sha256(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(4*1024*1024), b''):
            digest.update(block)
    return digest.hexdigest()


def destination(out, name):
    path = out / name
    if path.resolve().parent != out:
        raise RuntimeError(f'Refusing redirected figure output: {path}')
    return path


def level(image, scale=None):
    # Match the archived Figure 13 renderer, including the pre-colormap floor.
    magnitude = np.abs(image)
    reference = float(magnitude.max()) if scale is None else scale
    if not np.isfinite(reference) or reference <= 0:
        raise ValueError('The saved display reference must be finite and positive')
    return np.maximum(20*np.log10(np.maximum(magnitude, 1e-30)/reference), -50)


def render(out):
    visual = loadmat(VISUAL, squeeze_me=True,
                     variable_names=['truth_magnitude', 'range_m', 'angle_deg'])
    with np.load(IMAGES, allow_pickle=False) as data:
        r, a = data['ranges_m'], data['angles_deg']
        images = [data[method] for method in METHODS]
    scale = max(np.abs(image).max() for image in images)
    panels = [(visual['truth_magnitude'], visual['range_m'], visual['angle_deg'],
               'Reference (display only)', None)]
    panels.extend((image, r, a, label, scale) for image, label in zip(images, LABELS))

    plt.rcParams.update({
        'font.family': 'serif', 'font.serif': ['Times New Roman', 'DejaVu Serif'],
        'font.size': 10, 'axes.labelsize': 10, 'xtick.labelsize': 9.6,
        'ytick.labelsize': 9.6, 'legend.fontsize': 9.6, 'pdf.fonttype': 42,
        'ps.fonttype': 42, 'axes.spines.top': False, 'axes.spines.right': False,
    })
    fig, axes = plt.subplots(2, 2, figsize=(7.15, 4.15), sharex=True, sharey=True)
    try:
        for idx, (ax, (image, rr, aa, label, ref)) in enumerate(zip(axes.flat, panels)):
            if image.shape != (len(rr), len(aa)) or not np.isfinite(image).all():
                raise ValueError(f'Invalid saved range-angle image for {label}')
            # Use the actual angle coordinates: no resampling or deskewing.
            xx = rr[:, None]*np.sin(np.deg2rad(aa))[None, :]
            rp = rr[:, None]*np.cos(np.deg2rad(aa))[None, :]
            im = ax.pcolormesh(xx, rp, level(image, ref), shading='auto', cmap='cividis',
                               vmin=-40, vmax=0, rasterized=True)
            ax.text(.02, .96, f'({chr(97+idx)}) {label}', transform=ax.transAxes,
                    va='top', color='white', fontsize=9.6,
                    bbox=dict(facecolor='black', alpha=.4, edgecolor='none', pad=1))
            ax.set_aspect('equal', adjustable='box')
            ax.set_xlim(-95, 95)
            ax.set_ylim(26, 122)
            if idx//2:
                ax.set_xlabel('Lateral coordinate x (m)')
            if idx%2 == 0:
                ax.set_ylabel(r'Radial coordinate $R_\perp$ (m)')
        fig.subplots_adjust(left=.08, right=.89, bottom=.12, top=.98, wspace=.08, hspace=.12)
        cb = fig.add_axes([.91, .18, .017, .66])
        fig.colorbar(im, cax=cb, label='Display level (dB)')
        outputs = [destination(out, STEM+'.'+ext) for ext in ('pdf', 'png', 'svg')]
        for path in outputs:
            fig.savefig(path, dpi=600, bbox_inches='tight')
    finally:
        plt.close(fig)
    return outputs, float(scale)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path,
                        help='Output directory (default: outputs/figures in the package).')
    requested = parser.parse_args().out
    default = ROOT/'outputs/figures'
    out = (requested if requested is not None else default).resolve()
    if requested is None and out != default:
        parser.error('Refusing redirected default output directory.')
    if out.is_relative_to(ROOT) and not out.is_relative_to(ROOT/'outputs'):
        parser.error('Outputs inside the package must stay under outputs/.')
    inputs = (VISUAL, IMAGES)
    missing = [path.relative_to(ROOT).as_posix() for path in inputs if not path.is_file()]
    if missing:
        parser.error('Extract the shipwreck-data overlay at the package root. Missing: '+', '.join(missing))
    hashes = {path.relative_to(ROOT).as_posix(): sha256(path) for path in inputs}
    out.mkdir(parents=True, exist_ok=True)
    outputs, scale = render(out)
    if any(sha256(ROOT/name) != digest for name, digest in hashes.items()):
        raise RuntimeError('A released shipwreck input changed while plotting')
    receipt = {
        'figure': 13, 'source': 'v08 manuscript common-scale shipwreck panel',
        'input_sha256': hashes, 'inputs_unchanged': True, 'script_sha256': sha256(SCRIPT),
        'output_sha256': {path.name: sha256(path) for path in outputs},
        'nominal_Dense_raw_input_snr_db': 10, 'noise_normalization': 'common_per_element_noise',
        'coordinate_formula': 'x = r sin(theta), R_perp = r cos(theta)',
        'display_limits_db': [-40, 0], 'common_numerical_peak_amplitude': scale,
        'geometric_reference': 'display only, separately normalized',
        'no_simulation_or_metric_recalculation': True,
    }
    destination(out, 'figure13_replot_receipt.json').write_text(
        json.dumps(receipt, indent=2)+'\n', encoding='utf-8')
    print(f'Replotted Figure 13 into {out}')


if __name__ == '__main__':
    main()
