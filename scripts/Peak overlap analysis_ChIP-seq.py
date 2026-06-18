#!/usr/bin/env python3
"""
Peak overlap analysis for IPMC ChIP-seq experiment
Conditions: Control, ZL0454, TGFβ, TGFβ+ZL0454
Replicates: A (IPMC1-20), B (IPMC21-40)
Marks: H3K4me1, H3K4me3, H3K27ac, H3K122ac
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from matplotlib_venn import venn2, venn2_circles
import os

# ── Sample map ────────────────────────────────────────────────────────────────
# Format: IPMC_ID -> (replicate, condition, mark)
SAMPLE_MAP = {
    # Replicate A
    'IPMC1':  ('A', 'Control',       'input'),
    'IPMC2':  ('A', 'Control',       'H3K4me1'),
    'IPMC3':  ('A', 'Control',       'H3K4me3'),
    'IPMC4':  ('A', 'Control',       'H3K27ac'),
    'IPMC5':  ('A', 'Control',       'H3K122ac'),
    'IPMC6':  ('A', 'ZL0454',        'input'),
    'IPMC7':  ('A', 'ZL0454',        'H3K4me1'),
    'IPMC8':  ('A', 'ZL0454',        'H3K4me3'),
    'IPMC9':  ('A', 'ZL0454',        'H3K27ac'),
    'IPMC10': ('A', 'ZL0454',        'H3K122ac'),
    'IPMC11': ('A', 'TGFb',          'input'),
    'IPMC12': ('A', 'TGFb',          'H3K4me1'),
    'IPMC13': ('A', 'TGFb',          'H3K4me3'),
    'IPMC14': ('A', 'TGFb',          'H3K27ac'),
    'IPMC15': ('A', 'TGFb',          'H3K122ac'),
    'IPMC16': ('A', 'TGFb+ZL0454',   'input'),
    'IPMC17': ('A', 'TGFb+ZL0454',   'H3K4me1'),
    'IPMC18': ('A', 'TGFb+ZL0454',   'H3K4me3'),
    'IPMC19': ('A', 'TGFb+ZL0454',   'H3K27ac'),
    'IPMC20': ('A', 'TGFb+ZL0454',   'H3K122ac'),
    # Replicate B
    'IPMC21': ('B', 'Control',       'input'),
    'IPMC22': ('B', 'Control',       'H3K4me3'),
    'IPMC23': ('B', 'Control',       'H3K4me1'),
    'IPMC24': ('B', 'Control',       'H3K27ac'),
    'IPMC25': ('B', 'Control',       'H3K122ac'),
    'IPMC26': ('B', 'ZL0454',        'input'),
    'IPMC27': ('B', 'ZL0454',        'H3K4me3'),
    'IPMC28': ('B', 'ZL0454',        'H3K4me1'),
    'IPMC29': ('B', 'ZL0454',        'H3K27ac'),
    'IPMC30': ('B', 'ZL0454',        'H3K122ac'),
    'IPMC31': ('B', 'TGFb',          'input'),
    'IPMC32': ('B', 'TGFb',          'H3K4me3'),
    'IPMC33': ('B', 'TGFb',          'H3K4me1'),
    'IPMC34': ('B', 'TGFb',          'H3K27ac'),
    'IPMC35': ('B', 'TGFb',          'H3K122ac'),   # 0 reads – will be empty
    'IPMC36': ('B', 'TGFb+ZL0454',   'input'),
    'IPMC37': ('B', 'TGFb+ZL0454',   'H3K4me3'),
    'IPMC38': ('B', 'TGFb+ZL0454',   'H3K4me1'),
    'IPMC39': ('B', 'TGFb+ZL0454',   'H3K27ac'),
    'IPMC40': ('B', 'TGFb+ZL0454',   'H3K122ac'),
}

MARKS   = ['H3K4me1', 'H3K4me3', 'H3K27ac', 'H3K122ac']
REPS    = ['A', 'B']

# Comparisons of interest
COMPARISONS = [
    ('Control',     'ZL0454',      'Effect of ZL0454 (no TGFβ)'),
    ('TGFb',        'TGFb+ZL0454', 'Effect of ZL0454 on TGFβ response'),
    ('Control',     'TGFb',        'Effect of TGFβ alone'),
    ('ZL0454',      'TGFb+ZL0454', 'TGFβ effect in ZL0454-treated cells'),
]

PEAK_DIR = '.'   # ← adjust if broadPeak files live elsewhere

# ── Helpers ───────────────────────────────────────────────────────────────────

def broadpeak_path(ipmc_id):
    return os.path.join(PEAK_DIR, f'{ipmc_id}_broad_peaks.broadPeak')

def load_broadpeak(ipmc_id):
    fp = broadpeak_path(ipmc_id)
    try:
        df = pd.read_csv(fp, sep='\t', header=None,
                         names=['chr','start','end','name','score',
                                'strand','signalValue','pValue','qValue'])
        if df.empty:
            return None
        df['peak_id'] = df['chr'] + ':' + (df['start'] // 500).astype(str)
        return df
    except (FileNotFoundError, pd.errors.EmptyDataError):
        return None

def get_ipmc(rep, condition, mark):
    for iid, (r, c, m) in SAMPLE_MAP.items():
        if r == rep and c == condition and m == mark:
            return iid
    return None

def analyze_overlap(set_a, set_b, label_a, label_b, mark, rep, comp_label):
    overlap  = set_a & set_b
    a_only   = set_a - set_b
    b_only   = set_b - set_a
    n_a, n_b, n_ov = len(set_a), len(set_b), len(overlap)
    return {
        'comparison': comp_label,
        'replicate':  rep,
        'mark':       mark,
        'label_a':    label_a,
        'label_b':    label_b,
        'total_a':    n_a,
        'total_b':    n_b,
        'maintained': n_ov,
        'lost':       len(a_only),   # in A but not B  (lost when going A→B)
        'gained':     len(b_only),   # in B but not A  (gained when going A→B)
        'pct_overlap_a': n_ov / n_a * 100 if n_a else 0,
        'pct_overlap_b': n_ov / n_b * 100 if n_b else 0,
    }

# ── Run all comparisons ───────────────────────────────────────────────────────

all_results = []

for cond_a, cond_b, comp_label in COMPARISONS:
    print(f"\n{'='*70}")
    print(f"Comparison: {comp_label}  ({cond_a}  vs  {cond_b})")
    print(f"{'='*70}")
    for rep in REPS:
        print(f"\n  Replicate {rep}:")
        for mark in MARKS:
            id_a = get_ipmc(rep, cond_a, mark)
            id_b = get_ipmc(rep, cond_b, mark)
            if id_a is None or id_b is None:
                print(f"    {mark}: sample not found")
                continue
            peaks_a = load_broadpeak(id_a)
            peaks_b = load_broadpeak(id_b)
            if peaks_a is None or peaks_b is None:
                missing = id_a if peaks_a is None else id_b
                print(f"    {mark}: broadPeak file missing/empty for {missing}")
                continue
            res = analyze_overlap(
                set(peaks_a['peak_id']), set(peaks_b['peak_id']),
                cond_a, cond_b, mark, rep, comp_label)
            all_results.append(res)
            print(f"    {mark}: {cond_a} {res['total_a']:,}  |  "
                  f"{cond_b} {res['total_b']:,}  |  "
                  f"Maintained {res['maintained']:,}  |  "
                  f"Lost {res['lost']:,}  |  Gained {res['gained']:,}")

# ── Plot: Venn diagrams for each comparison ───────────────────────────────────

COLORS = {
    'Control':     '#3498db',
    'ZL0454':      '#9b59b6',
    'TGFb':        '#e67e22',
    'TGFb+ZL0454': '#e74c3c',
}

for cond_a, cond_b, comp_label in COMPARISONS:
    subset = [r for r in all_results
              if r['comparison'] == comp_label]
    if not subset:
        continue

    # One row per replicate, one column per mark
    reps_present = sorted(set(r['replicate'] for r in subset))
    n_rows = len(reps_present)
    n_cols = len(MARKS)

    fig, axes = plt.subplots(n_rows, n_cols,
                             figsize=(4.5 * n_cols, 4.5 * n_rows))
    if n_rows == 1:
        axes = axes[np.newaxis, :]

    for ri, rep in enumerate(reps_present):
        for ci, mark in enumerate(MARKS):
            ax = axes[ri, ci]
            row = next((r for r in subset
                        if r['replicate'] == rep and r['mark'] == mark), None)
            if row is None:
                ax.set_visible(False)
                continue

            col_a = COLORS.get(cond_a, '#3498db')
            col_b = COLORS.get(cond_b, '#e74c3c')

            v = venn2(subsets=(row['lost'], row['gained'], row['maintained']),
                      set_labels=(cond_a, cond_b),
                      ax=ax, alpha=0.65,
                      set_colors=(col_a, col_b))
            venn2_circles(subsets=(row['lost'], row['gained'], row['maintained']),
                          ax=ax, linewidth=1.2, color='black')
            ax.set_title(f'{mark}  (Rep {rep})', fontsize=12, fontweight='bold')

            direction = 'GAIN' if row['gained'] > row['lost'] else 'LOSS'
            net = abs(row['gained'] - row['lost'])
            ax.text(0.5, -0.12, f'Net {direction}: {net:,}',
                    ha='center', transform=ax.transAxes,
                    fontsize=10, fontstyle='italic', fontweight='bold')

    fig.suptitle(f'Peak Overlaps — {comp_label}',
                 fontsize=15, fontweight='bold', y=1.01)
    plt.tight_layout()
    safe = comp_label.replace(' ', '_').replace('(', '').replace(')', '').replace('+','plus')
    fig.savefig(f'venn_{safe}.png', dpi=300, bbox_inches='tight')
    fig.savefig(f'venn_{safe}.pdf', bbox_inches='tight')
    plt.close(fig)
    print(f"\nSaved venn_{safe}.png")

# ── Plot: grouped bar chart summary ──────────────────────────────────────────

if all_results:
    df = pd.DataFrame(all_results)

    for cond_a, cond_b, comp_label in COMPARISONS:
        sub = df[df['comparison'] == comp_label]
        if sub.empty:
            continue

        reps_present = sorted(sub['replicate'].unique())
        fig, axes = plt.subplots(1, len(reps_present),
                                 figsize=(8 * len(reps_present), 6),
                                 sharey=False)
        if len(reps_present) == 1:
            axes = [axes]

        for ax, rep in zip(axes, reps_present):
            s = sub[sub['replicate'] == rep].set_index('mark').reindex(MARKS).dropna()
            if s.empty:
                ax.set_visible(False)
                continue
            x = np.arange(len(s))
            w = 0.25
            ax.bar(x - w, s['maintained'],      w, label='Maintained', color='#95a5a6', alpha=0.85)
            ax.bar(x,     s['lost'],             w, label=f'Lost (→{cond_b})',  color='#e74c3c', alpha=0.85)
            ax.bar(x + w, s['gained'],           w, label=f'Gained (→{cond_b})',color='#27ae60', alpha=0.85)
            ax.set_xticks(x); ax.set_xticklabels(s.index, fontsize=11)
            ax.set_ylabel('Number of Peaks', fontsize=11)
            ax.set_title(f'Replicate {rep}', fontsize=12, fontweight='bold')
            ax.legend(fontsize=10); ax.grid(axis='y', alpha=0.3)

        fig.suptitle(f'Peak Categories — {comp_label}',
                     fontsize=14, fontweight='bold')
        plt.tight_layout()
        safe = comp_label.replace(' ', '_').replace('(', '').replace(')', '').replace('+','plus')
        fig.savefig(f'bar_{safe}.png', dpi=300, bbox_inches='tight')
        fig.savefig(f'bar_{safe}.pdf', bbox_inches='tight')
        plt.close(fig)
        print(f"Saved bar_{safe}.png")

    # ── Reads-in-peaks bar chart (QC) ────────────────────────────────────────
    reads_data = {
        'IPMC1':346430,'IPMC2':833681,'IPMC3':1659168,'IPMC4':164883,'IPMC5':24679,
        'IPMC6':439987,'IPMC7':390341,'IPMC8':48471,'IPMC9':48471,'IPMC10':19403,
        'IPMC11':404062,'IPMC12':663201,'IPMC13':2645232,'IPMC14':353092,'IPMC15':25856,
        'IPMC16':384940,'IPMC17':756639,'IPMC18':3664937,'IPMC19':126077,'IPMC20':21892,
        'IPMC21':425118,'IPMC22':5686282,'IPMC23':6637950,'IPMC24':222342,'IPMC25':15612,
        'IPMC26':400744,'IPMC27':7912868,'IPMC28':6473573,'IPMC29':224597,'IPMC30':29692,
        'IPMC31':226663,'IPMC32':5312456,'IPMC33':5993150,'IPMC34':105922,'IPMC35':0,
        'IPMC36':367595,'IPMC37':4204744,'IPMC38':5498772,'IPMC39':80252,'IPMC40':22204,
    }

    rows = []
    for iid, reads in reads_data.items():
        rep, cond, mark = SAMPLE_MAP[iid]
        rows.append({'sample': iid, 'replicate': rep, 'condition': cond,
                     'mark': mark, 'reads_in_peaks': reads})
    qc = pd.DataFrame(rows)

    COND_ORDER = ['Control', 'ZL0454', 'TGFb', 'TGFb+ZL0454']
    COND_COLORS = ['#3498db', '#9b59b6', '#e67e22', '#e74c3c']

    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    for ai, rep in enumerate(['A', 'B']):
        ax = axes[ai]
        sub = qc[(qc['replicate'] == rep) & (qc['mark'] != 'input')]
        marks_in_rep = [m for m in MARKS if m in sub['mark'].values]
        x = np.arange(len(marks_in_rep))
        w = 0.18
        for ci, cond in enumerate(COND_ORDER):
            vals = [sub[(sub['mark'] == m) & (sub['condition'] == cond)]['reads_in_peaks'].values
                    for m in marks_in_rep]
            vals = [v[0] if len(v) else 0 for v in vals]
            ax.bar(x + (ci - 1.5) * w, vals, w,
                   label=cond, color=COND_COLORS[ci], alpha=0.85)
        ax.set_xticks(x); ax.set_xticklabels(marks_in_rep, fontsize=11)
        ax.set_ylabel('Reads in Peaks', fontsize=11)
        ax.set_title(f'Replicate {rep} — Reads in Peaks QC', fontsize=12, fontweight='bold')
        ax.legend(fontsize=10); ax.grid(axis='y', alpha=0.3)
        ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f'{int(v):,}'))

    fig.suptitle('ChIP-seq Quality Control: Reads in Peaks per Condition',
                 fontsize=14, fontweight='bold')
    plt.tight_layout()
    fig.savefig('reads_in_peaks_QC.png', dpi=300, bbox_inches='tight')
    fig.savefig('reads_in_peaks_QC.pdf', bbox_inches='tight')
    plt.close(fig)
    print("Saved reads_in_peaks_QC.png")

print("\n" + "="*70)
print("Done. Output files:")
for cond_a, cond_b, comp_label in COMPARISONS:
    safe = comp_label.replace(' ', '_').replace('(','').replace(')','').replace('+','plus')
    print(f"  venn_{safe}.png/pdf")
    print(f"  bar_{safe}.png/pdf")
print("  reads_in_peaks_QC.png/pdf")
print("="*70)
