#!/usr/bin/env python3
"""
Script to analyze peak overlaps and create Venn diagrams
for TGFβ vs TGFβ+BRD4i comparison

This script requires the narrowPeak files to identify specific peaks
that are gained, lost, or maintained with BRD4 inhibition
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib_venn import venn2, venn2_circles
import numpy as np

# File paths - adjust these to your actual paths
file_paths = {
    'TB1_H3K4me3': 'TGFB-TB1-H3K4me3_narrow_qval_peaks.narrowPeak',
    'IT1_H3K4me3': 'TGF-IT1-H3K4me3_narrow_qval_peaks.narrowPeak',
    'TB2_H3K4me1': 'TGFB-TB2-H3K4me1_narrow_qval_peaks.narrowPeak',
    'IT2_H3K4me1': 'TGF-IT2-H3K4me1_narrow_qval_peaks.narrowPeak',
    'TB3_H3K27ac': 'TGFB-TB3-H3K27ac_narrow_qval_peaks.narrowPeak',
    'IT3_H3K27ac': 'TGF-IT3-H3K27ac_narrow_qval_peaks.narrowPeak',
    'TB4_H3K122ac': 'TGFB-TB4-H3K122ac_narrow_qval_peaks.narrowPeak',
    'IT4_H3K122ac': 'TGF-IT4-H3K122ac_narrow_qval_peaks.narrowPeak',
}

def load_narrowpeak(filepath):
    """Load narrowPeak file and create genomic coordinate IDs"""
    try:
        df = pd.read_csv(filepath, sep='\t', header=None,
                        names=['chr', 'start', 'end', 'name', 'score',
                               'strand', 'signalValue', 'pValue', 'qValue', 'peak'])
        # Create unique peak ID based on genomic location (within 500bp window)
        df['peak_id'] = df['chr'] + ':' + (df['start']//500).astype(str)
        return df
    except FileNotFoundError:
        print(f"Warning: Could not find {filepath}")
        return None

def analyze_peak_overlaps(tb_peaks, it_peaks, mark_name):
    """Analyze overlapping and unique peaks between conditions"""
    if tb_peaks is None or it_peaks is None:
        return None

    tb_set = set(tb_peaks['peak_id'])
    it_set = set(it_peaks['peak_id'])

    overlap = tb_set & it_set
    tb_only = tb_set - it_set
    it_only = it_set - tb_set

    results = {
        'mark': mark_name,
        'tb_total': len(tb_set),
        'it_total': len(it_set),
        'maintained': len(overlap),
        'lost_with_brd4i': len(tb_only),
        'gained_with_brd4i': len(it_only),
        'overlap_pct_tb': len(overlap)/len(tb_set)*100 if len(tb_set) > 0 else 0,
        'overlap_pct_it': len(overlap)/len(it_set)*100 if len(it_set) > 0 else 0,
    }

    return results

# Main analysis
print("Analyzing peak overlaps between TGFβ and TGFβ+BRD4i conditions...")
print("="*80)

all_results = []
marks = ['H3K4me3', 'H3K4me1', 'H3K27ac', 'H3K122ac']

fig, axes = plt.subplots(2, 2, figsize=(14, 12))
axes = axes.flatten()

for idx, mark in enumerate(marks):
    tb_key = f'TB{idx+1}_{mark}'
    it_key = f'IT{idx+1}_{mark}'

    tb_peaks = load_narrowpeak(file_paths[tb_key])
    it_peaks = load_narrowpeak(file_paths[it_key])

    results = analyze_peak_overlaps(tb_peaks, it_peaks, mark)

    if results:
        all_results.append(results)

        # Print results
        print(f"\n{mark}:")
        print(f"  TGFβ only: {results['tb_total']:,} peaks")
        print(f"  TGFβ + BRD4i: {results['it_total']:,} peaks")
        print(f"  Maintained: {results['maintained']:,} peaks ({results['overlap_pct_tb']:.1f}% of TB)")
        print(f"  Lost with BRD4i: {results['lost_with_brd4i']:,} peaks")
        print(f"  Gained with BRD4i: {results['gained_with_brd4i']:,} peaks")

        # Create Venn diagram
        ax = axes[idx]
        venn = venn2(subsets=(results['lost_with_brd4i'],
                             results['gained_with_brd4i'],
                             results['maintained']),
                    set_labels=('TGFβ only', 'TGFβ + BRD4i'),
                    ax=ax, alpha=0.7,
                    set_colors=('#3498db', '#e74c3c'))

        venn2_circles(subsets=(results['lost_with_brd4i'],
                              results['gained_with_brd4i'],
                              results['maintained']),
                     ax=ax, linewidth=1.5, color='black')

        ax.set_title(f'{mark} Peak Overlap', fontsize=14, fontweight='bold', pad=20)

        # Add text annotations
        if results['gained_with_brd4i'] > results['lost_with_brd4i']:
            direction = "GAIN"
            net_change = results['gained_with_brd4i'] - results['lost_with_brd4i']
        else:
            direction = "LOSS"
            net_change = results['lost_with_brd4i'] - results['gained_with_brd4i']

        ax.text(0.5, -0.15, f'Net {direction}: {net_change:,} peaks',
                ha='center', transform=ax.transAxes, fontsize=11,
                fontweight='bold', style='italic')

plt.suptitle('Peak Overlap Analysis: Effect of BRD4 Inhibition on TGFβ Response',
             fontsize=16, fontweight='bold', y=0.98)
plt.tight_layout()
plt.savefig('peak_overlap_venn_diagrams.png', dpi=300, bbox_inches='tight')
plt.savefig('peak_overlap_venn_diagrams.pdf', bbox_inches='tight')
plt.show()

# Create summary bar plot
if all_results:
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

    results_df = pd.DataFrame(all_results)
    x = np.arange(len(results_df))
    width = 0.25

    # Plot 1: Absolute numbers
    ax1.bar(x - width, results_df['maintained'], width,
            label='Maintained', color='#95a5a6', alpha=0.8)
    ax1.bar(x, results_df['lost_with_brd4i'], width,
            label='Lost with BRD4i', color='#e74c3c', alpha=0.8)
    ax1.bar(x + width, results_df['gained_with_brd4i'], width,
            label='Gained with BRD4i', color='#27ae60', alpha=0.8)

    ax1.set_xlabel('Histone Mark', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Number of Peaks', fontsize=12, fontweight='bold')
    ax1.set_title('Peak Categories by BRD4 Inhibition Status',
                  fontsize=14, fontweight='bold')
    ax1.set_xticks(x)
    ax1.set_xticklabels(results_df['mark'])
    ax1.legend(fontsize=11)
    ax1.grid(axis='y', alpha=0.3)

    # Plot 2: Percentages
    total_peaks = results_df['tb_total'] + results_df['it_total'] - results_df['maintained']
    maintained_pct = results_df['maintained'] / total_peaks * 100
    lost_pct = results_df['lost_with_brd4i'] / total_peaks * 100
    gained_pct = results_df['gained_with_brd4i'] / total_peaks * 100

    ax2.bar(x - width, maintained_pct, width,
            label='Maintained', color='#95a5a6', alpha=0.8)
    ax2.bar(x, lost_pct, width,
            label='Lost with BRD4i', color='#e74c3c', alpha=0.8)
    ax2.bar(x + width, gained_pct, width,
            label='Gained with BRD4i', color='#27ae60', alpha=0.8)

    ax2.set_xlabel('Histone Mark', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Percentage of All Peaks (%)', fontsize=12, fontweight='bold')
    ax2.set_title('Peak Distribution (% of Total Unique Peaks)',
                  fontsize=14, fontweight='bold')
    ax2.set_xticks(x)
    ax2.set_xticklabels(results_df['mark'])
    ax2.legend(fontsize=11)
    ax2.grid(axis='y', alpha=0.3)

    plt.tight_layout()
    plt.savefig('peak_categories_comparison.png', dpi=300, bbox_inches='tight')
    plt.savefig('peak_categories_comparison.pdf', bbox_inches='tight')
    plt.show()

    print("\n" + "="*80)
    print("Analysis complete! Generated files:")
    print("  - peak_overlap_venn_diagrams.png/pdf")
    print("  - peak_categories_comparison.png/pdf")
    print("="*80)
else:
    print("\nNote: To run the full overlap analysis, ensure narrowPeak files are in the working directory")
    print("The script can still generate summary statistics from the annotation data.")
