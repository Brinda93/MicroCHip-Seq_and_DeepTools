#Integrating two ChIP-seq batches and performing peak-gene overlap analysis with RNA-seq

#!/usr/bin/env python3
import subprocess, sys, re
subprocess.run([sys.executable, '-m', 'pip', 'install', 'matplotlib-venn', 'openpyxl', '-q'])

import os, io, warnings
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib_venn import venn2
from google.colab import files
warnings.filterwarnings('ignore')

# ── FILE_MAP ─────────────────────────────────────────────────────────────────
FILE_MAP = {
    'TGFB-TB1-H3K4me3_narrow_qval_peaks_HOMER_annot.txt':  ('TGFb',       'H3K4me3',  'Jan'),
    'TGFB-TB2-H3K4me1_narrow_qval_peaks_HOMER_annot.txt':  ('TGFb',       'H3K4me1',  'Jan'),
    'TGFB-TB3-H3K27ac_narrow_qval_peaks_HOMER_annot.txt':  ('TGFb',       'H3K27ac',  'Jan'),
    'TGFB-TB4-H3K122ac_narrow_qval_peaks_HOMER_annot.txt': ('TGFb',       'H3K122ac', 'Jan'),
    'TGF-IT1-H3K4me3_narrow_qval_peaks_HOMER_annot.txt':   ('TGFb_BRD4i', 'H3K4me3',  'Jan'),
    'TGF-IT2-H3K4me1_narrow_qval_peaks_HOMER_annot.txt':   ('TGFb_BRD4i', 'H3K4me1',  'Jan'),
    'TGF-IT3-H3K27ac_narrow_qval_peaks_HOMER_annot.txt':   ('TGFb_BRD4i', 'H3K27ac',  'Jan'),
    'TGF-IT4-H3K122ac_narrow_qval_peaks_HOMER_annot.txt':  ('TGFb_BRD4i', 'H3K122ac', 'Jan'),
    'IPMC32_peaks.annotated.txt':  ('TGFb',       'H3K4me3',  '3B'),
    'IPMC33_peaks.annotated.txt':  ('TGFb',       'H3K4me1',  '3B'),
    'IPMC34_peaks.annotated.txt':  ('TGFb',       'H3K27ac',  '3B'),
    'IPMC37_peaks.annotated.txt':  ('TGFb_BRD4i', 'H3K4me3',  '4B'),
    'IPMC38_peaks.annotated.txt':  ('TGFb_BRD4i', 'H3K4me1',  '4B'),
    'IPMC39_peaks.annotated.txt':  ('TGFb_BRD4i', 'H3K27ac',  '4B'),
    'IPMC40_peaks.annotated.txt':  ('TGFb_BRD4i', 'H3K122ac', '4B'),
}

MARKS    = ['H3K4me3', 'H3K4me1', 'H3K27ac', 'H3K122ac']
COLORS   = {'TGFb': '#3498db', 'TGFb_BRD4i': '#e74c3c'}
P_THRESH = 0.05   # raw p-value threshold
FC_THRESH = 0.0   # set e.g. 1.0 to also require |log2FC| > 1

# ── HELPERS ──────────────────────────────────────────────────────────────────
def strip_colab_suffix(fname):
    return re.sub(r'\s*\(\d+\)(\.[^.]+)?$', lambda m: m.group(1) or '', fname)

def load_homer_bytes(data_bytes):
    df = pd.read_csv(io.BytesIO(data_bytes), sep='\t', header=0, low_memory=False)
    df.rename(columns={df.columns[0]: 'PeakID'}, inplace=True)
    return df

def get_genes(df, annot_types=None):
    if df is None or df.empty: return set()
    sub = df.copy()
    if annot_types:
        sub = sub[sub['Annotation'].str.contains('|'.join(annot_types), case=False, na=False)]
    col = next((c for c in ['Gene Name', 'Nearest Gene', 'Gene Symbol'] if c in sub.columns), None)
    return set(sub[col].dropna().unique()) if col else set()

def load_deg_bytes(data_bytes):
    for sep in ['\t', ',']:
        try:
            df = pd.read_csv(io.BytesIO(data_bytes), sep=sep)
            if len(df.columns) > 3: break
        except Exception: continue
    # Standardise column names
    sig_col = next((c for c in df.columns if c.lower() in ('sign', 'diff_expressed')), None)
    if sig_col and sig_col != 'sign': df.rename(columns={sig_col: 'sign'}, inplace=True)
    fc  = next((c for c in df.columns if 'log2fold' in c.lower()), 'log2FoldChange')
    adj = next((c for c in df.columns if 'padj'     in c.lower()), 'padj')
    pv  = next((c for c in df.columns if c.lower()  == 'pvalue'),  'pvalue')
    if fc  != 'log2FoldChange': df.rename(columns={fc:  'log2FoldChange'}, inplace=True)
    if adj != 'padj':           df.rename(columns={adj: 'padj'},           inplace=True)
    if pv  != 'pvalue':         df.rename(columns={pv:  'pvalue'},         inplace=True)
    return df

def gene_table(genes, deg):
    if not genes:
        return pd.DataFrame(columns=['Gene', 'log2FoldChange', 'pvalue', 'padj', 'baseMean', 'sign'])
    cols = [c for c in ['gene_id', 'log2FoldChange', 'pvalue', 'padj', 'baseMean', 'sign'] if c in deg.columns]
    return (pd.DataFrame({'Gene': list(genes)})
              .merge(deg[cols], left_on='Gene', right_on='gene_id', how='left')
              .drop(columns='gene_id', errors='ignore')
              .sort_values('log2FoldChange', na_position='last'))

# ════════════════════════════════════════════════════════════════════════════
# STEP 1 — Upload annotation files
# ════════════════════════════════════════════════════════════════════════════
print("=" * 65)
print("ANALYSIS 1 — Non-inflamed (A0): TGFβ vs TGFβ+BRD4i")
print("Combining January batch + new 3B/4B batch")
print("=" * 65)
print("\nSTEP 1: Upload ALL 15 annotation files now...")
uploaded_annot = files.upload()

# ════════════════════════════════════════════════════════════════════════════
# STEP 2 — Parse annotation files
# ════════════════════════════════════════════════════════════════════════════
print("\nSTEP 2: Parsing annotation files...")
pool = {m: {c: {'all': set(), 'promoter': set()}
            for c in ['TGFb', 'TGFb_BRD4i']} for m in MARKS}

for fname, data_bytes in uploaded_annot.items():
    clean = strip_colab_suffix(fname)
    if clean not in FILE_MAP:
        print(f"  ⚠  Unrecognised (skipped): {fname}  [cleaned: {clean}]")
        continue
    cond, mark, batch = FILE_MAP[clean]
    df = load_homer_bytes(data_bytes)
    pool[mark][cond]['all']      |= get_genes(df)
    pool[mark][cond]['promoter'] |= get_genes(df, ['promoter'])
    print(f"  ✔ {clean}  →  {cond} | {mark} | {batch} | {len(df):,} peaks")

print("\nPooled gene counts:")
for m in MARKS:
    for c in ['TGFb', 'TGFb_BRD4i']:
        print(f"  {m:10s}  {c:12s}  all={len(pool[m][c]['all']):>6,}  "
              f"promoter={len(pool[m][c]['promoter']):>5,}")

# ════════════════════════════════════════════════════════════════════════════
# STEP 3 — Upload DEG file
# ════════════════════════════════════════════════════════════════════════════
print("\nSTEP 3: Upload DEG CSV file...")
uploaded_deg = files.upload()
deg_fname = list(uploaded_deg.keys())[0]
deg = load_deg_bytes(uploaded_deg[deg_fname])

# Use raw p-value < 0.05 for up/down gene sets
sig = deg[deg['pvalue'] < P_THRESH].copy()
if FC_THRESH > 0:
    sig = sig[sig['log2FoldChange'].abs() > FC_THRESH]

up_g   = set(sig[sig['log2FoldChange'] > 0]['gene_id'])
down_g = set(sig[sig['log2FoldChange'] < 0]['gene_id'])
print(f"  Loaded: {deg_fname}")
print(f"  Threshold: pvalue < {P_THRESH}" + (f" & |log2FC| > {FC_THRESH}" if FC_THRESH > 0 else ""))
print(f"  Upregulated: {len(up_g):,}   Downregulated: {len(down_g):,}")

# ════════════════════════════════════════════════════════════════════════════
# STEP 4 — Integrate ChIP-seq × RNA-seq
# ════════════════════════════════════════════════════════════════════════════
print("\nSTEP 4: Integrating ChIP-seq × RNA-seq...")
results = []
for m in MARKS:
    tb  = pool[m]['TGFb']['all'];       it  = pool[m]['TGFb_BRD4i']['all']
    tbp = pool[m]['TGFb']['promoter'];  itp = pool[m]['TGFb_BRD4i']['promoter']
    tb_only = tb - it;  it_only = it - tb;  shared = tb & it
    r = {
        'mark': m,
        'TB_genes':         tb,
        'IT_genes':         it,
        'TB_only':          tb_only,
        'IT_only':          it_only,
        'Shared':           shared,
        'direct_targets':   tb_only & down_g,
        'gained_activated': it_only & up_g,
        'indirect':         tb_only & up_g,
        'shared_up':        shared  & up_g,
        'shared_down':      shared  & down_g,
        'prom_direct':      (tbp - itp) & down_g,
        'prom_gained':      (itp - tbp) & up_g,
    }
    results.append(r)
    print(f"  {m:10s}  direct={len(r['direct_targets']):>4}  "
          f"gained+up={len(r['gained_activated']):>4}  "
          f"indirect={len(r['indirect']):>4}  "
          f"prom_direct={len(r['prom_direct']):>4}")

# ════════════════════════════════════════════════════════════════════════════
# STEP 5 — Save Excel
# ════════════════════════════════════════════════════════════════════════════
print("\nSTEP 5: Saving Excel...")
excel_file = 'Analysis1_ChIPseq_RNAseq_Integration.xlsx'
with pd.ExcelWriter(excel_file, engine='openpyxl') as writer:

    # Summary sheet
    summ = [{
        'Mark':             r['mark'],
        'TGFb_genes':       len(r['TB_genes']),
        'TGFb_BRD4i_genes': len(r['IT_genes']),
        'TB_only_lost':     len(r['TB_only']),
        'IT_only_gained':   len(r['IT_only']),
        'Shared':           len(r['Shared']),
        'Direct_targets':   len(r['direct_targets']),
        'Gained_Activated': len(r['gained_activated']),
        'Indirect':         len(r['indirect']),
        'Promoter_Direct':  len(r['prom_direct']),
        'Promoter_Gained':  len(r['prom_gained']),
    } for r in results]
    pd.DataFrame(summ).to_excel(writer, sheet_name='Summary', index=False)

    # Per-mark gene sheets
    for r in results:
        m = r['mark']
        for sheet, genes in [
            (f'{m}_DirectTargets',  r['direct_targets']),
            (f'{m}_GainedUp',       r['gained_activated']),
            (f'{m}_Indirect',       r['indirect']),
            (f'{m}_PromDirect',     r['prom_direct']),
            (f'{m}_PromGained',     r['prom_gained']),
            (f'{m}_TGFb_all',       r['TB_genes']),
            (f'{m}_BRD4i_all',      r['IT_genes']),
        ]:
            gene_table(genes, deg).to_excel(writer, sheet_name=sheet[:31], index=False)

    # Flat integrated table
    rows = []
    for r in results:
        for gene in r['direct_targets']:
            rows.append({'Gene': gene, 'Mark': r['mark'], 'Category': 'Direct_Target'})
        for gene in r['gained_activated']:
            rows.append({'Gene': gene, 'Mark': r['mark'], 'Category': 'Gained_Activated'})
        for gene in r['indirect']:
            rows.append({'Gene': gene, 'Mark': r['mark'], 'Category': 'Indirect'})

    if rows:
        flat = (pd.DataFrame(rows)
                  .merge(deg[['gene_id', 'log2FoldChange', 'pvalue', 'padj', 'baseMean', 'sign']],
                         left_on='Gene', right_on='gene_id', how='left')
                  .drop(columns='gene_id', errors='ignore'))
    else:
        flat = pd.DataFrame(columns=['Gene', 'Mark', 'Category', 'log2FoldChange', 'pvalue', 'padj', 'baseMean', 'sign'])
    flat.to_excel(writer, sheet_name='All_Integrated_Genes', index=False)

files.download(excel_file)
print(f"  ✔ Downloaded: {excel_file}")

# ════════════════════════════════════════════════════════════════════════════
# STEP 6 — Figures
# ════════════════════════════════════════════════════════════════════════════
print("\nSTEP 6: Creating figures...")
sns.set_style('whitegrid')

# Fig 1: Peak gene counts
fig, ax = plt.subplots(figsize=(10, 6))
x = np.arange(len(MARKS)); w = 0.35
tb_c = [len(next(r for r in results if r['mark']==m)['TB_genes']) for m in MARKS]
it_c = [len(next(r for r in results if r['mark']==m)['IT_genes']) for m in MARKS]
b1 = ax.bar(x-w/2, tb_c, w, label='TGFβ',        color=COLORS['TGFb'],       alpha=0.85, edgecolor='k', lw=1.2)
b2 = ax.bar(x+w/2, it_c, w, label='TGFβ+BRD4i',  color=COLORS['TGFb_BRD4i'],alpha=0.85, edgecolor='k', lw=1.2)


ax.bar_label(b1, labels=[f'{v:,}' for v in tb_c], fontsize=9, fontweight='bold')
ax.bar_label(b2, labels=[f'{v:,}' for v in it_c], fontsize=9, fontweight='bold')

ax.set_xticks(x); ax.set_xticklabels(MARKS, fontsize=12)
ax.set_ylabel('Genes with peaks', fontsize=12, fontweight='bold')
ax.set_title('Combined Peak Gene Counts — Non-inflamed (A0)', fontsize=13, fontweight='bold', pad=14)
ax.legend(fontsize=11); plt.tight_layout()
plt.savefig('fig1_peak_counts.png', dpi=300, bbox_inches='tight'); plt.show(); plt.close()
files.download('fig1_peak_counts.png'); print('  ✔ fig1_peak_counts.png')

# Fig 2: Integration summary
fig, ax = plt.subplots(figsize=(11, 6))
x = np.arange(len(MARKS)); w = 0.25
dt  = [len(next(r for r in results if r['mark']==m)['direct_targets'])   for m in MARKS]
gu  = [len(next(r for r in results if r['mark']==m)['gained_activated']) for m in MARKS]
ind = [len(next(r for r in results if r['mark']==m)['indirect'])         for m in MARKS]
b1 = ax.bar(x-w, dt,  w, label='Direct (lost+down)',  color='#e74c3c', alpha=0.85, edgecolor='k', lw=1.2)
b2 = ax.bar(x,   gu,  w, label='Gained+Up',            color='#27ae60', alpha=0.85, edgecolor='k', lw=1.2)
b3 = ax.bar(x+w, ind, w, label='Indirect (lost+up)',   color='#f39c12', alpha=0.85, edgecolor='k', lw=1.2)
for b in [b1, b2, b3]: ax.bar_label(b, fontsize=9, fontweight='bold')
ax.set_xticks(x); ax.set_xticklabels(MARKS, fontsize=12)
ax.set_ylabel('Number of genes', fontsize=12, fontweight='bold')
ax.set_title('BRD4-Dependent Gene Regulation — Non-inflamed (A0)', fontsize=13, fontweight='bold', pad=14)

ax.bar_label(b1, labels=[f'{v:,}' for v in dt],  fontsize=9, fontweight='bold')
ax.bar_label(b2, labels=[f'{v:,}' for v in gu],  fontsize=9, fontweight='bold')
ax.bar_label(b3, labels=[f'{v:,}' for v in ind], fontsize=9, fontweight='bold')

ax.legend(fontsize=11); plt.tight_layout()
plt.savefig('fig2_integration_summary.png', dpi=300, bbox_inches='tight'); plt.show(); plt.close()
files.download('fig2_integration_summary.png'); print('  ✔ fig2_integration_summary.png')

# Fig 3: Venn diagrams
fig, axes = plt.subplots(1, len(MARKS), figsize=(5*len(MARKS), 5))
for ax, r in zip(axes, results):
    venn2([r['TB_genes'], r['IT_genes']],
          set_labels=('TGFβ', 'TGFβ+BRD4i'),
          set_colors=(COLORS['TGFb'], COLORS['TGFb_BRD4i']),
          alpha=0.6, ax=ax)
    ax.set_title(r['mark'], fontsize=13, fontweight='bold')
fig.suptitle('Peak Gene Overlap — Non-inflamed (A0)', fontsize=14, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig('fig3_venn.png', dpi=300, bbox_inches='tight'); plt.show(); plt.close()
files.download('fig3_venn.png'); print('  ✔ fig3_venn.png')

# Fig 4–7: Volcano plots per mark
for r in results:
    m = r['mark']
    fig, ax = plt.subplots(figsize=(12, 9))

    # Background
    ax.scatter(deg['log2FoldChange'], -np.log10(deg['pvalue'].clip(lower=1e-300)),
               s=10, alpha=0.25, color='lightgray', rasterized=True, label='All genes')

    # Highlight each category
    for genes, color, label in [
        (r['direct_targets'],   '#e74c3c', f'Direct targets (n={len(r["direct_targets"])})'),
        (r['gained_activated'], '#27ae60', f'Gained+Up (n={len(r["gained_activated"])})'),
        (r['indirect'],         '#f39c12', f'Indirect (n={len(r["indirect"])})'),
    ]:
        sub = deg[deg['gene_id'].isin(genes)]
        if not sub.empty:
            ax.scatter(sub['log2FoldChange'], -np.log10(sub['pvalue'].clip(lower=1e-300)),
                       s=65, alpha=0.9, color=color, edgecolor='k', lw=0.5, label=label, zorder=5)

    # Label top genes from ALL categories (not just direct targets)
    label_candidates = pd.concat([
        deg[deg['gene_id'].isin(r['direct_targets'])].nsmallest(5, 'pvalue'),
        deg[deg['gene_id'].isin(r['gained_activated'])].nsmallest(5, 'pvalue'),
        deg[deg['gene_id'].isin(r['indirect'])].nsmallest(5, 'pvalue'),
    ]).drop_duplicates(subset='gene_id')

    # Color map for labels
    label_color_map = {}
    for gene in r['direct_targets']:   label_color_map[gene] = '#e74c3c'
    for gene in r['gained_activated']: label_color_map[gene] = '#27ae60'
    for gene in r['indirect']:         label_color_map[gene] = '#f39c12'

    from adjustText import adjust_text
    texts = []
    for _, row in label_candidates.iterrows():
        gid = row['gene_id']
        x   = row['log2FoldChange']
        y   = -np.log10(row['pvalue'] + 1e-300)
        fc  = label_color_map.get(gid, 'black')
        texts.append(ax.text(x, y, gid, fontsize=8, fontweight='bold', color='black',
                             bbox=dict(boxstyle='round,pad=0.25', fc=fc, alpha=0.4, ec='none')))

    try:
        adjust_text(texts, ax=ax,
                    arrowprops=dict(arrowstyle='-', color='gray', lw=0.6),
                    expand_points=(1.5, 1.5), expand_text=(1.4, 1.4))
    except Exception:
        pass  # adjustText optional — falls back to default positions

    ax.axhline(-np.log10(P_THRESH), color='k', ls='--', lw=1, alpha=0.5, label=f'pvalue={P_THRESH}')
    ax.axvline(0, color='k', ls='-', lw=0.8, alpha=0.4)
    ax.set_xlabel('log2FC (TGFβ+BRD4i vs TGFβ)', fontsize=12, fontweight='bold')
    ax.set_ylabel('-log10(pvalue)', fontsize=12, fontweight='bold')
    ax.set_title(f'{m} — Non-inflamed (A0)', fontsize=13, fontweight='bold', pad=14)
    ax.legend(fontsize=10, loc='upper right'); ax.grid(alpha=0.25); plt.tight_layout()
    fname = f'fig_volcano_{m}.png'
    plt.savefig(fname, dpi=300, bbox_inches='tight'); plt.show(); plt.close()
    files.download(fname); print(f'  ✔ {fname}')

    
# Fig 8: Promoter targets
fig, ax = plt.subplots(figsize=(9, 5))
pd_v = [len(next(r for r in results if r['mark']==m)['prom_direct'])  for m in MARKS]
pg_v = [len(next(r for r in results if r['mark']==m)['prom_gained']) for m in MARKS]
x = np.arange(len(MARKS)); w = 0.35
b1 = ax.bar(x-w/2, pd_v, w, label='Promoter Direct (lost+down)', color='#9b59b6', alpha=0.85, edgecolor='k', lw=1.2)
b2 = ax.bar(x+w/2, pg_v, w, label='Promoter Gained+Up',          color='#1abc9c', alpha=0.85, edgecolor='k', lw=1.2)
for b in [b1, b2]: ax.bar_label(b, fontsize=9, fontweight='bold')
ax.set_xticks(x); ax.set_xticklabels(MARKS, fontsize=12)
ax.set_ylabel('Number of genes', fontsize=12, fontweight='bold')
ax.set_title('Promoter-Associated BRD4 Targets — Non-inflamed (A0)', fontsize=13, fontweight='bold', pad=14)

ax.bar_label(b1, labels=[f'{v:,}' for v in pd_v], fontsize=9, fontweight='bold')
ax.bar_label(b2, labels=[f'{v:,}' for v in pg_v], fontsize=9, fontweight='bold')

ax.legend(fontsize=11); plt.tight_layout()
plt.savefig('fig8_promoter_targets.png', dpi=300, bbox_inches='tight'); plt.show(); plt.close()
files.download('fig8_promoter_targets.png'); print('  ✔ fig8_promoter_targets.png')

print("\n" + "=" * 65)
print("ANALYSIS 1 COMPLETE")
print("All Excel and PNG files downloaded to your computer.")
print("=" * 65)
