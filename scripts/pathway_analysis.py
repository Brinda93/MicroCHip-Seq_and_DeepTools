#!/usr/bin/env python3
"""
Pathway Enrichment Analysis for BRD4-Dependent TGFβ Targets
Tailored to ChIPseq_RNAseq_Integration_Analysis.xlsx output structure.

Marks analysed : H3K4me3, H3K4me1, H3K27ac, H3K122ac
Gene-set categories used:
  • DirectTargets      – lost mark + downregulated (BRD4-dependent activation)
  • GainedActivated    – gained mark + upregulated
  • IndirectActivated  – lost mark + upregulated
  • PromoterDirect     – promoter-level direct targets (lost + down)
  • PromoterGained     – promoter-level gained + up

Enrichment is run via the Enrichr REST API (no key required).
"""

# ── Imports ──────────────────────────────────────────────────────────────────
import pandas as pd
import matplotlib
matplotlib.use('Agg')          # safe for Colab / headless
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import numpy as np
import requests, json, time, os, warnings
warnings.filterwarnings('ignore')

from google.colab import files

print("=" * 80)
print("Pathway Enrichment Analysis — BRD4-Dependent TGFβ Targets")
print("=" * 80)

# ════════════════════════════════════════════════════════════════════════════
# STEP 1 — Upload integration Excel
# ════════════════════════════════════════════════════════════════════════════
print("\nSTEP 1: Upload ChIPseq_RNAseq_Integration_Analysis.xlsx")
uploaded = files.upload()
INTEGRATION_FILE = list(uploaded.keys())[0]

xl = pd.ExcelFile(INTEGRATION_FILE)
print(f"Sheets found: {xl.sheet_names}")

# ════════════════════════════════════════════════════════════════════════════
# STEP 2 — Extract gene lists
# ════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 80)
print("STEP 2: Extracting gene lists")
print("=" * 80)

# Map sheet suffixes → clean category labels
SHEET_CATEGORY_MAP = {
    'DirectTargets':    'DirectTargets',
    'GainedUp':         'GainedActivated',
    'IndirectUp':       'IndirectActivated',
    'PromDirectTargets':'PromoterDirect',
    'PromGainedUp':     'PromoterGained',
}

gene_sets: dict[str, list[str]] = {}   # key = "MARK_Category"

for sheet in xl.sheet_names:
    if sheet in ('Summary', 'All_Integrated_Genes'):
        continue
    for suffix, label in SHEET_CATEGORY_MAP.items():
        if suffix in sheet:
            mark = sheet.replace(suffix, '').rstrip('_')
            df   = pd.read_excel(INTEGRATION_FILE, sheet_name=sheet)
            if 'Gene' not in df.columns:
                print(f"  ⚠ No 'Gene' column in {sheet}, skipping")
                break
            genes = df['Gene'].dropna().astype(str).str.strip().str.upper()
            genes = [g for g in genes if g and g != 'NAN']
            if genes:
                key = f"{mark}_{label}"
                gene_sets[key] = genes
                print(f"  ✔ {key:45s} → {len(genes):>3} genes")
            break

# Also pull from All_Integrated_Genes if present (covers H3K122ac edge cases)
if 'All_Integrated_Genes' in xl.sheet_names:
    agi = pd.read_excel(INTEGRATION_FILE, sheet_name='All_Integrated_Genes')
    cat_map = {
        'Direct_Target':     'DirectTargets',
        'Gained_Activated':  'GainedActivated',
        'Indirect_Activated':'IndirectActivated',
    }
    for raw_cat, label in cat_map.items():
        sub = agi[agi['Category'] == raw_cat]
        for mark in sub['Mark'].unique():
            key = f"{mark}_{label}"
            if key not in gene_sets:
                genes = (sub[sub['Mark'] == mark]['Gene']
                         .dropna().astype(str).str.upper().str.strip().tolist())
                genes = [g for g in genes if g and g != 'NAN']
                if genes:
                    gene_sets[key] = genes
                    print(f"  ✔ {key:45s} → {len(genes):>3} genes  (from All_Integrated_Genes)")

if not gene_sets:
    raise SystemExit("No gene sets extracted — check sheet structure.")

print(f"\nTotal gene sets: {len(gene_sets)}")

# ════════════════════════════════════════════════════════════════════════════
# STEP 3 — Enrichr API helpers
# ════════════════════════════════════════════════════════════════════════════
ENRICHR_ADD   = 'https://maayanlab.cloud/Enrichr/addList'
ENRICHR_ENRICH = 'https://maayanlab.cloud/Enrichr/enrich'

LIBRARIES = [
    'GO_Biological_Process_2023',
    'GO_Molecular_Function_2023',
    'KEGG_2021_Human',
    'WikiPathway_2023_Human',
    'Reactome_2022',
    'MSigDB_Hallmark_2020',
]

def submit_to_enrichr(genes: list[str], description: str) -> str | None:
    try:
        r = requests.post(
            ENRICHR_ADD,
            files={'list': (None, '\n'.join(genes)),
                   'description': (None, description)},
            timeout=30)
        if r.ok:
            return r.json().get('userListId')
    except Exception as e:
        print(f"    Submit error: {e}")
    return None


def get_enrichr_results(list_id: str, library: str) -> pd.DataFrame:
    try:
        r = requests.get(
            f"{ENRICHR_ENRICH}?userListId={list_id}&backgroundType={library}",
            timeout=30)
        if r.ok:
            data = r.json().get(library, [])
            if data:
                df = pd.DataFrame(data, columns=[
                    'Rank', 'Term', 'P-value', 'Odds_Ratio', 'Combined_Score',
                    'Genes', 'Adj_Pvalue', 'Old_Pvalue', 'Old_Adj_Pvalue'])
                return df[df['Adj_Pvalue'] < 0.05].head(25)
    except Exception as e:
        print(f"    Fetch error: {e}")
    return pd.DataFrame()


# ════════════════════════════════════════════════════════════════════════════
# STEP 4 — Run enrichment
# ════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 80)
print("STEP 3: Running Enrichr (this may take a few minutes)")
print("=" * 80)

MIN_GENES = 2   # skip gene sets that are too small

all_results: dict[str, dict[str, pd.DataFrame]] = {}

for gs_name, genes in gene_sets.items():
    if len(genes) < MIN_GENES:
        print(f"\n  Skipping {gs_name}: only {len(genes)} gene(s)")
        continue

    print(f"\n  ► {gs_name}  ({len(genes)} genes)")
    list_id = submit_to_enrichr(genes, gs_name)
    if not list_id:
        print("    ✗ Submission failed")
        continue

    lib_results = {}
    for lib in LIBRARIES:
        df = get_enrichr_results(list_id, lib)
        if not df.empty:
            lib_results[lib] = df
            print(f"    {lib:40s} → {len(df):>2} sig. terms")
        else:
            print(f"    {lib:40s} → —")
        time.sleep(0.4)

    if lib_results:
        all_results[gs_name] = lib_results

print(f"\nGene sets with significant enrichment: {len(all_results)}")

# ════════════════════════════════════════════════════════════════════════════
# STEP 5 — Save Excel
# ════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 80)
print("STEP 4: Saving enrichment Excel")
print("=" * 80)

ENRICH_EXCEL = 'Pathway_Enrichment_Results.xlsx'

with pd.ExcelWriter(ENRICH_EXCEL, engine='openpyxl') as writer:

    # ── Summary sheet ────────────────────────────────────────────────────────
    rows = []
    for gs, libs in all_results.items():
        for lib, df in libs.items():
            rows.append({
                'Gene_Set':         gs,
                'Library':          lib,
                'Sig_Terms':        len(df),
                'Top_Term':         df.iloc[0]['Term'] if not df.empty else '',
                'Top_Adj_Pval':     df.iloc[0]['Adj_Pvalue'] if not df.empty else 1,
                'Top_Combined_Score': df.iloc[0]['Combined_Score'] if not df.empty else 0,
            })
    if rows:
        pd.DataFrame(rows).to_excel(writer, sheet_name='Enrichment_Summary', index=False)
        print("✔ Sheet: Enrichment_Summary")

    # ── Per gene-set / library sheets ────────────────────────────────────────
    for gs, libs in all_results.items():
        for lib, df in libs.items():
            if df.empty:
                continue
            # Abbreviate library name to keep sheet name ≤ 31 chars
            lib_short = (lib
                         .replace('GO_Biological_Process_2023', 'GOBP')
                         .replace('GO_Molecular_Function_2023', 'GOMF')
                         .replace('KEGG_2021_Human',            'KEGG')
                         .replace('WikiPathway_2023_Human',     'Wiki')
                         .replace('Reactome_2022',              'React')
                         .replace('MSigDB_Hallmark_2020',       'Hallmark'))
            sheet = f"{gs}_{lib_short}"[:31]
            df.to_excel(writer, sheet_name=sheet, index=False)
            print(f"✔ Sheet: {sheet}")

print(f"\n✔ Saved: {ENRICH_EXCEL}")
files.download(ENRICH_EXCEL)

# ════════════════════════════════════════════════════════════════════════════
# STEP 6 — Figures
# ════════════════════════════════════════════════════════════════════════════
if not all_results:
    print("\n⚠ No enrichment results — skipping figures.")
else:
    print("\n" + "=" * 80)
    print("STEP 5: Generating figures")
    print("=" * 80)

    # ── Colour palette (one colour per histone mark) ─────────────────────────
    MARK_COLOURS = {
        'H3K4me3':  '#3498db',
        'H3K4me1':  '#e74c3c',
        'H3K27ac':  '#2ecc71',
        'H3K122ac': '#9b59b6',
    }

    def mark_colour(gs_name: str) -> str:
        for mark, col in MARK_COLOURS.items():
            if gs_name.startswith(mark):
                return col
        return '#95a5a6'

    # ── Helper: horizontal bar plot ──────────────────────────────────────────
    def barh_enrichment(ax, df: pd.DataFrame, title: str,
                        n: int = 15, colour: str = '#3498db'):
        df = df.copy().nlargest(n, 'Combined_Score')
        terms  = [t[:65] + '…' if len(t) > 65 else t for t in df['Term']]
        scores = df['Combined_Score'].values
        ypos   = np.arange(len(terms))
        bars   = ax.barh(ypos, scores, color=colour, alpha=0.82,
                         edgecolor='k', linewidth=0.8)
        ax.set_yticks(ypos)
        ax.set_yticklabels(terms, fontsize=8)
        ax.set_xlabel('Combined Score', fontsize=10, fontweight='bold')
        ax.set_title(title, fontsize=10, fontweight='bold', pad=6)
        ax.grid(axis='x', alpha=0.3)
        for bar, val in zip(bars, scores):
            ax.text(bar.get_width() + 0.5, bar.get_y() + bar.get_height()/2,
                    f'{val:.1f}', va='center', fontsize=7)

    # ════════════════════════════════════════════════════════════════════════
    # Figure 1 — Per gene-set bar charts (one row per gene set, one col per lib)
    # ════════════════════════════════════════════════════════════════════════
    gs_list  = list(all_results.keys())
    lib_list = LIBRARIES
    n_gs     = len(gs_list)

    fig, axes = plt.subplots(n_gs, 1, figsize=(14, 6 * n_gs))
    if n_gs == 1:
        axes = [axes]

    for ax, gs_name in zip(axes, gs_list):
        # Combine top 5 from each library
        frames = []
        for lib, df in all_results[gs_name].items():
            top = df.head(5).copy()
            top['Library'] = lib
            frames.append(top)
        if not frames:
            ax.axis('off'); continue

        combined = (pd.concat(frames)
                      .sort_values('Combined_Score')
                      .tail(15))

        libs_present = combined['Library'].unique()
        cmap   = plt.cm.Set2(np.linspace(0, 1, max(len(libs_present), 1)))
        cmap_d = {l: cmap[i] for i, l in enumerate(libs_present)}
        colours = [cmap_d[l] for l in combined['Library']]

        ypos = np.arange(len(combined))
        ax.barh(ypos, combined['Combined_Score'], color=colours,
                edgecolor='k', linewidth=0.8, alpha=0.85)
        terms = [t[:70] + '…' if len(t) > 70 else t for t in combined['Term']]
        ax.set_yticks(ypos); ax.set_yticklabels(terms, fontsize=8)
        ax.set_xlabel('Combined Score', fontsize=10, fontweight='bold')
        ax.set_title(f'Top Enriched Pathways — {gs_name}',
                     fontsize=11, fontweight='bold')
        ax.grid(axis='x', alpha=0.3)
        patches = [mpatches.Patch(color=cmap_d[l],
                   label=l.replace('_2023','').replace('_2021','')
                          .replace('_2022','').replace('_2020','')
                          .replace('_Human',''))
                   for l in libs_present]
        ax.legend(handles=patches, loc='lower right', fontsize=7, framealpha=0.9)

    plt.suptitle('Pathway Enrichment — BRD4-Dependent TGFβ Gene Sets',
                 fontsize=14, fontweight='bold', y=1.002)
    plt.tight_layout()
    fname = 'enrichment_all_genesets.png'
    plt.savefig(fname, dpi=300, bbox_inches='tight')
    files.download(fname)
    plt.close()
    print(f"✔ Downloaded: {fname}")

    # ════════════════════════════════════════════════════════════════════════
    # Figure 2 — Individual plots per gene set (GO BP + KEGG + Reactome)
    # ════════════════════════════════════════════════════════════════════════
    FOCUS_LIBS = ['GO_Biological_Process_2023', 'KEGG_2021_Human', 'Reactome_2022']

    for gs_name, libs in all_results.items():
        focus = {l: libs[l] for l in FOCUS_LIBS if l in libs and not libs[l].empty}
        if not focus:
            continue

        ncols = len(focus)
        fig, axes = plt.subplots(1, ncols, figsize=(7 * ncols, 9))
        if ncols == 1:
            axes = [axes]

        col = mark_colour(gs_name)
        for ax, (lib, df) in zip(axes, focus.items()):
            lib_label = (lib.replace('GO_Biological_Process_2023','GO Biological Process')
                            .replace('KEGG_2021_Human','KEGG')
                            .replace('Reactome_2022','Reactome'))
            barh_enrichment(ax, df, lib_label, n=12, colour=col)

        fig.suptitle(f'Pathway Enrichment — {gs_name}',
                     fontsize=13, fontweight='bold', y=1.02)
        plt.tight_layout()
        fname = f'pathway_{gs_name}.png'
        plt.savefig(fname, dpi=300, bbox_inches='tight')
        files.download(fname)
        plt.close()
        print(f"✔ Downloaded: {fname}")

    # ════════════════════════════════════════════════════════════════════════
    # Figure 3 — GO BP dot plot: gene sets × top pathways
    # ════════════════════════════════════════════════════════════════════════
    GO_LIB = 'GO_Biological_Process_2023'
    plot_rows = []
    for gs_name, libs in all_results.items():
        if GO_LIB not in libs or libs[GO_LIB].empty:
            continue
        for _, row in libs[GO_LIB].head(10).iterrows():
            plot_rows.append({
                'Gene_Set':      gs_name,
                'Term':          row['Term'][:65],
                'Adj_Pval':      row['Adj_Pvalue'],
                'Combined_Score': row['Combined_Score'],
                'Gene_Count':    len(str(row['Genes']).split(';')),
            })

    if plot_rows:
        dot_df = pd.DataFrame(plot_rows)
        # Keep top 20 terms by max score
        top_terms = (dot_df.groupby('Term')['Combined_Score']
                           .max().nlargest(20).index.tolist())
        dot_df = dot_df[dot_df['Term'].isin(top_terms)].copy()

        gs_order   = dot_df['Gene_Set'].unique().tolist()
        term_order = (dot_df.groupby('Term')['Combined_Score']
                            .max().sort_values(ascending=False).index.tolist())

        dot_df['x'] = dot_df['Gene_Set'].apply(lambda v: gs_order.index(v))
        dot_df['y'] = dot_df['Term'].apply(lambda v: term_order.index(v))
        dot_df['neg_log_p'] = -np.log10(dot_df['Adj_Pval'].clip(lower=1e-50))

        fig, ax = plt.subplots(figsize=(max(8, len(gs_order)*2.2),
                                         max(8, len(term_order)*0.45)))
        sc = ax.scatter(dot_df['x'], dot_df['y'],
                        s=dot_df['Gene_Count'] * 18,
                        c=dot_df['neg_log_p'],
                        cmap='YlOrRd', alpha=0.85,
                        edgecolors='k', linewidth=0.6,
                        vmin=0, vmax=dot_df['neg_log_p'].quantile(0.95))
        cbar = plt.colorbar(sc, ax=ax, shrink=0.6)
        cbar.set_label('−log₁₀(Adj P-value)', fontsize=10, fontweight='bold')

        ax.set_xticks(range(len(gs_order)))
        ax.set_xticklabels(gs_order, rotation=35, ha='right', fontsize=9)
        ax.set_yticks(range(len(term_order)))
        ax.set_yticklabels(term_order, fontsize=8)
        ax.set_title('GO Biological Process — Enrichment Dot Plot',
                     fontsize=13, fontweight='bold', pad=12)
        ax.grid(alpha=0.25)

        # Size legend
        for cnt in [5, 10, 20]:
            ax.scatter([], [], s=cnt*18, c='#aaa', edgecolors='k',
                       linewidth=0.6, label=f'{cnt} genes')
        ax.legend(title='Gene Count', loc='lower right',
                  fontsize=8, title_fontsize=9, framealpha=0.9)

        plt.tight_layout()
        fname = 'enrichment_GOBP_dotplot.png'
        plt.savefig(fname, dpi=300, bbox_inches='tight')
        files.download(fname)
        plt.close()
        print(f"✔ Downloaded: {fname}")

    # ════════════════════════════════════════════════════════════════════════
    # Figure 4 — Heatmap: top pathways × gene sets (GO BP + KEGG + Reactome)
    # ════════════════════════════════════════════════════════════════════════
    heat_rows = []
    for gs_name, libs in all_results.items():
        for lib in ['GO_Biological_Process_2023', 'KEGG_2021_Human', 'Reactome_2022']:
            if lib not in libs or libs[lib].empty:
                continue
            for _, row in libs[lib].head(8).iterrows():
                heat_rows.append({
                    'Gene_Set': gs_name,
                    'Pathway':  row['Term'][:55],
                    'Score':    -np.log10(float(row['Adj_Pvalue']) + 1e-50),
                })

    if heat_rows:
        heat_df  = pd.DataFrame(heat_rows)
        pivot    = heat_df.pivot_table(index='Pathway', columns='Gene_Set',
                                       values='Score', fill_value=0)
        # Sort rows by row-max
        pivot = pivot.loc[pivot.max(axis=1).sort_values(ascending=False).index]

        fig, ax = plt.subplots(figsize=(max(8, len(pivot.columns)*1.8),
                                         max(6, len(pivot)*0.38)))
        sns.heatmap(pivot, cmap='RdYlBu_r', linewidths=0.4, linecolor='#ddd',
                    ax=ax, cbar_kws={'label': '−log₁₀(Adj P-value)', 'shrink': 0.7})
        ax.set_xlabel('Gene Set', fontsize=11, fontweight='bold')
        ax.set_ylabel('Pathway',  fontsize=11, fontweight='bold')
        ax.set_title('Pathway Enrichment Heatmap — Top Pathways Across Gene Sets',
                     fontsize=13, fontweight='bold', pad=12)
        plt.xticks(rotation=35, ha='right', fontsize=9)
        plt.yticks(fontsize=8)
        plt.tight_layout()
        fname = 'enrichment_heatmap.png'
        plt.savefig(fname, dpi=300, bbox_inches='tight')
        files.download(fname)
        plt.close()
        print(f"✔ Downloaded: {fname}")

    # ════════════════════════════════════════════════════════════════════════
    # Figure 5 — MSigDB Hallmark summary (one bar per gene set)
    # ════════════════════════════════════════════════════════════════════════
    HALL_LIB = 'MSigDB_Hallmark_2020'
    hall_gs = {gs: libs[HALL_LIB]
               for gs, libs in all_results.items()
               if HALL_LIB in libs and not libs[HALL_LIB].empty}

    if hall_gs:
        n = len(hall_gs)
        ncols = min(2, n)
        nrows = (n + ncols - 1) // ncols
        fig, axes = plt.subplots(nrows, ncols,
                                  figsize=(9 * ncols, 6 * nrows))
        axes = np.array(axes).flatten()

        for ax, (gs_name, df) in zip(axes, hall_gs.items()):
            col = mark_colour(gs_name)
            barh_enrichment(ax, df, gs_name, n=10, colour=col)

        # Hide unused axes
        for ax in axes[len(hall_gs):]:
            ax.axis('off')

        fig.suptitle('MSigDB Hallmark Enrichment — BRD4-Dependent Gene Sets',
                     fontsize=13, fontweight='bold', y=1.01)
        plt.tight_layout()
        fname = 'enrichment_hallmark.png'
        plt.savefig(fname, dpi=300, bbox_inches='tight')
        files.download(fname)
        plt.close()
        print(f"✔ Downloaded: {fname}")

# ════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 80)
print("PATHWAY ENRICHMENT ANALYSIS COMPLETE")
print("=" * 80)
print("\nOutput files downloaded:")
print("  Pathway_Enrichment_Results.xlsx      — all enrichment tables")
print("  enrichment_all_genesets.png          — combined bar chart")
print("  pathway_<MARK>_<Category>.png        — per gene-set bar charts")
print("  enrichment_GOBP_dotplot.png          — GO BP dot plot")
print("  enrichment_heatmap.png               — pathway × gene-set heatmap")
print("  enrichment_hallmark.png              — MSigDB Hallmark bars")
print("=" * 80)
