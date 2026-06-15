# Peak Overlap Analysis Script

## Overview

This script analyzes peak overlaps and creates Venn diagrams to visualize the effect of BRD4 inhibition on TGFβ-induced histone modifications. It compares chromatin immunoprecipitation-sequencing (ChIP-seq) peaks between TGFβ and TGFβ+BRD4i conditions across four histone marks:

- **H3K4me3** - Active promoter marks
- **H3K4me1** - Active enhancer marks
- **H3K27ac** - Acetylated chromatin (active regions)
- **H3K122ac** - H3 histone acetylation mark

## Features

- **Peak Overlap Analysis**: Identifies peaks that are:
  - Maintained between conditions
  - Lost with BRD4 inhibition
  - Gained with BRD4 inhibition

- **Visualizations**:
  - Venn diagrams for each histone mark showing peak overlap
  - Bar charts with absolute peak counts
  - Percentage-based comparison charts

- **Summary Statistics**: Generates detailed metrics for each histone mark including:
  - Total peak counts per condition
  - Number and percentage of overlapping peaks
  - Net changes in peak distribution

## Requirements

### Dependencies

```bash
pip install pandas matplotlib seaborn matplotlib-venn numpy
```

### Input Files

The script expects narrowPeak files in the working directory. By default, it looks for:

```
TGFB-TB1-H3K4me3_narrow_qval_peaks.narrowPeak
TGF-IT1-H3K4me3_narrow_qval_peaks.narrowPeak
TGFB-TB2-H3K4me1_narrow_qval_peaks.narrowPeak
TGF-IT2-H3K4me1_narrow_qval_peaks.narrowPeak
TGFB-TB3-H3K27ac_narrow_qval_peaks.narrowPeak
TGF-IT3-H3K27ac_narrow_qval_peaks.narrowPeak
TGFB-TB4-H3K122ac_narrow_qval_peaks.narrowPeak
TGF-IT4-H3K122ac_narrow_qval_peaks.narrowPeak
```

**narrowPeak Format** (tab-separated):
```
chrom | chromStart | chromEnd | name | score | strand | signalValue | pValue | qValue | peak
```

## Usage

### Basic Run

```bash
python scripts/peak_overlap_analysis.py
```

### Customizing File Paths

Edit the `file_paths` dictionary in the script to point to your narrowPeak files:

```python
file_paths = {
    'TB1_H3K4me3': '/path/to/your/TGFB-TB1-H3K4me3_narrow_qval_peaks.narrowPeak',
    'IT1_H3K4me3': '/path/to/your/TGF-IT1-H3K4me3_narrow_qval_peaks.narrowPeak',
    # ... etc
}
```

## Output Files

The script generates four output files in the working directory:

1. **peak_overlap_venn_diagrams.png** - High-resolution Venn diagrams (PNG format)
2. **peak_overlap_venn_diagrams.pdf** - Vector format Venn diagrams (PDF format)
3. **peak_categories_comparison.png** - Bar charts comparing peak categories (PNG)
4. **peak_categories_comparison.pdf** - Bar charts comparing peak categories (PDF)

### Output Details

#### Venn Diagrams
- Shows 2×2 grid layout (one for each histone mark)
- Left circle: TGFβ-only peaks
- Right circle: TGFβ+BRD4i-only peaks
- Center: Maintained peaks (present in both conditions)
- Bottom annotation: Net gain or loss with BRD4 inhibition

#### Comparison Charts
- **Left panel**: Absolute peak counts for each category
- **Right panel**: Percentage distribution of peak categories

## Script Functions

### `load_narrowpeak(filepath)`
Loads a narrowPeak file and creates peak IDs based on genomic coordinates.
- Bins peaks into 500bp windows for robust overlap detection
- Returns pandas DataFrame with peak information

### `analyze_peak_overlaps(tb_peaks, it_peaks, mark_name)`
Compares peaks between two conditions.
- Returns dictionary with overlap statistics
- Calculates percentages and net changes

## Interpretation Guide

### Peak Categories

- **Maintained Peaks**: Present in both TGFβ and TGFβ+BRD4i conditions
  - Indicates BRD4-independent regulation
  
- **Lost Peaks**: Present in TGFβ only, absent with BRD4i
  - Suggests BRD4-dependent peak establishment/maintenance
  
- **Gained Peaks**: Present only with BRD4i treatment
  - May indicate compensatory or alternative chromatin remodeling

### Net Change Direction

- **Net GAIN**: More peaks gained than lost with BRD4 inhibition
- **Net LOSS**: More peaks lost than gained with BRD4 inhibition

## Example Output

```
Analyzing peak overlaps between TGFβ and TGFβ+BRD4i conditions...
================================================================================

H3K4me3:
  TGFβ only: 12,345 peaks
  TGFβ + BRD4i: 11,234 peaks
  Maintained: 10,000 peaks (81.0% of TB)
  Lost with BRD4i: 2,345 peaks
  Gained with BRD4i: 1,234 peaks

...
```

## Troubleshooting

### FileNotFoundError
- **Issue**: Script cannot find narrowPeak files
- **Solution**: Verify file paths in `file_paths` dictionary match your actual files. Check working directory.

### Missing Dependencies
- **Issue**: `ImportError` for pandas, matplotlib, etc.
- **Solution**: Install required packages: `pip install -r requirements.txt` (create requirements.txt with dependencies listed)

### No Output Generated
- **Issue**: Script runs but produces no output files
- **Solution**: Check file permissions in working directory; ensure matplotlib can write PNG/PDF files

## Future Enhancements

- Add statistical significance testing for overlaps
- Include motif enrichment analysis for gained/lost peaks
- Support batch processing of multiple datasets
- Interactive visualization options (Plotly)

## Citation

If you use this script in your research, please cite:
- The original narrowPeak format papers
- matplotlib_venn library

## License

[Specify your license here]

## Contact

For questions or issues, please contact the repository maintainer.
