# spatial-cell-classification

Machine learning pipeline for cell detection and classification on microscopy images, with spatial transcriptomics analysis. Designed to be dataset-agnostic: the same workflow applies to any tissue, stain, or imaging modality.

## Overview

The pipeline detects and classifies cells in microscopy images and combines the results with spatial transcriptomics analysis. It is built around QuPath for image annotation and export, Python for classification, and scanpy for spatial analysis.

## Repository structure

```
scripts/            QuPath Groovy scripts (cell detection, pixel classification)
cell_detection/     cell detection and classification (Python, SquidPy, QuPath)
spatial/            spatial transcriptomics analysis (scanpy)
survival_analysis/  survival curve generation from cell abundance data
```

## Requirements

- QuPath (for image annotation and export)
- Python 3.9+
- scikit-learn, pandas, NumPy, matplotlib, SHAP
- SquidPy (for spatial analysis)
- scanpy (for spatial transcriptomics)
- lifelines (for survival analysis)

## Usage

```bash
# 1. Run cell detection in QuPath
#    QuPath -> Automate -> Run script: scripts/pipe.groovy

# 2. Export cell annotations
#    (see cell_detection/qupath_export.md)

# 3. Train the classifier
python cell_detection/train_classifier.py --data <annotations.csv>

# 4. Spatial analysis
python spatial/run_scanpy.py --input <expression_matrix.h5ad>

# 5. Survival analysis
python survival_analysis/generate_km_curves.py --input <cell_abundance.csv>
```

## Data

Source images, QuPath projects, and expression data are not included; they belong to the originating research projects and cannot be shared publicly. The pipeline is demonstrated with structure only.

## License

MIT
