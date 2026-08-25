# spatial-cell-classification

Machine learning pipeline for cell detection and classification on confocal microscopy images, with spatial transcriptomics analysis. Built as a course research project at the Department of Immunotechnology, Lund University (Medicon Village, 2025).

## Overview

The pipeline classifies cells in confocal microscopy images as part of a spatial scoring system in ovarian cancer research. It combines image-based cell detection with spatial transcriptomics analysis to support diagnostic-tool development.

## Repository structure

```
cell_detection/     cell detection and classification (Python, SquidPy, QuPath)
spatial/            spatial transcriptomics analysis (scanpy)
models/             trained model artifacts (not committed)
```

## Requirements

- Python 3.9+
- scikit-learn, pandas, NumPy, matplotlib, SHAP
- SquidPy, QuPath (for image export)
- scanpy (for spatial transcriptomics)

## Usage

```bash
# 1. Export cell annotations from QuPath
#    (see cell_detection/qupath_export.md)

# 2. Train the classifier
python cell_detection/train_classifier.py --data <annotations.csv>

# 3. Spatial analysis
python spatial/run_scanpy.py --input <expression_matrix.h5ad>
```

## Data

Patient-derived microscopy images and expression data are not included; they are confidential and cannot be shared publicly. The pipeline is demonstrated with structure only.

## License

MIT
