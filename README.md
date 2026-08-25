# spatial-cell-classification

Machine learning pipeline for cell detection and classification on microscopy images, with spatial transcriptomics analysis. Built for any spatial biology dataset.

## Overview

The pipeline detects and classifies cells in microscopy images and combines the results with spatial transcriptomics analysis. It is designed to be dataset-agnostic: the same workflow applies to different tissues, stains, and imaging modalities.

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

Source images and expression data are not included; they belong to the originating research projects and cannot be shared publicly. The pipeline is demonstrated with structure only.

## License

MIT
