# PyWGCNA with Affinity Propagation

An alternative module detection step for WGCNA: **affinity propagation** replaces hierarchical clustering of the TOM matrix. Applied to the 5xFAD Alzheimer's mouse model (GEO: GSE168137) using the PyWGCNA wrapper. See the paper PDF for full details.

## Dataset

Bulk RNA-seq of hippocampus and cortex from transgenic and wildtype 5xFAD mice at 4, 8, 12 and 18 months (192 samples).

## Pipeline

1. Pre-process with PyWGCNA (outlier removal, keep highly expressed genes: 2569 remain, log2-CPM, quantile normalisation)
2. Build the network (soft threshold beta = 15, TOM)
3. Detect modules with affinity propagation (preference = -5, damping = 0.5)
4. Correlate module eigengenes with traits (age, genotype, tissue, sex)
5. Plot module heatmaps, GO enrichment and hub genes
6. Venn diagrams of tissue vs age genes and hub genes

## Results

- 11 modules (20 to 407 genes)
- Age and genotype showed the strongest correlations; tissue and sex were weaker
- 567 shared tissue/age genes and 90 shared hub genes

## Run

```bash
pip install -r requirements.txt
jupyter notebook notebooks/analysis.ipynb
```
