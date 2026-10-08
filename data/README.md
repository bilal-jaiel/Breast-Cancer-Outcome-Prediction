# Data

The analysis uses `mrr_bio.Rdata`, the dataset provided for the M.R.R. course (Université Paris-Saclay, 2025). It is not redistributed in this repository.

To run the pipeline, place the file here:

```
data/raw/mrr_bio.Rdata
```

## Contents of `mrr_bio.Rdata`

| Object | Dimensions | Description |
|---|---|---|
| `clinical_data` | 1,231 × 24 | Clinical variables per patient (`S4Vectors::DataFrame`): age, tumour stage (AJCC), histology, prior treatment, sites of involvement, follow-up, `vital_status`, … |
| `GeneX` | 1,231 × 5,000 | mRNA expression of the 5,000 most variable genes, on the log2(counts + 1) scale |

Patient identifiers follow the TCGA barcode format (`TCGA-XX-XXXX-01A`), i.e. the cohort comes from The Cancer Genome Atlas breast cancer project (TCGA-BRCA).

Outcome (`vital_status`): 1,029 *Alive*, 201 *Dead*, 1 missing.

## Generated files

`notebooks/01_data_cleaning.Rmd` writes four objects to `data/processed/` (ignored by git), which `notebooks/02_modeling.Rmd` reads:

| File | Content |
|---|---|
| `cleaned_clinical_data_df.RData` | Cleaned and imputed clinical table (1,187 patients) |
| `df_clinical_data_brut.RData` | Clinical table before imputation (baseline model) |
| `GeneX.RData` | Expression matrix aligned with the cleaned patients |
| `GeneX_pca.RData` | Expression matrix used for the exploratory PCA/CCA |
