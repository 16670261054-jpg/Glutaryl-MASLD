# Glutaryl-MASLD

Analysis code for:

> Mi Y, Mo S. Hepatocellular decline of a glutarylation-related mitochondrial enzyme program in metabolic dysfunction-associated steatotic liver disease. (submitted to *Functional & Integrative Genomics*)

All analyses use publicly available data. The scripts reproduce every number, table and figure panel reported in the manuscript and its Online Resources.

## Contents

| Path | Purpose |
|---|---|
| `scripts/02_parse_glutarylome.py` | Extract glutaryl-lysine PSMs (Percolator q ≤ 0.01) from PRIDE pepXML files (PXD005895, PXD007881) |
| `scripts/03_build_geneset.py` | MGI one-to-one orthology mapping and assembly of the glutarylation-related gene (GRG) set |
| `scripts/10_mainline_gates.R` | Bulk differential expression (GSE135251), four-algorithm feature selection, frozen logistic model and external AUCs |
| `scripts/13_hep_program_and_mash.R` | Donor-level hepatocyte pseudobulk (voom–limma, camera), matched null distributions, MASLD-versus-normal AUCs |
| `scripts/14_batch1.R` | Cell-type localisation, hepatocyte states, GTRD transcription factor enrichment, bulk histology, GSE164760, calibration and DCA |
| `scripts/16_batch2.R` | Visium zonation analysis (GSE192741) |
| `scripts/18_fxr_agonist.R` | FXR-agonist datasets (GSE76161, GSE124053) |
| `scripts/21_protein_validation.R` | Liver proteome (PXD051911) association with Kleiner fibrosis stage and matched null distributions |
| `scripts/26_fxr_mito_null.R` | FXR-agonist datasets tested against mitochondria-matched random gene sets |
| `scripts/20_`, `22_`–`25_*.R` | Source data for figure panels |
| `scripts/fig*.R` | Figures 1–5 and Online Resource figures S1–S14 (read only the saved result tables) |
| `geneset/` | GRG lists used in the paper |

## Inputs

Set these environment variables before running (defaults in brackets):

| Variable | Content |
|---|---|
| `GLUT_PROJ` [`.`] | Repository root; downloaded files are expected under `data/` (PRIDE pepXML, MGI `HOM_MouseHumanSequence.rpt`, proteome supplementary tables, FXR GEO series matrices) |
| `GLUT_OUT` [`GLUT_PROJ`] | Output root for `results/` and `figures/` |
| `GLUT_BULK` [`data/bulk`] | Filtered, TMM-normalised log2 CPM matrices and phenotype tables for GSE135251, GSE162694 and GSE130970 (`<GSE>_log2cpm.rds`, `R1rep_scores_per_sample.tsv`) |
| `GLUT_SN` [`data/snRNA`] | Processed GSE289173 single-nucleus objects and donor × cell-type pseudobulk (`sc_data/`), GSE126848 matrices and GSE192741 Visium files, as produced for Mi et al. (2026) *Funct Integr Genomics* 26:305 |
| `GLUT_SM` | GSE289173 series matrix (donor metadata) |
| `geneset/msigdb_C3_TFT_GTRD.rds` | MSigDB C3 TFT:GTRD collection, e.g. `msigdbr::msigdbr(species = "Homo sapiens", collection = "C3", subcollection = "TFT:GTRD")` |

## Run order

```
python 02_parse_glutarylome.py data/PXD005895/glut.pep.xml geneset/wagner2017 Wagner2017_PXD005895 129.095
python 02_parse_glutarylome.py data/PXD007881/liver.pep.xml geneset/schmiesing2018 Schmiesing2018_PXD007881
python 03_build_geneset.py
Rscript 10 → 13 → 14 → 16 → 18 → 20 → 21 → 22 → 23 → 24 → 25 → 26
Rscript fig1.R … fig5.R, fig_supp.R, fig_supp2.R, fig_supp3.R
```

Random seeds are fixed in each script. Software: R 4.6.1 (limma 3.66, edgeR 4.8, Seurat 5.5, glmnet 5.0, e1071 1.7, xgboost 3.2, Boruta 10.0, pROC 1.19, metafor 5.0, fgsea 1.36, clusterProfiler 4.18, rms 8.1) and Python 3.11.

## License

MIT
