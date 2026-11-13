# 补充图第二批源数据（只报告）：GRG×细胞类型矩阵、列线图与训练集校准（rms）、蛋白组全蛋白 limma（Kleiner 模型，同 21 预处理）
if (nzchar(Sys.getenv("R_LIBS_GLUT"))) .libPaths(Sys.getenv("R_LIBS_GLUT"))
suppressPackageStartupMessages({ library(data.table); library(edgeR); library(rms); library(limma) })
PROJ <- Sys.getenv("GLUT_PROJ", "."); RES <- file.path(Sys.getenv("GLUT_OUT", PROJ), "results/24_supp2"); dir.create(RES, showWarnings = FALSE)
STIFF <- Sys.getenv("GLUT_BULK", "data/bulk"); LAC <- Sys.getenv("GLUT_SN", "data/snRNA")
write_tab <- function(df, name) { p <- file.path(RES, name); tmp <- paste0(p, ".tmp"); fwrite(as.data.frame(df), tmp, sep = "\t"); file.rename(tmp, p); message("written: ", p) }
grg <- readLines(file.path(PROJ, "geneset/GRG_full.txt")); hub <- readLines(file.path(Sys.getenv("GLUT_OUT", PROJ), "results/10_gates/02_hub.txt"))

## 1. GRG × 细胞类型（供者×类型伪 bulk logCPM，类型内供者均值，≥20 核）
ref <- readRDS(file.path(LAC, "sc_data/deconv_ref_pseudobulk.rds"))
lab <- ifelse(grepl("^Hep[0-9]+$", ref$label), "Hepatocyte", ref$label); key <- paste(ref$sample, lab, sep = "|")
agg <- t(rowsum(t(ref$mat), key)); nc <- tapply(ref$ncell, key, sum)[colnames(agg)]
lc <- cpm(DGEList(agg[, nc >= 20]), log = TRUE); ct <- sub("^[^|]+\\|", "", colnames(lc))
m <- sapply(split(seq_along(ct), ct), function(i) rowMeans(lc[intersect(grg, rownames(lc)), i, drop = FALSE]))
write_tab(data.table(gene = rownames(m), m), "S8_GRG_celltype_logCPM.tsv")

## 2. 列线图（训练集 GSE135251，4 hub 的 z 值，lrm）+ bootstrap 校准（B = 200）
per <- fread(file.path(STIFF, "R1rep_scores_per_sample.tsv")); per[cohort == "GSE135251_NASH", cohort := "GSE135251"]
L <- readRDS(file.path(STIFF, "GSE135251_log2cpm.rds")); P <- per[cohort == "GSE135251" & !is.na(fibrosis) & gsm %in% colnames(L)]
d <- data.frame(y = as.integer(P$fibrosis >= 2), scale(t(L[hub, P$gsm]))); dd <- datadist(d); options(datadist = "dd")
f <- lrm(as.formula(paste("y ~", paste(hub, collapse = " + "))), data = d, x = TRUE, y = TRUE)
set.seed(20261121); cal <- calibrate(f, method = "boot", B = 200)
saveRDS(list(fit = f, cal = cal, dd = dd), file.path(RES, "S7_nomogram_objects.rds"))
v <- validate(f, method = "boot", B = 200)
write_tab(data.table(index = rownames(v), v[, c("index.orig", "optimism", "index.corrected")]), "S7_validate_boot.tsv")
write_tab(data.table(pred = cal[, "predy"], apparent = cal[, "calibrated.orig"], bias_corrected = cal[, "calibrated.corrected"]), "S7_calibration_boot.tsv")

## 3. 蛋白组全蛋白 limma（Kleiner 纤维化模型，预处理与 21 相同）
DL <- file.path(PROJ, "data/commmed2025")
md <- fread(file.path(DL, "43856_2025_780_MOESM3_ESM.txt"), na.strings = "NA")[!is.na(liver_proteomics_filename)]
pq <- fread(file.path(DL, "43856_2025_780_MOESM4_ESM.txt"), na.strings = c("NA", ""))
X <- log2(as.matrix(pq[, md$liver_proteomics_filename, with = FALSE])); g1 <- sub(";.*", "", pq$Genes); nd <- rowSums(!is.na(X))
k <- !is.na(g1) & g1 != "" & nd >= 0.7 * ncol(X); X <- X[k, ]; g1 <- g1[k]; nd <- nd[k]
o <- order(g1, -nd); X <- X[o, ]; g1 <- g1[o]; dp <- !duplicated(g1); X <- X[dp, ]; rownames(X) <- g1[dp]
X <- sweep(X, 2, apply(X, 2, median, na.rm = TRUE)) + median(X, na.rm = TRUE)
md[, `:=`(batch = factor(substr(liver_proteomics_filename, 1, 8)), sex = factor(gender), fib = as.numeric(sub("F", "", kleiner_fibrosis_grade)))]
des <- model.matrix(~ batch + sex + fib, md); tt <- topTable(eBayes(lmFit(X, des)), coef = "fib", n = Inf, sort.by = "none")
chk <- fread(file.path(Sys.getenv("GLUT_OUT", PROJ), "results/21_protein/02_key_proteins.tsv"))[grepl("Kleiner", model)]
stopifnot(max(abs(tt[chk$gene, "t"] - chk$t)) < 1e-8)
write_tab(data.table(gene = rownames(tt), logFC = tt$logFC, t = tt$t, P = tt$P.Value, BH = tt$adj.P.Val, GRG = rownames(tt) %in% grg), "S11_protein_all_kleiner.tsv")
cat("DONE\n")
