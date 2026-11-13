# 出图源数据（只算出图要用、但前面脚本没落盘的量；每项与已落盘统计逐一核对）
if (nzchar(Sys.getenv("R_LIBS_GLUT"))) .libPaths(Sys.getenv("R_LIBS_GLUT"))
suppressPackageStartupMessages({ library(data.table); library(limma); library(edgeR); library(pROC); library(fgsea) })
PROJ <- Sys.getenv("GLUT_PROJ", "."); RES <- file.path(Sys.getenv("GLUT_OUT", PROJ), "results/20_figsrc"); dir.create(RES, showWarnings = FALSE)
GS <- file.path(PROJ, "geneset"); LAC <- Sys.getenv("GLUT_SN", "data/snRNA"); STIFF <- Sys.getenv("GLUT_BULK", "data/bulk")
SM <- Sys.getenv("GLUT_SM", "data/GSE289173_series_matrix.txt.gz")
write_tab <- function(df, name) { p <- file.path(RES, name); tmp <- paste0(p, ".tmp")
  fwrite(as.data.frame(df), tmp, sep = "\t", quote = FALSE); file.rename(tmp, p); message("written: ", p) }
grg <- readLines(file.path(GS, "GRG_full.txt")); hub <- readLines(file.path(Sys.getenv("GLUT_OUT", PROJ), "results/10_gates/02_hub.txt"))
chk <- list()

## 1. 肝细胞伪 bulk：每供者 GRG 评分与关键基因（与 13 相同的过滤与 voom）
ref <- readRDS(file.path(LAC, "sc_data/deconv_ref_pseudobulk.rds"))
lab <- ifelse(grepl("^Hep[0-9]+$", ref$label), "Hepatocyte", ref$label); key <- paste(ref$sample, lab, sep = "|")
agg <- t(rowsum(t(ref$mat), key)); nc <- tapply(ref$ncell, key, sum)[colnames(agg)]
sm <- readLines(gzfile(SM)); ch <- sm[grepl("^!Sample_characteristics_ch1", sm)]
row_of <- function(pat) gsub('"', '', strsplit(grep(pat, sm, value = TRUE)[1], "\t")[[1]][-1])
getch <- function(k) { ln <- ch[grepl(k, ch)][1]; trimws(sub(".*:", "", gsub('"', '', strsplit(ln, "\t")[[1]][-1]))) }
meta <- data.table(sample = row_of("^!Sample_title"), group = getch("disease group"), sex = getch("Sex: "))
meta[, severity := c(ctrl = 0, MASL = 1, eMASH = 2, aMASH = 3)[group]]; meta[, batch := ifelse(grepl("^FL", sample), "FL", "ctrl_series")]
k <- sub("^[^|]+\\|", "", colnames(agg)) == "Hepatocyte"
cn <- round(agg[, k]); colnames(cn) <- sub("\\|.*", "", colnames(cn)); md <- meta[match(colnames(cn), sample)]
md[, n_nuc := as.numeric(nc[k])]; keep <- md$n_nuc >= 20; cn <- cn[, keep]; md <- md[keep]
des <- model.matrix(~ factor(sex) + factor(batch) + severity, md)
dge <- DGEList(cn); dge <- dge[filterByExpr(dge, des), , keep.lib.sizes = FALSE]; dge <- calcNormFactors(dge); v <- voom(dge, des)
g_in <- intersect(grg, rownames(v)); Z <- t(scale(t(v$E[g_in, ])))
mito_ref <- readLines(file.path(Sys.getenv("GLUT_OUT", PROJ), "results/16_batch2/mito_reference_200.txt")); m_in <- intersect(mito_ref, rownames(v))
donor <- data.table(md[, .(sample, group, severity, sex, batch, n_nuc)], GRG_score = colMeans(Z), MITO_ref_score = colMeans(t(scale(t(v$E[m_in, ])))))
for (g in c(hub, "GCDH", "AASS", "SIRT5", "NR1H4")) if (g %in% rownames(v)) donor[[g]] <- v$E[g, ]
write_tab(donor, "F2_hep_donor_scores.tsv")
chk$n_GRG_hep <- data.table(item = "GRG genes tested in hepatocyte (vs 13: 171)", value = length(g_in))

## 2. NR1H4 GSEA 运行曲线（与 14 相同排序）
tt <- topTable(eBayes(lmFit(v, des)), coef = "severity", n = Inf, sort.by = "none"); st <- sort(setNames(tt$t, rownames(tt)), decreasing = TRUE)
tf <- readRDS(file.path(GS, "msigdb_C3_TFT_GTRD.rds")); fxr <- unique(tf$gene_symbol[tf$gs_name == "NR1H4_TARGET_GENES"])
ed <- plotEnrichmentData(fxr, st)
write_tab(as.data.table(ed$curve), "F5_NR1H4_curve.tsv"); write_tab(data.table(rank = ed$ticks$rank, stat = ed$ticks$stat, in_GRG = names(st)[ed$ticks$rank] %in% grg), "F5_NR1H4_ticks.tsv")
write_tab(data.table(rank = seq_along(st), stat = as.numeric(st)), "F5_rank_stats.tsv")
chk$es <- data.table(item = "NR1H4 ES min (fgsea A5 NES=-1.43; ES sign check)", value = min(ed$curve$ES))

## 3. F≥2 ROC 逐样本预测（复现 10 的 GLM 与 10 折 CV，核对 AUC）
per <- fread(file.path(STIFF, "R1rep_scores_per_sample.tsv")); per[cohort == "GSE135251_NASH", cohort := "GSE135251"]
ld <- function(cid) { L <- readRDS(file.path(STIFF, sprintf("%s_log2cpm.rds", cid))); P <- per[cohort == cid & !is.na(fibrosis) & gsm %in% colnames(L)]
  list(L = L[, P$gsm], y = factor(ifelse(P$fibrosis >= 2, "F2plus", "F0_1"), levels = c("F0_1", "F2plus")), id = P$gsm) }
Tr <- ld("GSE135251"); E1 <- ld("GSE162694"); E2 <- ld("GSE130970"); zs <- function(L, g) scale(t(L[g, , drop = FALSE]))
d <- data.frame(y = Tr$y, zs(Tr$L, hub)); m <- glm(y ~ ., d, family = binomial)
set.seed(20261112); f10 <- sample(rep(1:10, length.out = nrow(d))); cvp <- numeric(nrow(d))
for (kk in 1:10) { mk <- glm(y ~ ., d[f10 != kk, ], family = binomial); cvp[f10 == kk] <- predict(mk, d[f10 == kk, ], type = "response") }
pr <- rbind(data.table(set = "GSE135251_train_10foldCV", sample = Tr$id, y = as.integer(Tr$y == "F2plus"), score = cvp),
  data.table(set = "GSE162694_external", sample = E1$id, y = as.integer(E1$y == "F2plus"), score = predict(m, data.frame(zs(E1$L, hub)), type = "response")),
  data.table(set = "GSE130970_external", sample = E2$id, y = as.integer(E2$y == "F2plus"), score = predict(m, data.frame(zs(E2$L, hub)), type = "response")))
write_tab(pr, "F4_F2plus_predictions.tsv")
ref_auc <- fread(file.path(Sys.getenv("GLUT_OUT", PROJ), "results/10_gates/03_AUC.tsv"))
re <- pr[, .(AUC = as.numeric(auc(roc(y, score, levels = c(0, 1), direction = "<", quiet = TRUE)))), by = set]
chk$auc <- merge(re, ref_auc[, .(set, AUC_ref = AUC)], by = "set")[, .(item = paste("AUC", set), value = AUC - AUC_ref)]
rocs <- rbind(rbindlist(lapply(split(pr, pr$set), function(x) { r <- roc(x$y, x$score, levels = c(0, 1), direction = "<", quiet = TRUE)
  data.table(set = x$set[1], fpr = 1 - r$specificities, tpr = r$sensitivities) })),
  rbindlist(lapply(split(fread(file.path(Sys.getenv("GLUT_OUT", PROJ), "results/13_hep_mash/13_EMASH_predictions.tsv")), by = "comparison"), function(x) {
    r <- roc(x$y, x$score, levels = c(0, 1), direction = "<", quiet = TRUE); data.table(set = x$comparison[1], fpr = 1 - r$specificities, tpr = r$sensitivities) })))
write_tab(rocs, "F4_roc_curves.tsv")

## 4. Visium：两张代表切片的坐标 + 每切片按 zonation 十分位的中位评分
VD <- file.path(LAC, "01_数据/GSE192741_visium"); sp <- fread(file.path(Sys.getenv("GLUT_OUT", PROJ), "results/16_batch2/B_spatial_spot_scores.tsv.gz"))
gsm <- c(JBO14 = "GSM5764424", JBO15 = "GSM5764425", JBO18 = "GSM5764426", JBO19 = "GSM5764427", JBO22 = "GSM5764428")
pos <- rbindlist(lapply(names(gsm), function(s) { f <- list.files(VD, pattern = paste0(gsm[[s]], "_tissue_positions"), full.names = TRUE)
  p <- fread(f, header = FALSE, col.names = c("barcode", "in_tissue", "row", "col", "pxl_row", "pxl_col")); p[, sample := s] }))
sp[, barcode := sub("_[0-9]+$", "", spot)]; sp <- merge(sp, pos[, .(sample, barcode, pxl_row, pxl_col)], by = c("sample", "barcode"))
write_tab(sp[, .(sample, diet, barcode, pxl_row, pxl_col, zonation, GRG_full, MITO_ref200, GRG_enzyme_axis)], "F3_spatial_spots.tsv")
sp[, zdec := pmin(10L, as.integer(ceiling(10 * frank(zonation, ties.method = "average") / .N))), by = sample]  # zonation 在 0.99 处有大量并列，按平均秩分十分位
zb <- melt(sp, id.vars = c("sample", "diet", "zdec"), measure.vars = c("GRG_full", "MITO_ref200"))[, zonation := sp$zonation[rep(seq_len(nrow(sp)), 2)]][, .(median = median(value), n_spots = .N, zonation_median = median(zonation)), by = .(sample, diet, zdec, variable)]
zb[, z_within := (median - mean(median)) / sd(median), by = .(sample, variable)]
write_tab(zb, "F3_zonation_deciles.tsv")
CK <- rbindlist(chk, fill = TRUE); write_tab(CK, "00_checks.tsv"); print(CK)
