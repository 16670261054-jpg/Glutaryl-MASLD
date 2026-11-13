# 第二批（预注册变更记录"第二批判据"）：B 空间分带 + A4 沿既有拟时序的评分
if (nzchar(Sys.getenv("R_LIBS_GLUT"))) .libPaths(Sys.getenv("R_LIBS_GLUT"))
suppressPackageStartupMessages({ library(data.table); library(Seurat); library(Matrix); library(nnls); library(org.Hs.eg.db); library(AnnotationDbi) })
PROJ <- Sys.getenv("GLUT_PROJ", "."); RES <- file.path(Sys.getenv("GLUT_OUT", PROJ), "results/16_batch2"); dir.create(RES, showWarnings = FALSE)
GS <- file.path(PROJ, "geneset"); LAC <- Sys.getenv("GLUT_SN", "data/snRNA")
VDIR <- file.path(LAC, "01_数据/GSE192741_visium"); SEED <- 20261117; set.seed(SEED)
write_tab <- function(df, name) { p <- file.path(RES, name); tmp <- paste0(p, ".tmp")
  fwrite(as.data.frame(df), tmp, sep = "\t", quote = FALSE); file.rename(tmp, p); message("written: ", p) }
grg <- readLines(file.path(GS, "GRG_full.txt")); enz <- readLines(file.path(GS, "GRG_enzyme_axis.txt"))
mito <- unique(na.omit(suppressMessages(AnnotationDbi::select(org.Hs.eg.db, "GO:0005739", "SYMBOL", "GOALL"))$SYMBOL))
mito_ref <- sample(sort(setdiff(mito, grg)), 200)
writeLines(mito_ref, file.path(RES, "mito_reference_200.txt"))
SETS <- list(GRG_full = grg, GRG_enzyme_axis = enz, MITO_ref200 = mito_ref)
GENES <- c("GCDH", "AASS", "ACADSB", "SARDH")

## ================= B 空间 =================
ann <- fread(file.path(VDIR, "annot_humanVisium.csv")); ann[, barcode := sub("_[0-9]+$", "", spot)]
smap <- data.table(sample = c("JBO14", "JBO15", "JBO18", "JBO19", "JBO22"), gsm = paste0("GSM576442", 4:8), donor = c("H35", "H35", "H36", "H37", "H38"))
ref <- readRDS(file.path(LAC, "sc_data/deconv_ref_pseudobulk.rds"))
lab9 <- ifelse(grepl("^Hep[0-9]$", ref$label), "Hepatocyte", ref$label); sl9 <- paste(ref$sample, lab9, sep = "|")
agg9 <- t(rowsum(t(ref$mat), sl9)); ncl9 <- tapply(ref$ncell, sl9, sum)[colnames(agg9)]
cpm9 <- sweep(agg9, 2, pmax(colSums(agg9), 1), "/") * 1e6; l9 <- sub("^[^|]+\\|", "", colnames(agg9)); types <- sort(unique(l9))
prof9 <- sapply(types, function(L) rowMeans(cpm9[, l9 == L & ncl9 >= 20, drop = FALSE]))
sig9 <- unique(fread(file.path(LAC, "results/deconv_signature_genes_9type.tsv"))$gene)
sc1 <- function(dat, g) { g <- intersect(g, rownames(dat)); if (length(g) < 1) return(rep(NA_real_, ncol(dat))); Matrix::colMeans(dat[g, , drop = FALSE]) }
sp <- rbindlist(lapply(seq_len(nrow(smap)), function(i) { s <- smap$sample[i]
  m <- Read10X_h5(list.files(VDIR, pattern = paste0(smap$gsm[i], "_filtered.*h5$"), full.names = TRUE))
  a <- ann[sample == s]; hit <- intersect(a$barcode, colnames(m)); stopifnot(length(hit) / nrow(a) > 0.95)
  m <- m[, hit]; a <- a[match(hit, barcode)]
  ln <- m; ln@x <- log1p(ln@x / rep.int(colSums(m), diff(m@p)) * 1e4)
  cp <- m; cp@x <- cp@x / rep.int(colSums(m), diff(m@p)) * 1e6; cg <- intersect(sig9, rownames(m))
  fr <- t(apply(as.matrix(cp[cg, ]), 2, function(y) { co <- nnls(prof9[cg, types], y)$x; co / max(sum(co), 1e-12) })); colnames(fr) <- types
  d <- data.table(spot = a$spot, sample = s, donor = smap$donor[i], diet = a$diet, zonation = a$zonation, zonationGroup = a$zonationGroup, hep_frac = fr[, "Hepatocyte"])
  for (k in names(SETS)) d[[k]] <- sc1(ln, SETS[[k]])
  for (g in GENES) d[[g]] <- if (g %in% rownames(ln)) as.numeric(ln[g, ]) else NA_real_
  d }))
fwrite(sp, file.path(RES, "B_spatial_spot_scores.tsv.gz"), sep = "\t")
feats <- c(names(SETS), GENES)
zt <- rbindlist(lapply(feats, function(v) sp[, { ok <- is.finite(get(v)) & is.finite(zonation) & is.finite(hep_frac)
  ct <- suppressWarnings(cor.test(get(v)[ok], zonation[ok], method = "spearman", exact = FALSE))
  rx <- resid(lm(rank(get(v)[ok]) ~ rank(hep_frac[ok]))); ry <- resid(lm(rank(zonation[ok]) ~ rank(hep_frac[ok])))
  .(feature = v, diet = diet[1], n = sum(ok), rho = unname(ct$estimate), p = ct$p.value, rho_partial = cor(rx, ry)) }, by = sample]))
write_tab(zt, "B_zonation_per_section.tsv")
write_tab(zt[, .(n_sections = .N, n_central = sum(rho > 0), median_rho = median(rho), median_rho_partial = median(rho_partial)), by = feature], "B_zonation_summary.tsv")
med <- melt(sp[, lapply(.SD, median, na.rm = TRUE), by = .(sample, donor, diet), .SDcols = feats], id.vars = c("sample", "donor", "diet"))
dd <- med[, .(fatty_mean = mean(value[grepl("fat|NAFLD|obese", diet, ignore.case = TRUE)]), healthy_mean = mean(value[!grepl("fat|NAFLD|obese", diet, ignore.case = TRUE)]),
  n_fatty = sum(grepl("fat|NAFLD|obese", diet, ignore.case = TRUE)), n_healthy = sum(!grepl("fat|NAFLD|obese", diet, ignore.case = TRUE))), by = variable][, diff := fatty_mean - healthy_mean]
write_tab(med, "B_section_medians.tsv"); write_tab(dd, "B_fatty_vs_healthy.tsv")
print(unique(sp[, .(sample, diet)])); print(fread(file.path(RES, "B_zonation_summary.tsv"))); print(dd)

