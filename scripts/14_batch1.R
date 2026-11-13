# 第一批加分模块（预注册变更记录"第一批判据"）：A2 A3 A5 C1 C2 C3
if (nzchar(Sys.getenv("R_LIBS_GLUT"))) .libPaths(Sys.getenv("R_LIBS_GLUT"))
suppressPackageStartupMessages({ library(data.table); library(limma); library(edgeR); library(fgsea); library(pROC); library(org.Hs.eg.db); library(AnnotationDbi) })
PROJ <- Sys.getenv("GLUT_PROJ", "."); RES <- file.path(Sys.getenv("GLUT_OUT", PROJ), "results/14_batch1"); dir.create(RES, showWarnings = FALSE)
GS <- file.path(PROJ, "geneset"); G10 <- file.path(Sys.getenv("GLUT_OUT", PROJ), "results/10_gates")
LAC <- Sys.getenv("GLUT_SN", "data/snRNA"); STIFF <- Sys.getenv("GLUT_BULK", "data/bulk")
SM <- Sys.getenv("GLUT_SM", "data/GSE289173_series_matrix.txt.gz")
SEED <- 20261115; set.seed(SEED)
write_tab <- function(df, name) { p <- file.path(RES, name); tmp <- paste0(p, ".tmp")
  fwrite(as.data.frame(df), tmp, sep = "\t", quote = FALSE); file.rename(tmp, p); message("written: ", p) }
grg <- readLines(file.path(GS, "GRG_full.txt")); hub <- readLines(file.path(G10, "02_hub.txt"))
key_genes <- unique(c(hub, "GCDH", "AASS", "SIRT5", "KAT2A"))

## ---- snRNA 公共部分 ----
ref <- readRDS(file.path(LAC, "sc_data/deconv_ref_pseudobulk.rds"))
sm <- readLines(gzfile(SM)); ch <- sm[grepl("^!Sample_characteristics_ch1", sm)]
row_of <- function(pat) gsub('"', '', strsplit(grep(pat, sm, value = TRUE)[1], "\t")[[1]][-1])
getch <- function(k) { ln <- ch[grepl(k, ch)][1]; trimws(sub(".*:", "", gsub('"', '', strsplit(ln, "\t")[[1]][-1]))) }
meta <- data.table(sample = row_of("^!Sample_title"), group = getch("disease group"), sex = getch("Sex: "))
meta[, severity := c(ctrl = 0, MASL = 1, eMASH = 2, aMASH = 3)[group]]; meta[, batch := ifelse(grepl("^FL", sample), "FL", "ctrl_series")]
pb <- function(labels_keep, collapse_to = NULL) {
  lab <- ref$label; if (!is.null(collapse_to)) lab[lab %in% labels_keep] <- collapse_to
  key <- paste(ref$sample, lab, sep = "|"); agg <- t(rowsum(t(ref$mat), key)); nc <- tapply(ref$ncell, key, sum)[colnames(agg)]
  list(agg = agg, nc = nc) }
voom_ct <- function(agg, nc, ct, form = ~ factor(sex) + factor(batch) + severity, sub = NULL) {
  k <- sub("^[^|]+\\|", "", colnames(agg)) == ct
  cn <- round(agg[, k, drop = FALSE]); colnames(cn) <- sub("\\|.*", "", colnames(cn)); md <- meta[match(colnames(cn), sample)]
  md[, n_nuc := as.numeric(nc[k])]; keep <- md$n_nuc >= 20; if (!is.null(sub)) keep <- keep & eval(sub, md, parent.frame())
  cn <- cn[, keep, drop = FALSE]; md <- md[keep]
  if (ncol(cn) < 20) return(NULL)
  des <- model.matrix(form, md); dge <- DGEList(cn); dge <- dge[filterByExpr(dge, des), , keep.lib.sizes = FALSE]; dge <- calcNormFactors(dge)
  v <- voom(dge, des); list(v = v, des = des, md = md, fit = eBayes(lmFit(v, des))) }

## ---- A2 细胞定位 ----
P1 <- pb(grep("^Hep[0-9]+$", unique(ref$label), value = TRUE), "Hepatocyte")
cpm_all <- t(t(P1$agg) / colSums(P1$agg)) * 1e6
ct <- sub("^[^|]+\\|", "", colnames(cpm_all)); ok <- P1$nc >= 20
loc <- rbindlist(lapply(key_genes, function(g) { if (!g %in% rownames(cpm_all)) return(NULL)
  med <- tapply(cpm_all[g, ok], ct[ok], median)
  data.table(gene = g, celltype = names(med), median_CPM = as.numeric(med)) }))
loc[, hep_specificity := median_CPM[celltype == "Hepatocyte"] / max(median_CPM[celltype != "Hepatocyte"]), by = gene]
write_tab(loc, "A2_celltype_localization.tsv")

## ---- A3 肝细胞亚状态 ----
P2 <- pb(character(0))
states <- grep("^Hep[0-9]+$", unique(ref$label), value = TRUE)
a3 <- rbindlist(lapply(states, function(s) { r <- voom_ct(P2$agg, P2$nc, s); if (is.null(r)) return(data.table(state = s, n_donors = NA))
  idx <- list(GRG_full = which(rownames(r$v) %in% grg))
  cm <- as.data.table(camera(r$v, idx, r$des, contrast = ncol(r$des)), keep.rownames = "set")
  data.table(state = s, n_donors = ncol(r$v), n_genes = cm$NGenes, direction = cm$Direction, P = cm$PValue) }), fill = TRUE)
a3[!is.na(P), BH := p.adjust(P, "BH")]
hep_nc <- P2$nc[grepl("\\|Hep[0-9]+$", names(P2$nc))]
prop <- data.table(sample = sub("\\|.*", "", names(hep_nc)), state = sub("^[^|]+\\|", "", names(hep_nc)), n = as.numeric(hep_nc))
prop[, frac := n / sum(n), by = sample]; prop <- merge(prop, meta[, .(sample, severity)], by = "sample")
pr <- prop[, { ct <- suppressWarnings(cor.test(frac, severity, method = "spearman", exact = FALSE)); .(rho = unname(ct$estimate), P = ct$p.value) }, by = state]
pr[, BH := p.adjust(P, "BH")]
write_tab(a3, "A3_hepstate_GRG_camera.tsv"); write_tab(pr, "A3_hepstate_proportion_vs_severity.tsv"); write_tab(prop, "A3_hepstate_proportions.tsv")

## ---- A5 上游 TF ----
rh <- voom_ct(P1$agg, P1$nc, "Hepatocyte")
tt <- topTable(rh$fit, coef = "severity", n = Inf, sort.by = "none"); tv <- setNames(tt$t, rownames(tt))
tf <- readRDS(file.path(GS, "msigdb_C3_TFT_GTRD.rds")); tfl <- split(tf$gene_symbol, tf$gs_name)
set.seed(SEED); fg <- as.data.table(fgseaMultilevel(tfl, sort(tv, decreasing = TRUE), minSize = 15, maxSize = 500, eps = 0))
fg[, TF := sub("_TARGET_GENES$", "", pathway)]
fg[, `:=`(TF_t = tt[TF, "t"], TF_P = tt[TF, "P.Value"], overlap_GRG = sapply(leadingEdge, function(x) sum(x %in% grg)))]
fg[, candidate := padj < 0.05 & NES < 0 & !is.na(TF_t) & TF_t < 0 & TF_P < 0.05]
fg[, leadingEdge := sapply(leadingEdge, paste, collapse = ",")]
write_tab(fg[order(pval)], "A5_TF_GTRD_fgsea.tsv")

## ---- C1 病理分项 ----
zmean <- function(L, g) { g <- intersect(g, rownames(L)); colMeans(t(scale(t(L[g, , drop = FALSE])))) }
soft <- function(f) { ln <- readLines(f, warn = FALSE); gsm <- NA; rows <- list()
  for (x in ln) { if (grepl("^\\^SAMPLE = ", x)) { gsm <- sub("^\\^SAMPLE = ", "", x); rows[[gsm]] <- list(gsm = gsm) }
    else if (!is.na(gsm) && grepl("^!Sample_characteristics_ch1 = ", x)) { kv <- sub("^!Sample_characteristics_ch1 = ", "", x)
      rows[[gsm]][[trimws(sub(":.*", "", kv))]] <- trimws(sub("^[^:]*:", "", kv)) } }
  rbindlist(rows, fill = TRUE) }
p130 <- soft("D:/cc/MASH博士课题基因筛选/_tools/raw/cohort_scan/soft_GSE130970.txt")
L130 <- readRDS(file.path(STIFF, "GSE130970_log2cpm.rds")); p130 <- p130[gsm %in% colnames(L130)]
s130 <- zmean(L130[, p130$gsm], grg)
vars <- c(steatosis = "steatosis grade", inflammation = "lobular inflammation grade", ballooning = "cytological ballooning grade", fibrosis = "fibrosis stage", NAS = "nafld activity score")
c1 <- rbindlist(lapply(names(vars), function(v) { x <- as.numeric(p130[[vars[[v]]]]); ct <- cor.test(s130, x, method = "spearman", exact = FALSE)
  data.table(cohort = "GSE130970", variable = v, n = sum(!is.na(x)), rho = unname(ct$estimate), P = ct$p.value) }))
sm135 <- readLines(gzfile(file.path(LAC, "01_数据/GSE135251/GSE135251_series_matrix.txt.gz")))
g135 <- gsub('"', '', strsplit(grep("^!Sample_geo_accession", sm135, value = TRUE), "\t")[[1]][-1])
gv <- function(k) as.numeric(trimws(sub(".*:", "", gsub('"', '', strsplit(grep(paste0("^!Sample_characteristics_ch1.*", k), sm135, value = TRUE)[1], "\t")[[1]][-1]))))
L135 <- readRDS(file.path(STIFF, "GSE135251_log2cpm.rds")); d135 <- data.table(gsm = g135, NAS = gv("nas score"), fibrosis = gv("fibrosis stage"))[gsm %in% colnames(L135)]
s135 <- zmean(L135[, d135$gsm], grg)
c1 <- rbind(c1, rbindlist(lapply(c("NAS", "fibrosis"), function(v) { ct <- cor.test(s135, d135[[v]], method = "spearman", exact = FALSE)
  data.table(cohort = "GSE135251", variable = v, n = nrow(d135), rho = unname(ct$estimate), P = ct$p.value) })))
c1[, BH := p.adjust(P, "BH"), by = cohort]
write_tab(c1, "C1_pathology_components.tsv")
write_tab(rbind(data.table(cohort = "GSE130970", gsm = p130$gsm, score = s130), data.table(cohort = "GSE135251", gsm = d135$gsm, score = s135)), "C1_scores.tsv")

## ---- C2 疾病谱 GSE164760 ----
smx <- readLines(gzfile("D:/cc/生信分析/MASH数据集/BULK/GSE164760,BULK/GSE164760_series_matrix.txt.gz"))
i0 <- grep("^!series_matrix_table_begin", smx); i1 <- grep("^!series_matrix_table_end", smx)
ex <- fread(text = smx[(i0 + 1):(i1 - 1)]); pr_ <- ex[[1]]; ex <- as.matrix(ex[, -1]); rownames(ex) <- pr_
colnames(ex) <- gsub('"', '', strsplit(grep("^!Sample_geo_accession", smx, value = TRUE), "\t")[[1]][-1])
tis <- trimws(sub(".*tissue:", "", gsub('"', '', strsplit(grep("tissue:", smx, value = TRUE)[1], "\t")[[1]][-1])))
gpl <- fread("D:/cc/生信分析/MASH数据集/BULK/GSE164760,BULK/GPL13667-15572.txt", skip = 43, sep = "\t", quote = "")
ann <- gpl[, .(ID, symbol = `Gene Symbol`)][symbol != "" & !grepl("///", symbol)][ID %in% rownames(ex)]
exs <- log2(ex[ann$ID, ] + 1); o <- order(-rowMeans(exs)); exs <- exs[o, ]; sy <- ann$symbol[o]; kk <- !duplicated(sy); exs <- exs[kk, ]; rownames(exs) <- sy[kk]
d3 <- data.table(gsm = colnames(exs), tissue = tis, score = zmean(exs, grg))
d3[, ord := fcase(tissue == "Healthy liver", 0L, tissue == "NASH liver", 1L, tissue == "Non-tumoral NASH liver adjacent to HCC", 2L, tissue == "NASH-HCC tumor", 3L, default = NA_integer_)]
sp <- cor.test(d3[!is.na(ord), score], d3[!is.na(ord), ord], method = "spearman", exact = FALSE); kw <- kruskal.test(score ~ factor(ord), d3[!is.na(ord)])
write_tab(d3, "C2_GSE164760_scores.tsv")
write_tab(rbind(d3[!is.na(ord), .(n = .N, median = median(score)), by = .(ord, tissue)][order(ord)][, .(stat = paste0("median_", tissue), n, value = median)],
  data.table(stat = c("spearman_rho", "spearman_P", "KW_P"), n = sum(!is.na(d3$ord)), value = c(unname(sp$estimate), sp$p.value, kw$p.value))), "C2_GSE164760_summary.tsv")

## ---- C3 校准 / DCA（冻结模型） ----
per <- fread(file.path(STIFF, "R1rep_scores_per_sample.tsv")); per[cohort == "GSE135251_NASH", cohort := "GSE135251"]
P0 <- per[cohort == "GSE135251" & !is.na(fibrosis) & gsm %in% colnames(L135)]
zs <- function(L, g) scale(t(L[g, , drop = FALSE]))
m <- glm(y ~ ., data.frame(y = as.integer(P0$fibrosis >= 2), zs(L135[, P0$gsm], hub)), family = binomial)
L162 <- readRDS(file.path(STIFF, "GSE162694_log2cpm.rds")); p162 <- per[cohort == "GSE162694" & !is.na(fibrosis) & gsm %in% colnames(L162)]
e126 <- fread(file.path(LAC, "results/GSE126848_expr_log2cpm.tsv")); L126 <- as.matrix(e126[, -1]); rownames(L126) <- e126[[1]]
p126 <- fread(file.path(LAC, "results/GSE126848_pheno.tsv"))[group %in% c("MASH", "Normal-weight")]
sets <- list(GSE162694_F2plus = list(p = predict(m, data.frame(zs(L162[, p162$gsm], hub)), type = "response"), y = as.integer(p162$fibrosis >= 2)),
             GSE126848_MASH_vs_normal = list(p = predict(m, data.frame(zs(L126[, p126$gsm], hub)), type = "response"), y = as.integer(p126$group == "MASH")))
cal <- rbindlist(lapply(names(sets), function(k) { p <- pmin(pmax(sets[[k]]$p, 1e-6), 1 - 1e-6); y <- sets[[k]]$y; lp <- qlogis(p)
  sl <- coef(glm(y ~ lp, family = binomial))[2]; it <- coef(glm(y ~ offset(lp), family = binomial))[1]
  data.table(set = k, n = length(y), prevalence = mean(y), Brier = mean((p - y)^2), cal_intercept = unname(it), cal_slope = unname(sl)) }))
dca <- rbindlist(lapply(names(sets), function(k) { p <- sets[[k]]$p; y <- sets[[k]]$y; n <- length(y)
  rbindlist(lapply(seq(0.1, 0.6, 0.05), function(pt) { tp <- sum(p >= pt & y == 1); fp <- sum(p >= pt & y == 0)
    data.table(set = k, threshold = pt, NB_model = tp / n - fp / n * pt / (1 - pt), NB_all = mean(y) - (1 - mean(y)) * pt / (1 - pt), NB_none = 0) })) }))
write_tab(cal, "C3_calibration.tsv"); write_tab(dca, "C3_DCA.tsv")
print(loc[celltype == "Hepatocyte"]); print(a3); print(pr); print(head(fg[order(pval), .(TF, NES, padj, TF_t, TF_P, overlap_GRG, candidate)], 15)); print(c1); print(cal)
