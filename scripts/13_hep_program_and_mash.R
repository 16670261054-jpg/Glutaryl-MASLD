# 预注册变更记录 2026-11-12：S-hep（肝细胞内 GRG 程序 + N1'/N3' 零分布）与 E-MASH（冻结 4 基因模型用于 MASH vs 正常，只报告）
if (nzchar(Sys.getenv("R_LIBS_GLUT"))) .libPaths(Sys.getenv("R_LIBS_GLUT"))
suppressPackageStartupMessages({ library(data.table); library(limma); library(edgeR); library(pROC); library(org.Hs.eg.db); library(AnnotationDbi) })
PROJ <- Sys.getenv("GLUT_PROJ", "."); RES <- file.path(Sys.getenv("GLUT_OUT", PROJ), "results/13_hep_mash"); dir.create(RES, showWarnings = FALSE)
GS <- file.path(PROJ, "geneset"); G10 <- file.path(Sys.getenv("GLUT_OUT", PROJ), "results/10_gates")
LAC <- Sys.getenv("GLUT_SN", "data/snRNA"); STIFF <- Sys.getenv("GLUT_BULK", "data/bulk")
SM <- Sys.getenv("GLUT_SM", "data/GSE289173_series_matrix.txt.gz")
SEED <- 20261114; N <- 1000; set.seed(SEED)
write_tab <- function(df, name) { p <- file.path(RES, name); tmp <- paste0(p, ".tmp")
  fwrite(as.data.frame(df), tmp, sep = "\t", quote = FALSE); file.rename(tmp, p); message("written: ", p) }
rl <- function(f) readLines(file.path(GS, f))
sets <- list(GRG_full = rl("GRG_full.txt"), GRG_substrate = rl("GRG_substrate.txt"),
             GRG_substrate_both = rl("GRG_substrate_both.txt"), GRG_enzyme_axis = rl("GRG_enzyme_axis.txt"))

## ---- S-hep ----
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
dge <- DGEList(cn); dge <- dge[filterByExpr(dge, des), , keep.lib.sizes = FALSE]; dge <- calcNormFactors(dge)
v <- voom(dge, des); tt <- topTable(eBayes(lmFit(v, des)), coef = "severity", n = Inf, sort.by = "none")
tvec <- setNames(tt$t, rownames(tt)); mu <- rowMeans(v$E)
idx <- lapply(sets, function(g) which(rownames(v) %in% g))
cam <- as.data.table(camera(v, idx, des, contrast = ncol(des)), keep.rownames = "set")
mito <- unique(na.omit(suppressMessages(AnnotationDbi::select(org.Hs.eg.db, "GO:0005739", "SYMBOL", "GOALL"))$SYMBOL))
grg_all <- unique(unlist(sets))
pools <- list(N1 = setdiff(rownames(v), grg_all), N3 = setdiff(intersect(rownames(v), mito), grg_all)); nb <- c(N1 = 20, N3 = 5)
mk_draw <- function(pool, nbins, g) { br <- unique(quantile(mu[pool], 0:nbins / nbins))
  b <- findInterval(mu[pool], br, rightmost.closed = TRUE, all.inside = TRUE); names(b) <- pool
  hb <- findInterval(mu[g], br, rightmost.closed = TRUE, all.inside = TRUE); tb <- table(hb)
  function() unlist(lapply(names(tb), function(x) sample(names(b)[b == as.integer(x)], tb[[x]], replace = sum(b == as.integer(x)) < tb[[x]]))) }
DRAWS <- list()
nul <- rbindlist(lapply(names(sets), function(s) { g <- rownames(v)[idx[[s]]]; T0 <- mean(tvec[g])
  rbindlist(lapply(names(pools), function(p) { dr <- mk_draw(pools[[p]], nb[[p]], g); Tn <- replicate(N, mean(tvec[dr()])); DRAWS[[paste(s, p)]] <<- data.table(set = s, null = p, mean_t = Tn)
    data.table(set = s, null = p, pool = length(pools[[p]]), mean_t = T0, null_median = median(Tn), null_p = (sum(abs(Tn) >= abs(T0)) + 1) / (N + 1)) })) }))
S <- merge(cam, dcast(nul, set + mean_t ~ null, value.var = c("null_median", "null_p")), by = "set")
S[, `:=`(n_donors = ncol(cn), n_tested = sapply(idx, length)[set])]
S[, claim := fifelse(PValue < 0.05 & null_p_N1 < 0.05 & null_p_N3 < 0.05, "program_beyond_general_mito",
             fifelse(PValue < 0.05 & null_p_N1 < 0.05, "program_mainly_mito_enzymes", "not_supported"))]
write_tab(S, "13_Shep_hepatocyte_program.tsv")
write_tab(data.table(gene = rownames(tt), logFC = tt$logFC, t = tt$t, P = tt$P.Value, BH = p.adjust(tt$P.Value, "BH"))[gene %in% grg_all][order(P)], "13_Shep_GRG_genes.tsv")

## ---- E-MASH（冻结 F≥2 模型；只报告） ----
per <- fread(file.path(STIFF, "R1rep_scores_per_sample.tsv")); per[cohort == "GSE135251_NASH", cohort := "GSE135251"]
L0 <- readRDS(file.path(STIFF, "GSE135251_log2cpm.rds")); P0 <- per[cohort == "GSE135251" & !is.na(fibrosis) & gsm %in% colnames(L0)]
hub <- readLines(file.path(G10, "02_hub.txt"))
zs <- function(L, g) scale(t(L[g, , drop = FALSE]))
fit_m <- function(g) glm(y ~ ., data.frame(y = as.integer(P0$fibrosis >= 2), zs(L0[, P0$gsm], g)), family = binomial)
sm135 <- readLines(gzfile(file.path(LAC, "01_数据/GSE135251/GSE135251_series_matrix.txt.gz")))
gsm135 <- gsub('"', '', strsplit(grep("^!Sample_geo_accession", sm135, value = TRUE), "\t")[[1]][-1])
grp135 <- gsub('"', '', strsplit(grep("group in paper", sm135, value = TRUE), "\t")[[1]][-1]); grp135 <- trimws(sub(".*:", "", grp135))
e126 <- fread(file.path(LAC, "results/GSE126848_expr_log2cpm.tsv")); L126 <- as.matrix(e126[, -1]); rownames(L126) <- e126[[1]]
p126 <- fread(file.path(LAC, "results/GSE126848_pheno.tsv"))
L162 <- readRDS(file.path(STIFF, "GSE162694_log2cpm.rds")); p162 <- per[cohort == "GSE162694" & gsm %in% colnames(L162)]
comps <- list(
  GSE126848_NASH_vs_healthy = list(L = L126, ids = p126[group %in% c("MASH", "Normal-weight"), gsm], y = p126[group %in% c("MASH", "Normal-weight"), as.integer(group == "MASH")]),
  GSE126848_NAFL_vs_healthy = list(L = L126, ids = p126[group %in% c("NAFL", "Normal-weight"), gsm], y = p126[group %in% c("NAFL", "Normal-weight"), as.integer(group == "NAFL")]),
  GSE162694_NAFLD_vs_normal = list(L = L162, ids = p162$gsm, y = as.integer(!p162$normal)),
  GSE135251_NASH_vs_control = list(L = L0, ids = gsm135[grepl("^NASH|^control", grp135)], y = as.integer(grepl("^NASH", grp135[grepl("^NASH|^control", grp135)]))))
common <- Reduce(intersect, c(list(rownames(L0)), lapply(comps, function(x) rownames(x$L))))
mu0 <- rowMeans(L0[common, P0$gsm]); br <- quantile(mu0, 0:20 / 20); b0 <- findInterval(mu0, br, rightmost.closed = TRUE, all.inside = TRUE); names(b0) <- common
hb <- b0[hub]
fauc <- function(y, p) as.numeric(auc(roc(y, p, levels = c(0, 1), direction = "<", quiet = TRUE)))
eval_g <- function(g, m) sapply(comps, function(x) { Lx <- x$L[, x$ids]; fauc(x$y, predict(m, data.frame(zs(Lx, g)))) })
m_hub <- fit_m(hub); obs <- eval_g(hub, m_hub)
ci <- sapply(comps, function(x) { Lx <- x$L[, x$ids]; r <- roc(x$y, predict(m_hub, data.frame(zs(Lx, hub))), levels = c(0, 1), direction = "<", quiet = TRUE); as.numeric(ci.auc(r, method = "delong"))[c(1, 3)] })
nulA <- t(replicate(N, { g <- sapply(hb, function(x) sample(setdiff(names(b0)[b0 == x], hub), 1)); eval_g(g, fit_m(g)) }))
E <- data.table(comparison = names(comps), n = sapply(comps, function(x) length(x$y)), n_case = sapply(comps, function(x) sum(x$y)),
  AUC = obs, CI_low = ci[1, ], CI_high = ci[2, ], null_median = apply(nulA, 2, median), percentile = colMeans(sweep(nulA, 2, obs, "<")) * 100)
write_tab(E, "13_EMASH_frozen_model.tsv")
## 出图源数据（2026-11-13 追加；不改任何统计）
write_tab(rbindlist(DRAWS), "13_Shep_null_draws.tsv")
write_tab(data.table(set = names(sets), mean_t = sapply(names(sets), function(s) mean(tvec[rownames(v)[idx[[s]]]]))), "13_Shep_obs_meant.tsv")
pred <- rbindlist(lapply(names(comps), function(k) { x <- comps[[k]]; Lx <- x$L[, x$ids]
  data.table(comparison = k, sample = x$ids, y = x$y, score = as.numeric(predict(m_hub, data.frame(zs(Lx, hub))))) }))
write_tab(pred, "13_EMASH_predictions.tsv")
print(S); print(E)
