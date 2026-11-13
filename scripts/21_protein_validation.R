# P1 人肝蛋白组验证（预注册变更记录 P1，看数据分布前写死）
if (nzchar(Sys.getenv("R_LIBS_GLUT"))) .libPaths(Sys.getenv("R_LIBS_GLUT"))
suppressPackageStartupMessages({ library(data.table); library(limma); library(org.Hs.eg.db); library(AnnotationDbi) })
PROJ <- Sys.getenv("GLUT_PROJ", "."); RES <- file.path(Sys.getenv("GLUT_OUT", PROJ), "results/21_protein"); dir.create(RES, showWarnings = FALSE)
DL <- file.path(PROJ, "data/commmed2025"); SEED <- 20261120; N <- 1000
write_tab <- function(df, name) { p <- file.path(RES, name); tmp <- paste0(p, ".tmp"); fwrite(as.data.frame(df), tmp, sep = "\t"); file.rename(tmp, p); message("written: ", p) }
grg <- readLines(file.path(PROJ, "geneset/GRG_full.txt")); KEY <- c("GCDH", "AASS", "ACADSB", "SARDH", "ACOT13", "ACAA2")

md <- fread(file.path(DL, "43856_2025_780_MOESM3_ESM.txt"), na.strings = "NA")[!is.na(liver_proteomics_filename)]
pq <- fread(file.path(DL, "43856_2025_780_MOESM4_ESM.txt"), na.strings = c("NA", ""))
X <- as.matrix(pq[, md$liver_proteomics_filename, with = FALSE]); X <- log2(X)
gene1 <- sub(";.*", "", pq$Genes); nd <- rowSums(!is.na(X))
keep <- !is.na(gene1) & gene1 != "" & nd >= 0.7 * ncol(X); X <- X[keep, ]; gene1 <- gene1[keep]; nd <- nd[keep]
o <- order(gene1, -nd); X <- X[o, ]; gene1 <- gene1[o]; d <- !duplicated(gene1); X <- X[d, ]; rownames(X) <- gene1[d]
X <- sweep(X, 2, apply(X, 2, median, na.rm = TRUE)) + median(X, na.rm = TRUE)
md[, `:=`(batch = factor(substr(liver_proteomics_filename, 1, 8)), sex = factor(gender), fib = as.numeric(sub("F", "", kleiner_fibrosis_grade)))]
stopifnot(!anyNA(md$fib))
write_tab(data.table(n_samples = ncol(X), n_proteins = nrow(X), n_GRG = sum(rownames(X) %in% grg), key_present = paste(intersect(KEY, rownames(X)), collapse = ",")), "00_input.tsv")

mito <- unique(na.omit(suppressMessages(AnnotationDbi::select(org.Hs.eg.db, "GO:0005739", "SYMBOL", "GOALL"))$SYMBOL))
mu <- rowMeans(X, na.rm = TRUE)
mk_draw <- function(pool, nb, g) { br <- unique(quantile(mu[pool], 0:nb / nb)); b <- findInterval(mu[pool], br, rightmost.closed = TRUE, all.inside = TRUE); names(b) <- pool
  hb <- findInterval(mu[g], br, rightmost.closed = TRUE, all.inside = TRUE); tb <- table(hb)
  function() unlist(lapply(names(tb), function(x) { cand <- names(b)[b == as.integer(x)]; sample(cand, tb[[x]], replace = length(cand) < tb[[x]]) })) }
run <- function(var, label) {
  des <- model.matrix(as.formula(paste("~ batch + sex +", var)), md)
  fit <- eBayes(lmFit(X, des)); tt <- topTable(fit, coef = var, n = Inf, sort.by = "none"); tv <- setNames(tt$t, rownames(tt))
  Xc <- X[rowSums(is.na(X)) == 0, ]; gi <- which(rownames(X) %in% grg); cm <- camera(Xc, list(GRG_full = which(rownames(Xc) %in% grg)), des, contrast = var)  # camera 不接受 NA：只用全检出蛋白
  g <- rownames(X)[gi]; T0 <- mean(tv[g], na.rm = TRUE); set.seed(SEED)
  pools <- list(N1 = setdiff(rownames(X), grg), N3 = setdiff(intersect(rownames(X), mito), grg)); nbs <- c(N1 = 20, N3 = 5)
  nul <- sapply(names(pools), function(p) { dr <- mk_draw(pools[[p]], nbs[[p]], g); Tn <- replicate(N, mean(tv[dr()], na.rm = TRUE)); c(median = median(Tn), P = (sum(abs(Tn) >= abs(T0)) + 1) / (N + 1)) })
  pc <- data.table(gene = "COL1A1", t = tt["COL1A1", "t"], P = tt["COL1A1", "P.Value"])
  kg <- data.table(model = label, gene = KEY, logFC = tt[KEY, "logFC"], t = tt[KEY, "t"], P = tt[KEY, "P.Value"])[, BH := p.adjust(P, "BH")]
  list(sum = data.table(model = label, n = ncol(X), n_GRG = length(gi), n_GRG_camera = sum(rownames(Xc) %in% grg), n_prot_camera = nrow(Xc), direction = cm$Direction, camera_P = cm$PValue, mean_t = T0,
       null_median_N1 = nul["median", "N1"], null_P_N1 = nul["P", "N1"], null_median_N3 = nul["median", "N3"], null_P_N3 = nul["P", "N3"],
       COL1A1_t = pc$t[1], COL1A1_P = pc$P[1]),
       key = kg, genes = data.table(model = label, gene = rownames(tt), logFC = tt$logFC, t = tt$t, P = tt$P.Value)[gene %in% grg]) }
R2 <- run("fib", "Kleiner_fibrosis")
S <- R2$sum
S[, posctrl_ok := COL1A1_t > 0 & COL1A1_P < 0.05]
S[, verdict := fifelse(!posctrl_ok, "VOID_posctrl", fifelse(camera_P < 0.05 & direction == "Down" & null_P_N1 < 0.05,
                 fifelse(null_P_N3 < 0.05, "protein_replicated_beyond_mito", "protein_replicated"), "not_replicated"))]
write_tab(S, "01_summary.tsv"); write_tab(R2$key, "02_key_proteins.tsv"); write_tab(R2$genes, "03_GRG_proteins.tsv")
donor <- data.table(md[, .(patient_name, fib, sex, batch)], GRG_score = colMeans(t(scale(t(X[rownames(X) %in% grg, ]))), na.rm = TRUE))
for (g in KEY) donor[[g]] <- X[g, ]
write_tab(donor, "04_sample_scores.tsv"); print(S); print(R2$key)
