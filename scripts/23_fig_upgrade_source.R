# 出图升级源数据（预注册变更记录"出图升级"①–⑥；只报告）
if (nzchar(Sys.getenv("R_LIBS_GLUT"))) .libPaths(Sys.getenv("R_LIBS_GLUT"))
suppressPackageStartupMessages({ library(data.table); library(Matrix); library(e1071); library(glmnet); library(xgboost); library(Boruta)
  library(metafor); library(kernelshap); library(clusterProfiler); library(org.Hs.eg.db) })
PROJ <- Sys.getenv("GLUT_PROJ", "."); RES <- file.path(Sys.getenv("GLUT_OUT", PROJ), "results/23_upgrade"); dir.create(RES, showWarnings = FALSE)
STIFF <- Sys.getenv("GLUT_BULK", "data/bulk"); G10 <- file.path(Sys.getenv("GLUT_OUT", PROJ), "results/10_gates")
write_tab <- function(df, name) { p <- file.path(RES, name); tmp <- paste0(p, ".tmp"); fwrite(as.data.frame(df), tmp, sep = "\t"); file.rename(tmp, p); message("written: ", p) }
grg <- readLines(file.path(PROJ, "geneset/GRG_full.txt")); hub <- readLines(file.path(G10, "02_hub.txt"))
per <- fread(file.path(STIFF, "R1rep_scores_per_sample.tsv")); per[cohort == "GSE135251_NASH", cohort := "GSE135251"]
ld <- function(cid) { L <- readRDS(file.path(STIFF, sprintf("%s_log2cpm.rds", cid))); P <- per[cohort == cid & !is.na(fibrosis) & gsm %in% colnames(L)]
  list(L = L[, P$gsm], y = factor(ifelse(P$fibrosis >= 2, "F2plus", "F0_1"), levels = c("F0_1", "F2plus"))) }
Tr <- ld("GSE135251"); E1 <- ld("GSE162694"); E2 <- ld("GSE130970")

## ④ SVM-RFE 曲线（复现脚本 10 的 ML，同种子）
SEED <- 20261112
fu <- fread(file.path(G10, "01_GRG_funnel.tsv")); cand <- fu[candidate == TRUE, gene]
X <- scale(t(Tr$L[cand, ])); y <- Tr$y
svmrfe_rank <- function(X, y) { s <- colnames(X); r <- character(0)
  while (length(s) > 1) { f <- svm(X[, s, drop = FALSE], y, kernel = "linear", cost = 1, scale = FALSE)
    w <- as.numeric((t(f$coefs) %*% f$SV)^2); worst <- s[which.min(w)]; r <- c(worst, r); s <- setdiff(s, worst) }
  c(s, r) }
set.seed(SEED); cvf <- cv.glmnet(X, y, family = "binomial", alpha = 1, nfolds = 10)   # 与脚本 10 相同的随机数消耗顺序
set.seed(SEED); rk <- svmrfe_rank(X, y); folds <- sample(rep(1:10, length.out = nrow(X)))
cv_err <- sapply(seq_along(rk), function(k) mean(sapply(1:10, function(f) { tr <- folds != f
  m <- svm(X[tr, rk[1:k], drop = FALSE], y[tr], kernel = "linear", cost = 1, scale = FALSE); mean(predict(m, X[!tr, rk[1:k], drop = FALSE]) != y[!tr]) })))
vote <- fread(file.path(G10, "02_ML_vote.tsv")); stopifnot(setequal(rk[1:which.min(cv_err)], vote[SVM_RFE == TRUE, gene]))
write_tab(data.table(k = seq_along(cv_err), cv_error = cv_err, added_gene = rk), "U4_svmrfe_curve.tsv")
## LASSO 路径
co <- as.matrix(coef(cvf$glmnet.fit)); lp <- data.table(lambda = cvf$glmnet.fit$lambda, t(co[-1, ]))
write_tab(melt(lp, id.vars = "lambda", variable.name = "gene", value.name = "coef"), "U4_lasso_path.tsv")
write_tab(data.table(lambda = cvf$lambda, cvm = cvf$cvm, cvsd = cvf$cvsd, lambda_1se = cvf$lambda.1se, lambda_min = cvf$lambda.min), "U4_lasso_cv.tsv")

## ⑤ Hedges g + DL meta
hg <- function(x1, x0) { n1 <- length(x1); n0 <- length(x0); sp <- sqrt(((n1 - 1) * var(x1) + (n0 - 1) * var(x0)) / (n1 + n0 - 2))
  g <- (mean(x1) - mean(x0)) / sp * (1 - 3 / (4 * (n1 + n0) - 9)); c(g = g, se = sqrt((n1 + n0) / (n1 * n0) + g^2 / (2 * (n1 + n0))), n1 = n1, n0 = n0) }
fr <- rbindlist(lapply(list(GSE135251 = Tr, GSE162694 = E1, GSE130970 = E2), function(d) rbindlist(lapply(hub, function(g) {
  h <- hg(d$L[g, d$y == "F2plus"], d$L[g, d$y == "F0_1"]); data.table(gene = g, g = h[["g"]], se = h[["se"]], n_F2 = h[["n1"]], n_F01 = h[["n0"]]) }))), idcol = "cohort")
mt <- rbindlist(lapply(hub, function(gg) { d <- fr[gene == gg]; stopifnot(nrow(d) == 3); m <- rma(yi = d$g, sei = d$se, method = "DL")
  data.table(gene = gg, g = as.numeric(m$b), se = m$se, lo = m$ci.lb, hi = m$ci.ub, P = m$pval, I2 = m$I2, tau2 = m$tau2) }))
write_tab(fr, "U5_forest_percohort.tsv"); write_tab(mt, "U5_forest_meta.tsv")

## ⑥ kernelSHAP（GLM，logit 尺度）
d <- data.frame(y = y, scale(t(Tr$L[hub, ]))); m <- glm(y ~ ., d, family = binomial)
Xs <- d[, hub]; set.seed(SEED); bg <- Xs[sample(nrow(Xs), 50), ]
ks <- kernelshap(m, Xs, bg_X = bg, pred_fun = function(o, x) predict(o, x, type = "link"), verbose = FALSE)
sh <- melt(data.table(id = seq_len(nrow(Xs)), ks$S), id.vars = "id", variable.name = "gene", value.name = "shap")
sh <- merge(sh, melt(data.table(id = seq_len(nrow(Xs)), Xs), id.vars = "id", variable.name = "gene", value.name = "z"), by = c("id", "gene"))
write_tab(sh, "U6_shap.tsv")

## ③ GO BP 过表达
ego <- enrichGO(gene = intersect(grg, rownames(Tr$L)), universe = rownames(Tr$L), OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
                pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 10, maxGSSize = 500)
eg <- as.data.table(as.data.frame(ego)); setorder(eg, p.adjust); write_tab(eg, "U3_GRG_GO_BP.tsv")
dg <- fread(file.path(G10, "01_DEG_all.tsv")); write_tab(dg[gene %in% grg, .(gene, logFC, adj.P.Val)], "U3_GRG_logFC.tsv")

## ①② 肝细胞：UMAP + 逐细胞 GRG 评分 + 亚状态 Ro/e
suppressPackageStartupMessages(library(Seurat))
hep <- readRDS("D:/cc/生信代码/Paper_Lactylation_MASH/sc_data/hep_states_tmp.rds")
um <- Embeddings(hep, "umap"); X <- GetAssayData(hep, layer = "data"); g_in <- intersect(grg, rownames(X))
cs <- data.table(cell = colnames(hep), umap_1 = um[, 1], umap_2 = um[, 2], hepstate = hep$hepstate, group = hep$group, sample = hep$sample,
                 GRG = Matrix::colMeans(X[g_in, ]))
for (g in c("GCDH", "AASS", "ACADSB")) cs[[g]] <- as.numeric(X[g, ])
rm(hep, X); invisible(gc())
write_tab(cs, "U1_hep_cells.tsv.gz")
tb <- as.matrix(table(cs$hepstate, cs$group)); ex <- outer(rowSums(tb), colSums(tb)) / sum(tb)
roe <- as.data.table(as.table(tb / ex)); setnames(roe, c("hepstate", "group", "Roe")); roe[, n := as.vector(tb)]
write_tab(roe, "U2_hepstate_Roe.tsv")
cat("DONE\n")
