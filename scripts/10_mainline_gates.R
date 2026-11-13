# 主线 + 止损门禁 G1–G3（预注册 00_预注册_主流程_20261112.md D1–D3、D5 部分、E 节）
#   D1 DEG（GSE135251，limma-trend，F≥2 vs F0–1，BH FDR<0.05）∩ GRG_full → G1（≥5）
#   D2 四 ML（LASSO 1se / SVM-RFE / XGBoost gain 前10 / Boruta Confirmed）→ hub（4/4；<3 则 ≥3/4）
#   D3 GLM 冻结系数 → GSE162694、GSE130970 外部 AUC（DeLong CI）→ G2（≥1 个 CI 下限>0.5）
#   G3 hub 在 GSE289173 肝细胞伪 bulk 中 voom t 与训练 logFC 同号的比例 ≥ 0.5
if (nzchar(Sys.getenv("R_LIBS_GLUT"))) .libPaths(Sys.getenv("R_LIBS_GLUT"))
suppressPackageStartupMessages({ library(data.table); library(limma); library(edgeR); library(glmnet); library(e1071)
  library(xgboost); library(Boruta); library(pROC) })
PROJ <- Sys.getenv("GLUT_PROJ", "."); RES <- file.path(Sys.getenv("GLUT_OUT", PROJ), "results", "10_gates")
dir.create(RES, showWarnings = FALSE, recursive = TRUE)
STIFF <- Sys.getenv("GLUT_BULK", "data/bulk"); LAC <- Sys.getenv("GLUT_SN", "data/snRNA")
SM <- Sys.getenv("GLUT_SM", "data/GSE289173_series_matrix.txt.gz")
SEED <- 20261112
write_tab <- function(df, name) { p <- file.path(RES, name); tmp <- paste0(p, ".tmp")
  fwrite(as.data.frame(df), tmp, sep = "\t", quote = FALSE); file.rename(tmp, p); message("written: ", p) }
gate <- list()

## ---- 数据 ----
per <- fread(file.path(STIFF, "R1rep_scores_per_sample.tsv")); per[cohort == "GSE135251_NASH", cohort := "GSE135251"]
load_c <- function(cid) { L <- readRDS(file.path(STIFF, sprintf("%s_log2cpm.rds", cid)))
  P <- per[cohort == cid & !is.na(fibrosis) & gsm %in% colnames(L)]
  list(L = L[, P$gsm], y = factor(ifelse(P$fibrosis >= 2, "F2plus", "F0_1"), levels = c("F0_1", "F2plus")), P = P) }
Tr <- load_c("GSE135251"); E1 <- load_c("GSE162694"); E2 <- load_c("GSE130970")
write_tab(data.table(cohort = c("GSE135251", "GSE162694", "GSE130970"),
  n = c(length(Tr$y), length(E1$y), length(E2$y)),
  n_F2plus = c(sum(Tr$y == "F2plus"), sum(E1$y == "F2plus"), sum(E2$y == "F2plus"))), "00_cohorts.tsv")
grg <- readLines(file.path(PROJ, "geneset/GRG_full.txt"))

## ---- D1 DEG ----
des <- model.matrix(~ Tr$y)
fit <- eBayes(lmFit(Tr$L, des), trend = TRUE)
tt <- as.data.table(topTable(fit, coef = 2, n = Inf, sort.by = "none"), keep.rownames = "gene")
write_tab(tt[order(P.Value)], "01_DEG_all.tsv")
deg <- tt[adj.P.Val < 0.05]
common <- Reduce(intersect, list(rownames(Tr$L), rownames(E1$L), rownames(E2$L)))
cand_all <- intersect(deg$gene, grg); cand <- intersect(cand_all, common)
write_tab(data.table(gene = grg, in_train = grg %in% rownames(Tr$L), in_all3 = grg %in% common,
  DEG = grg %in% deg$gene)[, candidate := gene %in% cand][], "01_GRG_funnel.tsv")
gate$G1 <- data.table(gate = "G1", value = length(cand), rule = ">=5", PASS = length(cand) >= 5,
  note = sprintf("DEG=%d; DEG∩GRG=%d; 三队列共有=%d", nrow(deg), length(cand_all), length(cand)))
if (!gate$G1$PASS) { write_tab(rbindlist(gate, fill = TRUE), "99_gates.tsv"); stop("G1 未过，停") }

## ---- D2 四 ML ----
zs <- function(L, g) { X <- t(L[g, , drop = FALSE]); scale(X) }
X <- zs(Tr$L, cand); y <- Tr$y
set.seed(SEED); cvf <- cv.glmnet(X, y, family = "binomial", alpha = 1, nfolds = 10)
co <- as.matrix(coef(cvf, s = "lambda.1se")); lasso <- setdiff(rownames(co)[co[, 1] != 0], "(Intercept)")
svmrfe_rank <- function(X, y) { s <- colnames(X); r <- character(0)
  while (length(s) > 1) { f <- svm(X[, s, drop = FALSE], y, kernel = "linear", cost = 1, scale = FALSE)
    w <- as.numeric((t(f$coefs) %*% f$SV)^2); worst <- s[which.min(w)]; r <- c(worst, r); s <- setdiff(s, worst) }
  c(s, r) }
set.seed(SEED); rk <- svmrfe_rank(X, y); folds <- sample(rep(1:10, length.out = nrow(X)))
cv_err <- sapply(seq_along(rk), function(k) mean(sapply(1:10, function(f) { tr <- folds != f
  m <- svm(X[tr, rk[1:k], drop = FALSE], y[tr], kernel = "linear", cost = 1, scale = FALSE)
  mean(predict(m, X[!tr, rk[1:k], drop = FALSE]) != y[!tr]) })))
svmsel <- rk[1:which.min(cv_err)]
set.seed(SEED); bst <- xgb.train(params = list(objective = "binary:logistic", max_depth = 3, learning_rate = 0.1, nthread = 1),
  data = xgb.DMatrix(X, label = as.integer(y) - 1), nrounds = 100, verbose = 0)
imp <- as.data.table(xgb.importance(model = bst)); xgbsel <- head(imp[order(-Gain), Feature], 10)
set.seed(SEED); bo <- Boruta(X, y, maxRuns = 500); borsel <- names(bo$finalDecision)[bo$finalDecision == "Confirmed"]
vote <- data.table(gene = colnames(X))[, `:=`(LASSO = gene %in% lasso, SVM_RFE = gene %in% svmsel,
  XGBoost = gene %in% xgbsel, Boruta = gene %in% borsel)][, votes := LASSO + SVM_RFE + XGBoost + Boruta][order(-votes)]
vote <- merge(vote, tt[, .(gene, logFC, adj.P.Val)], by = "gene")[order(-votes, adj.P.Val)]
h4 <- vote[votes == 4, gene]; h3 <- vote[votes >= 3, gene]
hub <- if (length(h4) >= 3) h4 else h3; rule <- if (length(h4) >= 3) "4of4" else ">=3of4"
write_tab(vote, "02_ML_vote.tsv"); writeLines(hub, file.path(RES, "02_hub.txt"))
write_tab(data.table(metric = c("candidates", "LASSO", "SVM_RFE", "XGBoost_top10", "Boruta", "hub_4of4", "hub_3of4", "rule", "hub"),
  value = c(length(cand), length(lasso), length(svmsel), length(xgbsel), length(borsel), length(h4), length(h3), rule, paste(hub, collapse = ","))), "02_ML_summary.tsv")
if (length(hub) < 1) { write_tab(rbindlist(gate, fill = TRUE), "99_gates.tsv"); stop("hub 为空，停") }

## ---- D3 GLM + 外部 AUC ----
d <- data.frame(y = y, X[, hub, drop = FALSE]); m <- glm(y ~ ., d, family = binomial)
write_tab(as.data.table(summary(m)$coefficients, keep.rownames = "term"), "03_glm_coef.tsv")
set.seed(SEED); f10 <- sample(rep(1:10, length.out = nrow(d))); cvp <- numeric(nrow(d))
for (k in 1:10) { mk <- glm(y ~ ., d[f10 != k, ], family = binomial); cvp[f10 == k] <- predict(mk, d[f10 == k, ], type = "response") }
ev <- function(lab, yy, p) { r <- roc(yy, p, levels = c("F0_1", "F2plus"), direction = "<", quiet = TRUE); ci <- ci.auc(r, method = "delong")
  data.table(set = lab, n = length(yy), n_F2plus = sum(yy == "F2plus"), AUC = as.numeric(auc(r)), CI_low = ci[1], CI_high = ci[3]) }
ext <- function(E, lab) { Xe <- zs(E$L, hub); ev(lab, E$y, predict(m, data.frame(Xe), type = "response")) }
auc_tab <- rbind(ev("GSE135251_train_10foldCV", y, cvp), ext(E1, "GSE162694_external"), ext(E2, "GSE130970_external"))
write_tab(auc_tab, "03_AUC.tsv")
ok2 <- auc_tab[grepl("external", set), sum(CI_low > 0.5)]
gate$G2 <- data.table(gate = "G2", value = ok2, rule = ">=1 external CI_low>0.5", PASS = ok2 >= 1,
  note = paste(sprintf("%s %.3f[%.3f-%.3f]", auc_tab$set, auc_tab$AUC, auc_tab$CI_low, auc_tab$CI_high), collapse = "; "))

## ---- G3 snRNA 肝细胞同向 ----
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
dsn <- model.matrix(~ factor(sex) + factor(batch) + severity, md)
dge <- DGEList(cn); dge <- dge[filterByExpr(dge, dsn), , keep.lib.sizes = FALSE]; dge <- calcNormFactors(dge)
tsn <- topTable(eBayes(lmFit(voom(dge, dsn), dsn)), coef = "severity", n = Inf, sort.by = "none")
g3 <- data.table(gene = hub, train_logFC = tt[match(hub, gene), logFC], sn_t = tsn[hub, "t"], sn_P = tsn[hub, "P.Value"])
g3[, same_dir := sign(train_logFC) == sign(sn_t)]
write_tab(g3, "04_snRNA_hep_direction.tsv")
frac <- g3[!is.na(sn_t), mean(same_dir)]
gate$G3 <- data.table(gate = "G3", value = frac, rule = ">=0.5", PASS = !is.na(frac) & frac >= 0.5,
  note = sprintf("hub=%d; snRNA 可测=%d", length(hub), sum(!is.na(g3$sn_t))))
gt <- rbindlist(gate, fill = TRUE); write_tab(gt, "99_gates.tsv"); print(gt); print(auc_tab); print(vote[votes >= 3]); print(g3)
