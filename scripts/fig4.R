# Fig. 4 冻结 4 基因模型的区分度（只读盘画图）
source(file.path(Sys.getenv("GLUT_PROJ", "."), "scripts/fig_config.R"))
rc <- fread(file.path(R, "20_figsrc/F4_roc_curves.tsv"))
au1 <- fread(file.path(R, "10_gates/03_AUC.tsv"))[, .(set, n, n_case = n_F2plus, AUC, CI_low, CI_high, task = "F\u22652 vs F0\u20131")]
au2 <- fread(file.path(R, "13_hep_mash/13_EMASH_frozen_model.tsv"))[, .(set = comparison, n, n_case, AUC, CI_low, CI_high, task = "Disease vs normal")]
au <- rbind(au1, au2)
NM <- c(GSE135251_train_10foldCV = "GSE135251 (training, 10-fold CV)", GSE162694_external = "GSE162694 (external)", GSE130970_external = "GSE130970 (external)",
        GSE126848_NASH_vs_healthy = "GSE126848 MASH vs normal weight", GSE126848_NAFL_vs_healthy = "GSE126848 MASL vs normal weight",
        GSE162694_NAFLD_vs_normal = "GSE162694 MASLD vs normal histology", GSE135251_NASH_vs_control = "GSE135251 MASH vs control")
au[, lab := sprintf("%s\nAUC %.2f (%.2f\u2013%.2f), n = %d", NM[set], AUC, CI_low, pmin(CI_high, 1), n)]
rc <- merge(rc, au[, .(set, lab)], by = "set"); setorder(rc, set, fpr, tpr); rc[, lab := factor(lab, levels = au$lab)]
COLS1 <- setNames(c("grey45", COL_GRG, "#F4A582"), au$lab[1:3]); COLS2 <- setNames(c("#B2182B", "#F4A582", "#4393C3", "#2166AC"), au$lab[4:7])
roc_p <- function(d, cols, ttl) ggplot(d, aes(fpr, tpr, colour = lab)) + geom_abline(linetype = "dashed", linewidth = 0.3, colour = "grey60") +
  geom_step(linewidth = 0.55, direction = "vh") + coord_equal() + scale_colour_manual(values = cols, name = NULL) +
  scale_x_continuous(breaks = c(0, 0.5, 1)) + scale_y_continuous(breaks = c(0, 0.5, 1)) +
  labs(x = "1 \u2212 specificity", y = "Sensitivity", title = ttl) +
  theme(legend.position = "right", legend.text = element_text(size = 7, lineheight = 0.85), legend.key.height = unit(6.5, "mm"))
pa <- roc_p(rc[set %in% au1$set], COLS1, "Significant fibrosis (F\u22652 vs F0\u20131)")
pb <- roc_p(rc[set %in% au2$set], COLS2, "Disease vs normal liver (model frozen from a)")

co <- fread(file.path(R, "10_gates/03_glm_coef.tsv"))[term != "(Intercept)"]
co[, `:=`(OR = exp(Estimate), lo = exp(Estimate - 1.96 * `Std. Error`), hi = exp(Estimate + 1.96 * `Std. Error`), lab = stars(`Pr(>|z|)`))]
co[, term := factor(term, levels = rev(term))]
pc <- ggplot(co, aes(OR, term)) + geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.3, colour = "grey50") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.2, linewidth = 0.4) + geom_point(shape = 21, size = 2, fill = COL_DN, stroke = 0.3) +
  geom_text(aes(x = 1.15, label = lab), size = GT, family = "Arial", hjust = 0) + scale_x_log10(limits = c(0.18, 1.6), breaks = c(0.25, 0.5, 1)) +
  labs(x = "Odds ratio per SD (F\u22652, training)", y = NULL) + theme(axis.text.y = element_text(face = "italic"))

## d–f 人肝蛋白组（Commun Med 2025，58 例；只报 Kleiner 纤维化模型）
suppressPackageStartupMessages(library(ggtext)); SUB <- function(x) paste0("<sub><span style='font-size:7pt'>", x, "</span></sub>")
P <- file.path(R, "21_protein"); ps <- fread(file.path(P, "04_sample_scores.tsv")); sm <- fread(file.path(P, "01_summary.tsv"))[grepl("Kleiner", model)]
kp <- fread(file.path(P, "02_key_proteins.tsv"))[grepl("Kleiner", model)]
ps[, Fs := factor(paste0("F", fib), levels = paste0("F", 0:3))]; set.seed(SEED); ps[, xj := as.numeric(Fs) + runif(.N, -0.05, 0.05)]
COL_F <- c(F0 = "#4393C3", F1 = "#92C5DE", F2 = "#F4A582", F3 = "#B2182B")
pd <- ggplot(ps, aes(Fs, GRG_score)) + geom_boxplot(aes(fill = Fs), width = 0.5, outlier.shape = NA, linewidth = 0.3, alpha = 0.8) +
  geom_point(aes(x = xj, fill = Fs), shape = 21, size = 1.1, stroke = 0.2) +
  annotate("text", x = 2.5, y = max(ps$GRG_score) * 1.15, label = paste0("camera ", stars(sm$camera_P)), size = GT, family = "Arial") +
  scale_fill_manual(values = COL_F, guide = "none") +
  labs(x = "Kleiner fibrosis stage", y = "Liver GRG protein score (z)", title = sprintf("Human liver proteome (n = %d)", sm$n)) +
  scale_x_discrete(labels = sprintf("%s
(n = %d)", levels(ps$Fs), as.integer(table(ps$Fs)))) + theme(axis.text.x = element_text(lineheight = 0.85))
nz <- data.table(ref = factor(c("GRG", "Expression-
matched", "Mitochondrial
(non-GRG)"), levels = c("Mitochondrial
(non-GRG)", "Expression-
matched", "GRG")),
                 t = c(sm$mean_t, sm$null_median_N1, sm$null_median_N3), lab = c("", stars(sm$null_P_N1), stars(sm$null_P_N3)))
pe <- ggplot(nz, aes(t, ref)) + geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey50") +
  geom_segment(aes(x = 0, xend = t, yend = ref), linewidth = 0.4, colour = "grey60") +
  geom_point(aes(fill = ref), shape = 21, size = 2.2, stroke = 0.3) + geom_text(aes(label = lab), nudge_y = 0.3, size = GT, family = "Arial") +
  scale_fill_manual(values = c(GRG = COL_GRG, `Expression-
matched` = "grey80", `Mitochondrial
(non-GRG)` = "grey45"), guide = "none") +
  scale_x_continuous(breaks = c(-2, -1, 0)) + labs(x = "Mean moderated t
(fibrosis)", y = NULL) + theme(axis.text.y = element_text(lineheight = 0.85))
kp[, gene := factor(gene, levels = rev(kp[order(t), gene]))]
pf <- ggplot(kp, aes(logFC, gene)) + geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey50") +
  geom_segment(aes(x = 0, xend = logFC, yend = gene), linewidth = 0.4, colour = "grey55") +
  geom_point(shape = 21, size = 2, fill = COL_DN, stroke = 0.3) +
  geom_text(aes(x = min(logFC) * 1.3, label = stars(BH)), size = GT, family = "Arial", hjust = 0) +
  scale_x_continuous(limits = c(min(kp$logFC) * 1.35, 0.03)) +
  labs(x = paste0("Protein log", SUB(2), " FC per stage"), y = NULL) + theme(axis.title.x = element_markdown(), axis.text.y = element_text(face = "italic"))

## SHAP 蜂群（训练集 GLM，logit 尺度）
suppressPackageStartupMessages({ library(ggbeeswarm); library(ggdist) })
sh <- fread(file.path(R, "23_upgrade/U6_shap.tsv")); ordg <- sh[, .(m = mean(abs(shap))), by = gene][order(m), gene]
sh[, gene := factor(gene, levels = ordg)]
pshap <- ggplot(sh, aes(shap, gene, colour = pmax(pmin(z, 2.5), -2.5))) + geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3, colour = "grey50") +
  geom_quasirandom(size = 0.6, alpha = 0.85, orientation = "y", width = 0.3) +
  scale_colour_gradient2(low = COL_DN, mid = "grey88", high = COL_UP, name = "Expression
(z)") +
  labs(x = "SHAP value (log-odds of F≥2)", y = NULL, title = "Training set (kernelSHAP)") +
  theme(axis.text.y = element_text(face = "italic"), legend.key.height = unit(4, "mm"), legend.key.width = unit(2, "mm"))
## 蛋白评分：雨云
pd <- ggplot(ps, aes(Fs, GRG_score, fill = Fs)) +
  stat_slab(scale = 0.5, alpha = 0.85, colour = NA, position = position_nudge(x = 0.12), normalize = "xy") +
  geom_boxplot(width = 0.14, outlier.shape = NA, linewidth = 0.3, position = position_nudge(x = -0.02)) +
  geom_point(aes(x = xj - 0.2), shape = 21, size = 1, stroke = 0.2) +
  annotate("text", x = 2.5, y = max(ps$GRG_score) * 1.15, label = paste0("camera ", stars(sm$camera_P)), size = GT, family = "Arial") +
  scale_fill_manual(values = COL_F, guide = "none") +
  labs(x = "Kleiner fibrosis stage", y = "Liver GRG protein score (z)", title = sprintf("Human liver proteome (n = %d)", sm$n)) +
  scale_x_discrete(labels = sprintf("%s
(n = %d)", levels(ps$Fs), as.integer(table(ps$Fs)))) + theme(axis.text.x = element_text(lineheight = 0.85))

rows <- list(list(height = 60, cells = list(list(width = 104, obj = pa, tag = "a"), list(width = 70, obj = pshap, tag = "b"))),
             list(height = 62, cells = list(list(width = 118, obj = pb, tag = "c"), list(width = 56, obj = pc, tag = "d"))),
             list(height = 56, cells = list(list(width = 62, obj = pd, tag = "e"), list(width = 56, obj = pe, tag = "f"), list(width = 56, obj = pf, tag = "g"))))
save_rows(rows, "Fig4")
