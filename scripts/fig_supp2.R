# 补充图 S6–S11（只读盘画图；图上只星号）
source(file.path(Sys.getenv("GLUT_PROJ", "."), "scripts/fig_config.R"))
suppressPackageStartupMessages({ library(ggtext); library(ggrepel); library(ggdist); library(rms) })
SUB <- function(x) paste0("<sub><span style='font-size:7pt'>", x, "</span></sub>")
U <- file.path(R, "23_upgrade"); V <- file.path(R, "24_supp2"); B1 <- file.path(R, "14_batch1")
hub <- readLines(file.path(R, "10_gates/02_hub.txt")); grg <- readLines(file.path(PROJ, "geneset/GRG_full.txt"))

## ---------- S6 LASSO 细节 ----------
lc <- fread(file.path(U, "U4_lasso_cv.tsv"))
s6a <- ggplot(lc, aes(log(lambda), cvm)) + geom_errorbar(aes(ymin = cvm - cvsd, ymax = cvm + cvsd), colour = "grey70", width = 0, linewidth = 0.3) +
  geom_point(colour = COL_GRG, size = 0.9) + geom_vline(xintercept = log(lc$lambda_min[1]), linetype = "dashed", linewidth = 0.3) +
  geom_vline(xintercept = log(lc$lambda_1se[1]), linetype = "dotted", linewidth = 0.4) +
  annotate("text", x = log(lc$lambda_1se[1]), y = max(lc$cvm), label = "\u03bb.1se (used)  ", hjust = 1, size = GT, family = "Arial") +
  labs(x = "log(\u03bb)", y = "Binomial deviance (10-fold CV)")
lp <- fread(file.path(U, "U4_lasso_path.tsv")); lp[, hub := gene %in% hub]
lab <- lp[lambda == min(lambda) & hub == TRUE]
s6b <- ggplot(lp, aes(log(lambda), coef, group = gene)) + geom_line(data = lp[hub == FALSE], colour = "grey80", linewidth = 0.3) +
  geom_line(data = lp[hub == TRUE], aes(colour = gene), linewidth = 0.6) + geom_hline(yintercept = 0, linewidth = 0.25) +
  geom_vline(xintercept = log(lc$lambda_1se[1]), linetype = "dotted", linewidth = 0.4) +
  geom_text_repel(data = lab, aes(label = gene, colour = gene), size = GT, family = "Arial", fontface = "italic", direction = "y", nudge_x = -0.3, segment.size = 0.2) +
  scale_colour_manual(values = c("#B2182B", "#D6604D", "#2166AC", "#4393C3"), guide = "none") + labs(x = "log(\u03bb)", y = "Coefficient")
save_rows(list(list(height = 66, cells = list(list(width = 84, obj = s6a, tag = "a"), list(width = 90, obj = s6b, tag = "b")))), "FigS6")

## ---------- S7 列线图 + bootstrap 校准 ----------
ob <- readRDS(file.path(V, "S7_nomogram_objects.rds")); dd <- ob$dd; options(datadist = "dd")
nom <- nomogram(ob$fit, ACADSB = seq(-4, 2, 1), ACOT13 = seq(-3, 4, 1), ACAA2 = seq(-3, 4, 1), SARDH = seq(-2, 2, 1), fun = plogis, funlabel = "Risk of F\u22652", fun.at = c(0.1, 0.3, 0.5, 0.7, 0.9), lp = FALSE)
draw_nom <- function() { par(mar = c(1, 1, 1, 1), family = "Arial", ps = 7); plot(nom, xfrac = 0.3, cex.axis = 0.6, cex.var = 0.7, lmgp = 0.15, tcl = -0.15) }
cb <- fread(file.path(V, "S7_calibration_boot.tsv")); vb <- fread(file.path(V, "S7_validate_boot.tsv"))
cl <- melt(cb, id.vars = "pred", variable.name = "type", value.name = "obs")
cl[, type := factor(type, levels = c("apparent", "bias_corrected"), labels = c("Apparent", "Bias-corrected (200 bootstraps)"))]
cidx <- vb[index == "Dxy", .(o = index.orig / 2 + 0.5, c = index.corrected / 2 + 0.5)]
s7b <- ggplot(cl, aes(pred, obs, colour = type, linetype = type)) + geom_abline(linetype = "dashed", colour = "grey60", linewidth = 0.3) + geom_line(linewidth = 0.6) +
  scale_colour_manual(values = c(Apparent = "grey40", `Bias-corrected (200 bootstraps)` = COL_GRG), name = NULL) +
  scale_linetype_manual(values = c(Apparent = "dashed", `Bias-corrected (200 bootstraps)` = "solid"), name = NULL) +
  annotate("text", x = 0.05, y = 0.95, label = sprintf("C-index %.2f (optimism-corrected %.2f)", cidx$o, cidx$c), hjust = 0, size = GT, family = "Arial") +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) + labs(x = "Predicted probability", y = "Observed proportion", title = "GSE135251 (training)") +
  theme(legend.position = "bottom", legend.direction = "vertical")
save_rows(list(list(height = 82, cells = list(list(width = 104, obj = draw_nom, tag = "a"), list(width = 70, obj = s7b, tag = "b")))), "FigS7")

## ---------- S8 GRG 功能组成与细胞类型分布 ----------
eg <- fread(file.path(U, "U3_GRG_GO_BP.tsv"))[order(p.adjust)][1:15]
eg[, `:=`(ratio = Count / as.numeric(sub(".*/", "", GeneRatio)), term = factor(stringr::str_wrap(Description, 38), levels = rev(stringr::str_wrap(Description, 38))))]
s8a <- ggplot(eg, aes(ratio, term)) + geom_segment(aes(x = 0, xend = ratio, yend = term), colour = "grey70", linewidth = 0.4) +
  geom_point(aes(size = Count, fill = -log10(p.adjust)), shape = 21, stroke = 0.3) +
  scale_fill_gradient(low = "#FDDBC7", high = "#B2182B", name = "\u2212log10\nFDR") + scale_size_continuous(range = c(1.5, 4), name = "Genes") +
  labs(x = "Gene ratio", y = NULL, title = "GO biological process (GRG vs expressed genes)") + theme(axis.text.y = element_text(lineheight = 0.85))
cm <- fread(file.path(V, "S8_GRG_celltype_logCPM.tsv")); M <- as.matrix(cm[, -1]); rownames(M) <- cm$gene
Z <- t(scale(t(M))); Z <- Z[complete.cases(Z), ]; ho <- hclust(dist(Z))$order
zl <- as.data.table(as.table(Z)); setnames(zl, c("gene", "celltype", "z")); zl[, gene := factor(gene, levels = rownames(Z)[ho])]
zl[, celltype := factor(celltype, levels = names(COL_CELLTYPE))]
s8b <- ggplot(zl, aes(celltype, gene, fill = pmax(pmin(z, 2.5), -2.5))) + geom_tile() +
  scale_fill_gradient2(low = COL_DN, mid = "white", high = COL_UP, name = "z") +
  labs(x = NULL, y = sprintf("GRG genes (n = %d, clustered)", nrow(Z))) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line = element_blank(), axis.text.x = element_text(angle = 40, hjust = 1))
save_rows(list(list(height = 92, cells = list(list(width = 112, obj = s8a, tag = "a"), list(width = 62, obj = s8b, tag = "b")))), "FigS8")

## ---------- S9 肝细胞亚状态供者级比例 ----------
pp <- fread(file.path(B1, "A3_hepstate_proportions.tsv")); pr <- fread(file.path(B1, "A3_hepstate_proportion_vs_severity.tsv"))
pp[, group := factor(names(COL_SEV)[severity + 1], levels = names(COL_SEV))]
pr[, lab := sprintf("\u03c1 = %.2f %s", rho, stars(BH))]
s9 <- ggplot(pp, aes(group, frac, fill = group)) + geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = 0.25) +
  geom_point(position = position_jitter(width = 0.12, seed = SEED), shape = 21, size = 0.8, stroke = 0.15) +
  geom_text(data = pr, aes(x = 2.5, y = Inf, label = lab), vjust = 1.4, size = GT, family = "Arial", inherit.aes = FALSE) +
  scale_fill_manual(values = COL_SEV, guide = "none") + scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0.02, 0.2))) +
  facet_wrap(~ state, nrow = 2, scales = "free_y") + labs(x = NULL, y = "Fraction of donor hepatocytes") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), strip.text = element_text(face = "bold"))
save_rows(list(list(height = 100, cells = list(list(width = 174, obj = s9)))), "FigS9")

## ---------- S10 上游 TF 全景 + FXR 配体下 GRG 基因 ----------
tf <- fread(file.path(B1, "A5_TF_GTRD_fgsea.tsv"))[padj < 0.05][order(NES)][, TF := sub("_FUSION.*", " fusion", TF)]
tf[, TF := factor(TF, levels = TF)]; tf[, hl := TF == "NR1H4"]
s10a <- ggplot(tf, aes(NES, TF)) + geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey50") +
  geom_segment(aes(x = 0, xend = NES, yend = TF), colour = "grey70", linewidth = 0.35) +
  geom_point(aes(fill = hl, size = size), shape = 21, stroke = 0.3) +
  scale_fill_manual(values = c(`TRUE` = COL_GRG, `FALSE` = "grey60"), guide = "none") + scale_size_continuous(range = c(1, 3), name = "Targets") +
  labs(x = "NES (hepatocyte severity ranking)", y = NULL, title = "GTRD TF targets (FDR < 0.05)") + theme(plot.title.position = "plot") + theme(axis.text.y = element_text(size = 7))
fx <- fread(file.path(R, "18_fxr/18_fxr_GRG_genes.tsv"))[dataset == "GSE124053_HepaRG_FXR-L_24h"]
pc <- fread(file.path(R, "18_fxr/18_fxr_posctrl.tsv"))[dataset == "GSE124053_HepaRG_FXR-L_24h"]
vv <- rbind(fx[, .(gene, logFC, P, cls = "GRG")], pc[, .(gene, logFC, P, cls = "FXR target")])
s10b <- ggplot(vv, aes(logFC, -log10(P))) + geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey50") +
  geom_point(aes(fill = cls, size = cls), shape = 21, stroke = 0.2, alpha = 0.85) +
  geom_text_repel(data = vv[cls == "FXR target" | gene %in% c(hub, "GCDH", "AASS")], aes(label = gene), size = GT, family = "Arial", fontface = "italic",
                  segment.size = 0.2, max.overlaps = Inf) +
  scale_fill_manual(values = c(GRG = COL_DN, `FXR target` = COL_UP), name = NULL) + scale_size_manual(values = c(GRG = 1.2, `FXR target` = 2.2), guide = "none") +
  labs(x = paste0("log", SUB(2), " FC (FXR ligand 24 h vs vehicle)"), y = paste0("\u2212log", SUB(10), " <i>P</i>"), title = "HepaRG (GSE124053)") +
  theme(axis.title.x = element_markdown(), axis.title.y = element_markdown(), legend.position = c(0.02, 0.98), legend.justification = c(0, 1))
save_rows(list(list(height = 100, cells = list(list(width = 82, obj = s10a, tag = "a"), list(width = 92, obj = s10b, tag = "b")))), "FigS10")

## ---------- S11 人肝蛋白组：全蛋白火山 + 关键蛋白按分期 ----------
pa <- fread(file.path(V, "S11_protein_all_kleiner.tsv")); setorder(pa, GRG)
pa[, cls := fifelse(GRG, "GRG protein", "Other protein")]
s11a <- ggplot(pa, aes(logFC, -log10(P), colour = cls)) + geom_point(size = 0.5, alpha = 0.6, stroke = 0) +
  geom_point(data = pa[GRG == TRUE], size = 1.1, alpha = 0.9, stroke = 0) +
  geom_text_repel(data = pa[gene %in% c(hub, "GCDH", "AASS")], aes(label = gene), colour = "black", size = GT, family = "Arial", fontface = "italic",
                  segment.size = 0.2, max.overlaps = Inf, box.padding = 0.4) +
  scale_colour_manual(values = c(`GRG protein` = COL_DN, `Other protein` = "grey80"), name = NULL) +
  labs(x = paste0("log", SUB(2), " FC per fibrosis stage"), y = paste0("\u2212log", SUB(10), " <i>P</i>"), title = sprintf("Human liver proteome (%s proteins)", format(nrow(pa), big.mark = ","))) +
  theme(axis.title.x = element_markdown(), axis.title.y = element_markdown(), legend.position = c(0.98, 0.98), legend.justification = c(1, 1))
ps <- fread(file.path(R, "21_protein/04_sample_scores.tsv")); kp <- fread(file.path(R, "21_protein/02_key_proteins.tsv"))[grepl("Kleiner", model)]
pl <- melt(ps, id.vars = c("patient_name", "fib"), measure.vars = c("GCDH", "AASS", "ACADSB", "SARDH", "ACOT13", "ACAA2"), variable.name = "gene")
pl[, z := (value - mean(value, na.rm = TRUE)) / sd(value, na.rm = TRUE), by = gene]; pl[, Fs := factor(paste0("F", fib))]
kl <- kp[, .(gene = factor(gene, levels = levels(pl$gene)), lab = stars(BH))]
s11b <- ggplot(pl, aes(Fs, z, fill = Fs)) + geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = 0.25) +
  geom_point(position = position_jitter(width = 0.12, seed = SEED), size = 0.5, alpha = 0.6, stroke = 0) +
  geom_text(data = kl, aes(x = 2.5, y = Inf, label = lab), vjust = 1.3, size = GT, family = "Arial", inherit.aes = FALSE) +
  scale_fill_manual(values = c(F0 = "#4393C3", F1 = "#92C5DE", F2 = "#F4A582", F3 = "#B2182B"), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) + facet_wrap(~ gene, nrow = 2) +
  labs(x = "Kleiner fibrosis stage", y = "Protein abundance (z)") + theme(strip.text = element_text(face = "bold.italic"))
save_rows(list(list(height = 82, cells = list(list(width = 80, obj = s11a, tag = "a"), list(width = 94, obj = s11b, tag = "b")))), "FigS11")
