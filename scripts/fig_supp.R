# 补充图 S1–S5（只读盘画图）。统计数字全部取自已落盘的 TSV，图上只星号
source(file.path(Sys.getenv("GLUT_PROJ", "."), "scripts/fig_config.R"))
suppressPackageStartupMessages({ library(ggtext); library(pROC) })
GS <- file.path(PROJ, "geneset"); B1 <- file.path(R, "14_batch1")

## ---------- Fig. S1 基因集构建细节 ----------
mm <- fread(file.path(GS, "GRG_substrate_mouse2human.tsv"))
ov <- data.table(cat = factor(c("Wagner 2017 only", "Schmiesing 2018 only", "Both"), levels = c("Wagner 2017 only", "Schmiesing 2018 only", "Both")),
                 n = c(sum(mm$source == "Wagner2017"), sum(mm$source == "Schmiesing2018"), sum(mm$both == 1)))
s1a <- ggplot(ov, aes(cat, n, fill = cat)) + geom_col(width = 0.6) + geom_text(aes(label = n), vjust = -0.4, size = GT, family = "Arial") +
  scale_fill_manual(values = c("#92C5DE", "#4393C3", COL_GRG), guide = "none") + scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "Mouse genes", title = "Overlap of liver glutarylomes") + theme(axis.text.x = element_text(angle = 25, hjust = 1))
oc <- mm[, .N, by = .(ortholog, both)][, both := factor(both, labels = c("One dataset", "Both datasets"))]
oc[, ortholog := factor(ortholog, levels = c("1:1", "multi", "none"), labels = c("1:1 orthologue", "Multiple orthologues", "No orthologue"))]
s1b <- ggplot(oc, aes(ortholog, N, fill = both)) + geom_col(width = 0.6) +
  scale_fill_manual(values = c(`One dataset` = "#92C5DE", `Both datasets` = COL_GRG), name = NULL) + scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = NULL, y = "Mouse genes", title = "MGI mouse\u2013human mapping") + theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = c(0.98, 0.98), legend.justification = c(1, 1))
pg <- rbind(fread(file.path(GS, "wagner2017_genes_q01.tsv")), fread(file.path(GS, "schmiesing2018_genes_q01.tsv")))[organism == "Mus musculus"]
pg[, dataset := factor(dataset, labels = c("Schmiesing 2018", "Wagner 2017"))]
s1c <- ggplot(pg, aes(n_unique_peptides, fill = dataset)) + geom_histogram(binwidth = 1, position = "identity", alpha = 0.65, colour = NA) +
  scale_fill_manual(values = c(`Wagner 2017` = "#92C5DE", `Schmiesing 2018` = "#4393C3"), name = NULL) + scale_x_continuous(limits = c(0, 30)) +
  labs(x = "Kglu peptides per protein", y = "Proteins") + theme(legend.position = c(0.98, 0.98), legend.justification = c(1, 1))
save_rows(list(list(height = 62, cells = list(list(width = 52, obj = s1a, tag = "a"), list(width = 62, obj = s1b, tag = "b"), list(width = 60, obj = s1c, tag = "c")))), "FigS1")

## ---------- Fig. S2 4 个 hub 在三套 bulk 中按分期 ----------
s2 <- fread(file.path(R, "22_supp/S2_hub_by_stage.tsv"))
s2[, cohort := factor(cohort, levels = c("GSE135251", "GSE162694", "GSE130970"), labels = c("GSE135251 (training)", "GSE162694 (external)", "GSE130970 (external)"))]
s2[, gene := factor(gene, levels = readLines(file.path(R, "10_gates/02_hub.txt")))]
s2[, Fs := factor(paste0("F", fibrosis))]
p2 <- ggplot(s2, aes(Fs, z, fill = Fs)) + geom_boxplot(width = 0.6, outlier.shape = NA, linewidth = 0.25) +
  geom_jitter(width = 0.15, size = 0.3, alpha = 0.4, stroke = 0) +
  scale_fill_manual(values = setNames(c("#4393C3", "#92C5DE", "#F7F7F7", "#F4A582", "#B2182B"), paste0("F", 0:4)), guide = "none") +
  facet_grid(gene ~ cohort) + labs(x = "Fibrosis stage", y = "Expression (z within cohort)") +
  theme(strip.text.y = element_text(face = "italic"), strip.background = element_rect(fill = "grey95", colour = NA))
save_rows(list(list(height = 150, cells = list(list(width = 174, obj = p2)))), "FigS2")

## ---------- Fig. S3 全部 5 张 Visium 切片 ----------
sp <- fread(file.path(R, "20_figsrc/F3_spatial_spots.tsv")); rng <- quantile(sp$GRG_full, c(0.01, 0.99))
lab <- c(JBO14 = "JBO14 (H35, steatotic)", JBO15 = "JBO15 (H35, steatotic)", JBO19 = "JBO19 (H37, steatotic)", JBO18 = "JBO18 (H36, healthy)", JBO22 = "JBO22 (H38, healthy)")
mk <- function(var, pal, lims, nm) lapply(names(lab), function(s) ggplot(sp[sample == s], aes(pxl_col, -pxl_row, colour = get(var))) +
  geom_point(size = 0.35, stroke = 0) + coord_equal() + scale_colour_gradientn(colours = pal, limits = lims, oob = scales::squish, name = nm) +
  labs(title = lab[[s]]) + theme_void(base_family = "Arial") +
  theme(plot.title = element_text(size = 7, hjust = 0.5), legend.text = element_text(size = 7), legend.title = element_text(size = 7),
        legend.key.width = unit(2, "mm"), legend.key.height = unit(4, "mm"), plot.tag = element_text(size = 10, face = "bold", family = "Arial")))
l1 <- mk("GRG_full", c("#F7FBFF", "#FDDBC7", "#F4A582", "#D6604D", "#67001F"), rng, "GRG"); l1[[1]] <- l1[[1]] + labs(tag = "a")
l2 <- mk("zonation", c("#2166AC", "#F7F7F7", "#B2182B"), c(0, 1), "Zonation"); l2[[1]] <- l2[[1]] + labs(tag = "b")
g1 <- wrap_plots(l1, nrow = 1, guides = "collect"); g2 <- wrap_plots(l2, nrow = 1, guides = "collect")
save_rows(list(list(height = 48, cells = list(list(width = 174, obj = g1))),
               list(height = 48, cells = list(list(width = 174, obj = g2)))), "FigS3")

## ---------- Fig. S4 校准与 DCA ----------
pr <- rbind(fread(file.path(R, "20_figsrc/F4_F2plus_predictions.tsv"))[set == "GSE162694_external", .(set = "GSE162694 (F\u22652)", y, p = score)],
            fread(file.path(R, "13_hep_mash/13_EMASH_predictions.tsv"))[comparison == "GSE126848_NASH_vs_healthy", .(set = "GSE126848 (MASH vs normal)", y, p = plogis(score))])
pr[, bin := cut(p, unique(quantile(p, 0:5 / 5)), include.lowest = TRUE, labels = FALSE), by = set]
cb <- pr[, .(pred = mean(p), obs = mean(y), n = .N), by = .(set, bin)]
cb[, `:=`(lo = pmax(0, obs - 1.96 * sqrt(obs * (1 - obs) / n)), hi = pmin(1, obs + 1.96 * sqrt(obs * (1 - obs) / n)))]
ca <- fread(file.path(B1, "C3_calibration.tsv"))[, set := c(GSE162694_F2plus = "GSE162694 (F\u22652)", GSE126848_MASH_vs_normal = "GSE126848 (MASH vs normal)")[set]]
s4a <- ggplot(cb, aes(pred, obs, colour = set)) + geom_abline(linetype = "dashed", linewidth = 0.3, colour = "grey60") +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.02, linewidth = 0.3) + geom_line(linewidth = 0.4) + geom_point(size = 1.4) +
  scale_colour_manual(values = c("GSE162694 (F\u22652)" = COL_GRG, "GSE126848 (MASH vs normal)" = "#4393C3"), name = NULL) + coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "Predicted probability (quintile mean)", y = "Observed proportion") + theme(legend.position = "bottom", legend.direction = "vertical")
dc <- fread(file.path(B1, "C3_DCA.tsv"))[, set := c(GSE162694_F2plus = "GSE162694 (F\u22652)", GSE126848_MASH_vs_normal = "GSE126848 (MASH vs normal)")[set]]
dl <- melt(dc, id.vars = c("set", "threshold"), variable.name = "strategy")[, strategy := factor(strategy, levels = c("NB_model", "NB_all", "NB_none"), labels = c("4-gene model", "Treat all", "Treat none"))]
s4b <- ggplot(dl, aes(threshold, value, colour = strategy, linetype = strategy)) + geom_line(linewidth = 0.5) + facet_wrap(~ set, nrow = 1) +
  scale_colour_manual(values = c(`4-gene model` = COL_GRG, `Treat all` = "grey40", `Treat none` = "grey70"), name = NULL) +
  scale_linetype_manual(values = c(`4-gene model` = "solid", `Treat all` = "dashed", `Treat none` = "dotted"), name = NULL) +
  coord_cartesian(ylim = c(-0.1, max(dl$value) * 1.05)) + labs(x = "Threshold probability", y = "Net benefit") + theme(legend.position = "bottom")
save_rows(list(list(height = 78, cells = list(list(width = 66, obj = s4a, tag = "a"), list(width = 108, obj = s4b, tag = "b")))), "FigS4")

## ---------- Fig. S5 bulk 病理分项与疾病谱 ----------
s5 <- fread(file.path(R, "22_supp/S5_pathology_long.tsv")); c1 <- fread(file.path(B1, "C1_pathology_components.tsv"))
s5[, variable := factor(variable, levels = c("steatosis", "inflammation", "ballooning", "fibrosis", "NAS"), labels = c("Steatosis", "Inflammation", "Ballooning", "Fibrosis", "NAS"))]
c1[, variable := factor(variable, levels = c("steatosis", "inflammation", "ballooning", "fibrosis", "NAS"), labels = c("Steatosis", "Inflammation", "Ballooning", "Fibrosis", "NAS"))]
c1[, lab := sprintf("\u03c1 = %.2f %s", rho, stars(BH))]
s5a <- ggplot(s5, aes(factor(value), score)) + geom_boxplot(width = 0.6, outlier.shape = NA, linewidth = 0.25, fill = "grey92") +
  geom_jitter(width = 0.15, size = 0.35, alpha = 0.5, stroke = 0, colour = COL_GRG) +
  geom_text(data = c1, aes(x = -Inf, y = Inf, label = lab), hjust = -0.08, vjust = 1.3, size = GT, family = "Arial", inherit.aes = FALSE) +
  facet_grid(cohort ~ variable, scales = "free_x", space = "free_x") + labs(x = "Histological grade / stage", y = "GRG score (z)") +
  theme(strip.background = element_rect(fill = "grey95", colour = NA))
c2 <- fread(file.path(B1, "C2_GSE164760_scores.tsv"))[!is.na(ord)]
c2[, tissue := factor(tissue, levels = c("Healthy liver", "NASH liver", "Non-tumoral NASH liver adjacent to HCC", "NASH-HCC tumor"),
                      labels = c("Healthy", "MASH", "MASH-HCC\nnon-tumour", "MASH-HCC\ntumour"))]
cs <- fread(file.path(B1, "C2_GSE164760_summary.tsv"))
s5b <- ggplot(c2, aes(tissue, score, fill = tissue)) + geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = 0.3) +
  geom_jitter(width = 0.15, size = 0.6, alpha = 0.6, stroke = 0) +
  annotate("text", x = 2.5, y = max(c2$score) * 1.1, label = paste0("Kruskal\u2013Wallis ", stars(cs[stat == "KW_P", value])), size = GT, family = "Arial") +
  scale_fill_manual(values = c("#4393C3", "#F4A582", "#D6604D", "#67001F"), guide = "none") +
  scale_x_discrete(labels = function(x) sprintf("%s\n(n = %d)", x, as.integer(table(c2$tissue)[x]))) +
  labs(x = NULL, y = "GRG score (z)", title = "GSE164760") + theme(axis.text.x = element_text(lineheight = 0.85))
save_rows(list(list(height = 92, cells = list(list(width = 174, obj = s5a, tag = "a"))),
               list(height = 64, cells = list(list(width = 90, obj = s5b, tag = "b")))), "FigS5")
