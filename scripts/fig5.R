# Fig. 5 FXR：靶基因富集（相关）与激动剂干预（只读盘画图）
source(file.path(Sys.getenv("GLUT_PROJ", "."), "scripts/fig_config.R"))
suppressPackageStartupMessages({ library(ggdist); library(ggtext) })
SUB <- function(x) paste0("<sub><span style='font-size:7pt'>", x, "</span></sub>")

## a NR1H4 靶基因 GSEA（肝细胞，严重度排序）
cv <- fread(file.path(R, "20_figsrc/F5_NR1H4_curve.tsv")); tk <- fread(file.path(R, "20_figsrc/F5_NR1H4_ticks.tsv"))
fg <- fread(file.path(R, "14_batch1/A5_TF_GTRD_fgsea.tsv"))[TF == "NR1H4"]
rs <- fread(file.path(R, "20_figsrc/F5_rank_stats.tsv"))
a1 <- ggplot(cv, aes(rank, ES)) + geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey50") +
  geom_area(fill = COL_DN, alpha = 0.18) + geom_line(colour = COL_DN, linewidth = 0.7) +
  annotate("text", x = max(cv$rank) * 0.05, y = min(cv$ES) * 0.9, label = sprintf("NES = %.2f
FDR %s", fg$NES, stars(fg$padj)), hjust = 0, vjust = 0,
           size = GT, family = "Arial", lineheight = 0.9) +
  scale_x_continuous(expand = c(0, 0)) + labs(x = NULL, y = "Enrichment
score", title = "NR1H4 (FXR) target genes (GTRD) in hepatocytes") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank())
a2 <- ggplot(tk) + geom_segment(aes(x = rank, xend = rank, y = 0, yend = 1, colour = in_GRG), linewidth = 0.25) +
  scale_colour_manual(values = c(`TRUE` = COL_GRG, `FALSE` = "grey35"), labels = c(`TRUE` = "GRG member", `FALSE` = "Other target"), name = NULL) +
  scale_x_continuous(limits = range(rs$rank), expand = c(0, 0)) + theme_void(base_family = "Arial") +
  theme(legend.position = "right", legend.text = element_text(size = 7), legend.key.width = unit(3, "mm"))
a3 <- ggplot(rs[seq(1, .N, by = 5)], aes(rank, stat)) + geom_area(aes(fill = stat > 0)) + geom_hline(yintercept = 0, linewidth = 0.2) +
  scale_fill_manual(values = c(`TRUE` = COL_UP, `FALSE` = COL_DN), guide = "none") + scale_x_continuous(expand = c(0, 0)) +
  labs(x = "Genes ranked by severity t (up → down)", y = "t")
pa <- (a1 + labs(tag = "a")) / a2 / a3 + plot_layout(heights = c(3, 0.45, 1.1), guides = "collect") & theme(legend.position = "bottom")

## b 肝细胞 NR1H4 表达
dn <- fread(file.path(R, "20_figsrc/F2_hep_donor_scores.tsv")); dn[, group := factor(group, levels = names(COL_SEV))]
set.seed(SEED); dn[, xj := as.numeric(group) - 0.2 + runif(.N, -0.06, 0.06)]
pb <- ggplot(dn, aes(group, NR1H4, fill = group)) +
  stat_slab(scale = 0.55, alpha = 0.85, colour = NA, position = position_nudge(x = 0.12), normalize = "xy") +
  geom_boxplot(width = 0.14, outlier.shape = NA, linewidth = 0.3, position = position_nudge(x = -0.02)) +
  geom_point(aes(x = xj), shape = 21, size = 1, stroke = 0.2, alpha = 0.85) +
  annotate("text", x = 2.5, y = max(dn$NR1H4) + 0.15, label = paste0("severity trend ", stars(fg$TF_P)), size = GT, family = "Arial") +
  scale_fill_manual(values = COL_SEV, guide = "none") + labs(x = NULL, y = paste0("<i>NR1H4</i> in hepatocytes (log", SUB(2), " CPM)")) +
  theme(axis.title.y = element_markdown())

## c 激动剂处理后 GRG 基因 logFC 分布 + camera
sm <- fread(file.path(R, "18_fxr/18_fxr_summary.tsv")); gg <- fread(file.path(R, "18_fxr/18_fxr_GRG_genes.tsv"))
DS <- c(GSE76161_PCLS_OCA_vs_DMSO = "Human liver slices\nOCA 24 h (n = 3 pairs)", `GSE124053_HepaRG_FXR-L_4h` = "HepaRG FXR ligand 4 h",
        `GSE124053_HepaRG_FXR-L_24h` = "HepaRG FXR ligand 24 h", GSE124053_HepaRG_CDCA_4h = "HepaRG CDCA 4 h", GSE124053_HepaRG_CDCA_24h = "HepaRG CDCA 24 h")
gg <- gg[dataset %in% names(DS)]; gg[, ds := factor(DS[dataset], levels = rev(DS))]
st <- sm[dataset %in% names(DS), .(ds = factor(DS[dataset], levels = rev(DS)), lab = paste0(fifelse(direction == "Down", "\u2193 ", "\u2191 "), stars(camera_P)))]
pc <- ggplot(gg, aes(logFC, ds)) + geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey50") +
  stat_slab(fill = COL_DN, alpha = 0.55, scale = 0.7, colour = NA, normalize = "groups") +
  geom_boxplot(width = 0.18, outlier.shape = NA, linewidth = 0.3, fill = "white") +
  geom_text(data = st, aes(x = 0.95, label = lab), size = GT, family = "Arial", hjust = 1, vjust = -0.6) +
  coord_cartesian(xlim = c(-1, 1)) + labs(x = paste0("GRG gene log", SUB(2), " FC (agonist vs vehicle)"), y = NULL) +
  theme(axis.title.x = element_markdown(), axis.text.y = element_text(lineheight = 0.85))

## d 阳性对照与关键 GRG 基因
pcg <- fread(file.path(R, "18_fxr/18_fxr_posctrl.tsv"))[dataset %in% names(DS)]
key <- c("GCDH", "AASS", "ACADSB", "SARDH", "ACOT13", "ACAA2")
hm <- rbind(pcg[, .(dataset, gene, logFC, P, grp = "FXR targets\n(pos. controls)")], gg[gene %in% key, .(dataset, gene, logFC, P, grp = "GRG key genes")])
hm[, ds := factor(DS[dataset], levels = rev(DS))]; hm[, gene := factor(gene, levels = c("NR0B2", "ABCB11", "SLC51B", key))]
hm[, lab := fifelse(is.na(logFC), "NA", stars(P))]
pd <- ggplot(hm, aes(gene, ds, fill = logFC)) + geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(data = hm[is.na(logFC)], aes(label = lab), size = GT, family = "Arial") + facet_grid(~ grp, scales = "free_x", space = "free_x") +
  scale_fill_gradient2(low = COL_DN, mid = "white", high = COL_UP, limits = c(-2, 2), oob = scales::squish, na.value = "grey90",
                       name = paste0("log2 FC")) +
  labs(x = NULL, y = NULL) + theme(axis.text.x = element_text(angle = 40, hjust = 1, face = "italic"), axis.text.y = element_text(lineheight = 0.85),
  axis.ticks.y = element_blank(), axis.line = element_blank(), legend.position = "right")

rows <- list(list(height = 64, cells = list(list(width = 96, obj = pa), list(width = 78, obj = pb, tag = "b"))),
             list(height = 64, cells = list(list(width = 78, obj = pc, tag = "c"), list(width = 96, obj = pd + theme(axis.text.y = element_text(size = 7, lineheight = 0.85)), tag = "d"))))
save_rows(rows, "Fig5")
