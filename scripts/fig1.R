# Fig. 1（升级版）GRG 构建、功能组成与 F≥2 筛选（只读盘画图）；v1 备份 _fig1_v1_backup.R
source(file.path(Sys.getenv("GLUT_PROJ", "."), "scripts/fig_config.R"))
suppressPackageStartupMessages({ library(ggrepel); library(ggtext); library(circlize) })
SUB <- function(x) paste0("<sub><span style='font-size:7pt'>", x, "</span></sub>")
GS <- file.path(PROJ, "geneset"); U <- file.path(R, "23_upgrade"); G10 <- file.path(R, "10_gates")
grg <- readLines(file.path(GS, "GRG_full.txt")); hub <- readLines(file.path(G10, "02_hub.txt"))

## a 基因-GO 弦图：前 30 个 BP 中按 Jaccard<0.4 贪婪去冗余取 6 个；基因限 F≥2 DEG（FDR<0.05）
eg <- fread(file.path(U, "U3_GRG_GO_BP.tsv"))[order(p.adjust)][1:30]; gl <- strsplit(eg$geneID, "/")
pick <- integer(0); for (i in seq_len(nrow(eg))) { if (length(pick) == 6) break
  if (all(sapply(pick, function(j) length(intersect(gl[[i]], gl[[j]])) / length(union(gl[[i]], gl[[j]])) < 0.4))) pick <- c(pick, i) }
fc <- fread(file.path(U, "U3_GRG_logFC.tsv"))
lk <- rbindlist(lapply(pick, function(i) data.table(term = eg$Description[i], gene = gl[[i]])))[gene %in% fc[adj.P.Val < 0.05, gene]]
gcount <- lk[, .N, by = gene][order(-N)]; keep_g <- head(gcount$gene, 26); lk <- lk[gene %in% keep_g]
fwrite(lk, file.path(FIGP, "Fig1a_chord_links.tsv"), sep = "\t")
terms <- unique(lk$term); genes <- unique(lk$gene)
col_fc <- colorRamp2(c(-0.5, 0, 0.5), c(COL_DN, "white", COL_UP))
col_trm <- setNames(c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F")[seq_along(terms)], terms)
grid_col <- c(col_trm, setNames(col_fc(fc$logFC[match(genes, fc$gene)]), genes))
short <- function(x) { x <- sub("^generation of precursor metabolites and energy$", "precursor metabolites
and energy", x)
  x <- sub("^tricarboxylic acid cycle$", "TCA cycle", x); x <- sub(" metabolic process$", "
metabolism", x); x <- sub(" catabolic process$", "
catabolism", x)
  sub("^proteinogenic ", "", x) }
wrap_trm <- setNames(short(terms), terms); ord <- c(terms, rev(genes))  # 图上缩写；完整 GO 名称见 Fig1a_chord_links.tsv 与图注
gaps <- setNames(rep(1.5, length(ord)), ord); gaps[c(terms[length(terms)], genes[1])] <- 8
draw_chord <- function() {
  par(mar = c(0, 0, 0, 0), family = "Arial"); circos.clear()
  circos.par(gap.after = gaps[ord], start.degree = 175, cell.padding = c(0, 0, 0, 0), canvas.xlim = c(-1.5, 1.5), canvas.ylim = c(-1.5, 1.5))
  chordDiagram(lk[, .(term, gene, v = 1)], order = ord, grid.col = grid_col, col = scales::alpha(col_trm[lk$term], 0.65),
               annotationTrack = "grid", annotationTrackHeight = 0.05, preAllocateTracks = list(track.height = 0.18))
  circos.track(track.index = 1, panel.fun = function(x, y) { s <- get.cell.meta.data("sector.index"); ig <- s %in% genes
    circos.text(CELL_META$xcenter, CELL_META$ylim[1] + mm_y(0.8), if (ig) s else wrap_trm[s], facing = "clockwise", niceFacing = TRUE,
                adj = c(0, 0.5), cex = 0.6, font = if (ig) 3 else 1, family = "Arial") }, bg.border = NA)
  circos.clear()
}

## b 火山图
dg <- fread(file.path(G10, "01_DEG_all.tsv"))
dg[, cls := fifelse(adj.P.Val < 0.05 & gene %in% grg, fifelse(logFC > 0, "GRG up", "GRG down"), fifelse(adj.P.Val < 0.05, "Other DEG", "NS"))]
dg[, cls := factor(cls, levels = c("NS", "Other DEG", "GRG up", "GRG down"))]; setorder(dg, cls)
pb <- ggplot(dg, aes(logFC, -log10(P.Value), colour = cls)) + geom_point(size = 0.5, alpha = 0.7, stroke = 0) +
  geom_point(data = dg[cls %in% c("GRG up", "GRG down")], size = 1.2, alpha = 0.95, stroke = 0) +
  geom_point(data = dg[gene %in% hub], shape = 21, size = 2.2, stroke = 0.5, colour = "black", fill = COL_DN) +
  geom_text_repel(data = dg[gene %in% hub], aes(label = gene), size = GT, family = "Arial", fontface = "italic", colour = "black",
                  min.segment.length = 0, segment.size = 0.25, box.padding = 0.5, nudge_x = -0.35, max.overlaps = Inf) +
  scale_colour_manual(values = c(NS = "grey88", `Other DEG` = "grey60", `GRG up` = COL_UP, `GRG down` = COL_DN), name = NULL) +
  labs(x = paste0("log", SUB(2), " fold change (F \u2265 2 vs F0\u20131)"), y = paste0("\u2212log", SUB(10), " <i>P</i>"), title = "GSE135251 (n = 216)") +
  theme(axis.title.x = element_markdown(), axis.title.y = element_markdown(), legend.position = c(0.02, 0.98), legend.justification = c(0, 1),
        legend.background = element_blank())

## c 漏斗（横向条）
fu <- fread(file.path(G10, "01_GRG_funnel.tsv")); ml <- fread(file.path(G10, "02_ML_vote.tsv"))
fcn <- data.table(step = c("GRG", "In training\nmatrix", "\u2229 DEG\n(FDR < 0.05)", "Shared by\n3 cohorts", "4/4 ML"),
                  n = c(nrow(fu), sum(fu$in_train), sum(fu$DEG), sum(fu$candidate), sum(ml$votes == 4)))
fcn[, step := factor(step, levels = rev(step))]
pc <- ggplot(fcn, aes(n, step)) + geom_col(width = 0.62, fill = c(COL_GRG, rep("#4393C3", 4))) +
  geom_text(aes(label = n), hjust = -0.25, size = GT, family = "Arial") + scale_x_continuous(breaks = c(0, 100, 200), expand = expansion(mult = c(0, 0.25))) +
  labs(x = "Genes", y = NULL) + theme(axis.text.y = element_text(lineheight = 0.85))

## d 四算法投票点阵（≥2 票）
v <- ml[votes >= 2][order(-votes, adj.P.Val)]; v[, gene := factor(gene, levels = rev(gene))]
vl <- melt(v, id.vars = c("gene", "votes"), measure.vars = c("LASSO", "SVM_RFE", "XGBoost", "Boruta"))
vl[, variable := factor(sub("_", "-", variable), levels = c("LASSO", "SVM-RFE", "XGBoost", "Boruta"))]
pd <- ggplot(vl, aes(variable, gene)) + geom_point(aes(colour = value, size = value)) +
  geom_text(data = v, aes(x = 4.75, y = gene, label = votes), size = GT, family = "Arial", inherit.aes = FALSE) +
  scale_colour_manual(values = c(`TRUE` = COL_GRG, `FALSE` = "grey85"), guide = "none") + scale_size_manual(values = c(`TRUE` = 2.6, `FALSE` = 1.1), guide = "none") +
  coord_cartesian(xlim = c(0.6, 4.9), clip = "off") + labs(x = NULL, y = NULL) +
  theme(axis.text.y = element_text(face = "italic"), axis.text.x = element_text(angle = 35, hjust = 1), axis.line.x = element_blank(), axis.ticks.x = element_blank())

## e SVM-RFE 曲线
sv <- fread(file.path(U, "U4_svmrfe_curve.tsv")); ko <- sv[which.min(cv_error)]
pe <- ggplot(sv, aes(k, cv_error)) + geom_line(colour = "grey40", linewidth = 0.4) + geom_point(size = 0.8, colour = "grey40") +
  geom_point(data = ko, colour = COL_GRG, size = 2.2) +
  annotate("text", x = ko$k + 3, y = ko$cv_error, label = sprintf("k = %d, error = %.2f", ko$k, ko$cv_error), hjust = 0, size = GT, family = "Arial", colour = COL_GRG) +
  labs(x = "Number of features", y = "10-fold CV error (SVM-RFE)")

## f 森林图：三队列 Hedges g + DL 合并
fr <- fread(file.path(U, "U5_forest_percohort.tsv")); mt <- fread(file.path(U, "U5_forest_meta.tsv"))
fd <- rbind(fr[, .(gene, cohort, g, lo = g - 1.96 * se, hi = g + 1.96 * se, type = "cohort")], mt[, .(gene, cohort = "Pooled (DL)", g, lo, hi, type = "pooled")])
fd[, cohort := factor(cohort, levels = rev(c("Pooled (DL)", "GSE135251", "GSE162694", "GSE130970")))]; fd[, gene := factor(gene, levels = hub)]
ms <- mt[, .(gene = factor(gene, levels = hub), lab = sprintf("I\u00b2 = %.0f%%  %s", I2, stars(P)))]
pf <- ggplot(fd, aes(g, cohort)) + geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3, colour = "grey50") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.25, linewidth = 0.35) +
  geom_point(aes(shape = type, fill = type, size = type), stroke = 0.3) +
  geom_text(data = ms, aes(x = -Inf, y = Inf, label = lab), hjust = -0.05, vjust = -0.3, size = GT, family = "Arial", inherit.aes = FALSE) +
  scale_shape_manual(values = c(cohort = 22, pooled = 23), guide = "none") + scale_fill_manual(values = c(cohort = COL_UP, pooled = COL_DN), guide = "none") +
  scale_size_manual(values = c(cohort = 1.6, pooled = 2.4), guide = "none") +
  facet_wrap(~ gene, nrow = 1) + labs(x = "Standardised mean difference (Hedges g, F\u22652 vs F0\u20131)", y = NULL) +
  coord_cartesian(clip = "off") + theme(strip.text = element_text(face = "bold.italic", margin = margin(b = 9)), panel.spacing = unit(2, "mm"))

rows <- list(list(height = 76, cells = list(list(width = 82, obj = draw_chord, tag = "a"), list(width = 92, obj = pb, tag = "b"))),
             list(height = 62, cells = list(list(width = 50, obj = pc, tag = "c"), list(width = 58, obj = pd, tag = "d"), list(width = 66, obj = pe, tag = "e"))),
             list(height = 52, cells = list(list(width = 174, obj = pf, tag = "f"))))
save_rows(rows, "Fig1")
