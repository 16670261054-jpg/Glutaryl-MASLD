# Fig. 2（升级版）snRNA 肝细胞内 GRG 程序（只读盘画图）；v1 备份 _fig2_v1_backup.R
# 注：不重复乳酸化稿的"状态着色 UMAP / 构成流线图"（同一数据，避免自我重复），改为逐细胞 GRG 评分 UMAP + Ro/e
source(file.path(Sys.getenv("GLUT_PROJ", "."), "scripts/fig_config.R"))
suppressPackageStartupMessages({ library(ggdist); library(ggtext); library(ggrepel) })
SUB <- function(x) paste0("<sub><span style='font-size:7pt'>", x, "</span></sub>")
H <- file.path(R, "13_hep_mash"); B1 <- file.path(R, "14_batch1"); U <- file.path(R, "23_upgrade")
S <- fread(file.path(H, "13_Shep_hepatocyte_program.tsv"))

## a GRG 评分 UMAP（按组分面）：六边形分箱内逐核评分的均值（全部肝细胞核，不抽样）；色标 inferno，与乳酸化稿评分图同一配色体系
suppressPackageStartupMessages(library(hexbin))
cs <- fread(file.path(U, "U1_hep_cells.tsv.gz")); cs[, group := factor(group, levels = names(COL_SEV))]
lim <- quantile(cs$GRG, c(0.02, 0.98), na.rm = TRUE)
xr <- quantile(cs$umap_1, c(0.005, 0.995)); yr <- quantile(cs$umap_2, c(0.005, 0.995))
lab <- cs[, { gx <- cut(umap_1, 25); gy <- cut(umap_2, 25); k <- paste(gx, gy); tk <- names(sort(table(k), decreasing = TRUE))[1]
  .(x = median(umap_1[k == tk]), y = median(umap_2[k == tk])) }, by = hepstate]
lab <- lab[x %between% xr & y %between% yr]  # 标签按全部核定位；落在裁剪范围外的状态（Hep6 卫星簇）不标
cs <- cs[umap_1 %between% xr & umap_2 %between% yr]  # 只裁掉坐标最外 0.5% 的离散核，使主体云团占满面板
ncell <- cs[, .N, by = group][, setNames(sprintf("%s (%s nuclei)", group, trimws(format(N, big.mark = ","))), group)]
pa <- ggplot(cs, aes(umap_1, umap_2, z = GRG)) +
  stat_summary_hex(fun = mean, bins = 70, colour = NA) +
  geom_text_repel(data = lab, aes(x, y, label = hepstate), colour = "white", bg.colour = "black", bg.r = 0.12, size = GT, family = "Arial",
                  fontface = "bold", inherit.aes = FALSE, seed = SEED, box.padding = 0.1, min.segment.length = Inf) +
  scale_fill_gradientn(colours = PAL_SC, limits = lim, oob = scales::squish, name = "Mean GRG
score") +
  facet_wrap(~ group, nrow = 1, labeller = as_labeller(ncell)) + coord_equal() + labs(x = "UMAP 1", y = "UMAP 2") +
  theme(strip.text = element_text(face = "bold"), legend.key.height = unit(5, "mm"), legend.key.width = unit(2.2, "mm"),
        axis.text = element_blank(), axis.ticks = element_blank(), panel.background = element_rect(fill = "white", colour = NA))

## b 每供者肝细胞 GRG 评分（伪 bulk）
dn <- fread(file.path(R, "20_figsrc/F2_hep_donor_scores.tsv")); dn[, group := factor(group, levels = names(COL_SEV))]
set.seed(SEED); dn[, xj := as.numeric(group) - 0.2 + runif(.N, -0.06, 0.06)]
pb <- ggplot(dn, aes(group, GRG_score, fill = group)) +
  stat_slab(scale = 0.55, alpha = 0.85, colour = NA, position = position_nudge(x = 0.12), normalize = "xy") +
  geom_boxplot(width = 0.14, outlier.shape = NA, linewidth = 0.3, position = position_nudge(x = -0.02), alpha = 0.95) +
  geom_point(aes(x = xj), shape = 21, size = 1, stroke = 0.2, alpha = 0.85) +
  annotate("text", x = 2.5, y = max(dn$GRG_score) * 1.12, label = paste0("camera ", stars(S[set == "GRG_full", PValue])), size = GT, family = "Arial") +
  scale_fill_manual(values = COL_SEV, guide = "none") + labs(x = NULL, y = "Hepatocyte GRG score (z)", title = sprintf("Donor pseudobulk (n = %d)", nrow(dn)))

## c 匹配零分布
dr <- fread(file.path(H, "13_Shep_null_draws.tsv"))[set == "GRG_full"]; ob <- fread(file.path(H, "13_Shep_obs_meant.tsv"))[set == "GRG_full", mean_t]
dr[, null := factor(null, levels = c("N1", "N3"), labels = c("Expression-matched genes", "Mitochondrial genes (non-GRG)"))]
st <- data.table(null = factor(levels(dr$null), levels = levels(dr$null)), lab = stars(c(S[set == "GRG_full", null_p_N1], S[set == "GRG_full", null_p_N3])))
pc <- ggplot(dr, aes(mean_t)) + stat_slab(aes(fill = after_stat(x < ob)), slab_colour = NA, scale = 0.9, normalize = "panels") +
  stat_pointinterval(aes(y = 0), point_size = 1, interval_size_range = c(0.4, 0.9), colour = "grey30") +
  geom_vline(xintercept = ob, colour = COL_GRG, linewidth = 0.6) +
  geom_text(data = st, aes(x = ob, y = 0.95, label = paste("GRG", lab)), hjust = -0.1, size = GT, colour = COL_GRG, family = "Arial") +
  scale_fill_manual(values = c(`TRUE` = COL_GRG, `FALSE` = "grey75"), guide = "none") +
  facet_wrap(~ null, nrow = 1) + scale_y_continuous(breaks = NULL) +
  labs(x = "Mean moderated t (severity), 1,000 random sets", y = "Density") + theme(axis.line.y = element_blank())

## d 亚状态 Ro/e（单元格数字）+ 状态内 GRG camera
ro <- fread(file.path(U, "U2_hepstate_Roe.tsv")); ro[, group := factor(group, levels = names(COL_SEV))]
a3 <- fread(file.path(B1, "A3_hepstate_GRG_camera.tsv"))
ro <- rbind(ro[, .(hepstate, col = group, val = Roe, lab = sprintf("%.2f", Roe), kind = "Roe")],
            a3[, .(hepstate = state, col = "GRG\nwithin state", val = NA_real_, lab = fifelse(is.na(BH), "n.t.", paste0(fifelse(direction == "Down", "\u2193", "\u2191"), stars(BH))), kind = "cam")])
ro[, col := factor(col, levels = c(names(COL_SEV), "GRG\nwithin state"))]; ro[, hepstate := factor(hepstate, levels = rev(sort(unique(hepstate))))]
pd <- ggplot(ro, aes(col, hepstate)) + geom_tile(data = ro[kind == "Roe"], aes(fill = pmin(val, 3)), colour = "white", linewidth = 0.5) +
  geom_tile(data = ro[kind == "cam"], fill = "grey96", colour = "white", linewidth = 0.5) +
  geom_text(aes(label = lab), size = GT, family = "Arial") +
  scale_fill_gradient2(low = COL_DN, mid = "white", high = COL_UP, midpoint = 1, limits = c(0, 3), name = "Ro/e") +
  labs(x = NULL, y = "Hepatocyte state") + theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_text(lineheight = 0.85))

## e 肝细胞内前 15 个 GRG 基因
gn <- fread(file.path(H, "13_Shep_GRG_genes.tsv"))[order(P)][1:15]; gn[, gene := factor(gene, levels = rev(gene))]
pe <- ggplot(gn, aes(logFC, gene)) + geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey50") +
  geom_segment(aes(x = 0, xend = logFC, yend = gene), linewidth = 0.4, colour = "grey55") +
  geom_point(aes(size = -log10(P), fill = logFC), shape = 21, stroke = 0.3) +
  geom_text(aes(x = min(logFC) * 1.28, label = stars(BH)), size = GT, family = "Arial", hjust = 0) +
  scale_fill_gradient(low = COL_DN, high = "#D1E5F0", guide = "none") + scale_size_continuous(range = c(1.2, 3), name = "\u2212log10 P") +
  scale_x_continuous(limits = c(min(gn$logFC) * 1.32, 0.04)) + labs(x = paste0("log", SUB(2), " FC per severity step"), y = NULL) +
  theme(axis.title.x = element_markdown(), axis.text.y = element_text(face = "italic"), legend.position = "right")

## f 关键基因的细胞类型表达
lc <- fread(file.path(B1, "A2_celltype_localization.tsv")); lc[, z := (median_CPM - mean(median_CPM)) / sd(median_CPM), by = gene]
lc[, gene := factor(gene, levels = rev(c("GCDH", "AASS", "SIRT5", "KAT2A", "ACADSB", "SARDH", "ACOT13", "ACAA2")))]
lc[, celltype := factor(celltype, levels = names(COL_CELLTYPE))]
pf <- ggplot(lc, aes(celltype, gene)) + geom_point(aes(size = log10(median_CPM + 1), fill = z), shape = 21, stroke = 0.25) +
  scale_fill_gradient2(low = COL_DN, mid = "white", high = COL_UP, name = "z") + scale_size_continuous(range = c(0.5, 3.6), name = "log10 CPM") +
  labs(x = NULL, y = NULL) + theme(axis.text.x = element_text(angle = 40, hjust = 1), axis.text.y = element_text(face = "italic"))

## g 四个集合：观测 vs 两种零分布中位数
LAB <- c(GRG_full = "GRG (all)", GRG_substrate = "Substrates", GRG_substrate_both = "Substrates\n(both datasets)", GRG_enzyme_axis = "Enzyme axis")
cs4 <- melt(S[, .(set, Observed = mean_t, `Expression-matched` = null_median_N1, Mitochondrial = null_median_N3)], id.vars = "set")
cs4[, set := factor(LAB[set], levels = rev(LAB))]
pg <- ggplot(cs4, aes(value, set)) + geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey50") +
  geom_line(aes(group = set), linewidth = 0.4, colour = "grey70") + geom_point(aes(fill = variable), shape = 21, size = 2, stroke = 0.3) +
  scale_fill_manual(values = c(Observed = COL_GRG, `Expression-matched` = "grey80", Mitochondrial = "grey45"), name = NULL) +
  labs(x = "Mean moderated t", y = NULL) + theme(legend.position = "bottom", legend.direction = "vertical", axis.text.y = element_text(lineheight = 0.85))

rows <- list(list(height = 54, cells = list(list(width = 174, obj = pa, tag = "a"))),
             list(height = 50, cells = list(list(width = 60, obj = pb, tag = "b"), list(width = 114, obj = pc, tag = "c"))),
             list(height = 58, cells = list(list(width = 82, obj = pd, tag = "d"), list(width = 92, obj = pe, tag = "e"))),
             list(height = 60, cells = list(list(width = 108, obj = pf, tag = "f"), list(width = 66, obj = pg, tag = "g"))))
save_rows(rows, "Fig2")
