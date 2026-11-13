# 单细胞补充图 S12–S14（只读盘画图；预注册变更记录 SC1–SC5）
source(file.path(Sys.getenv("GLUT_PROJ", "."), "scripts/fig_config.R"))
suppressPackageStartupMessages({ library(ggrepel); library(ggdist); library(ggtext) })
S <- file.path(R, "25_sc_more"); U <- file.path(R, "23_upgrade")
lab_pos <- function(d, by) d[, { gx <- cut(umap_1, 25); gy <- cut(umap_2, 25); k <- paste(gx, gy); tk <- names(sort(table(k), decreasing = TRUE))[1]
  .(x = median(umap_1[k == tk]), y = median(umap_2[k == tk])) }, by = by]

## ---------- S12 全图谱：细胞类型、GRG 评分、各细胞类型的 GRG 变化 ----------
a <- fread(file.path(S, "SC1_atlas_cells.tsv.gz"))[!is.na(GRG)]; a[, group := factor(group, levels = names(COL_SEV))]
a[, celltype := factor(celltype, levels = names(COL_CELLTYPE))]
set.seed(SEED); ap <- a[sample(.N, min(.N, 60000))]
lb <- lab_pos(ap, "celltype")
s12a <- ggplot(ap, aes(umap_1, umap_2, colour = celltype)) + geom_point(size = 0.08, stroke = 0, alpha = 0.6) +
  geom_text_repel(data = lb, aes(x, y, label = celltype), colour = "black", size = GT, family = "Arial", inherit.aes = FALSE, seed = SEED,
                  box.padding = 0.15, min.segment.length = 0.3, segment.size = 0.2) +
  scale_colour_manual(values = COL_CELLTYPE, guide = "none") + coord_equal() + labs(x = "UMAP 1", y = "UMAP 2", title = sprintf("GSE289173 (%s nuclei)", format(nrow(a), big.mark = ","))) +
  theme(axis.text = element_blank(), axis.ticks = element_blank())
setorder(ap, GRG); lim <- quantile(a$GRG, c(0.02, 0.98))
s12b <- ggplot(ap, aes(umap_1, umap_2, colour = GRG)) + geom_point(size = 0.1, stroke = 0) +
  scale_colour_gradientn(colours = PAL_SC, limits = lim, oob = scales::squish, name = "GRG\nscore") + coord_equal() +
  labs(x = "UMAP 1", y = "UMAP 2", title = "GRG score per nucleus") +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), legend.key.height = unit(5, "mm"), legend.key.width = unit(2.2, "mm"))
dc <- fread(file.path(S, "SC1_donor_celltype_GRG.tsv")); dc[, group := factor(group, levels = names(COL_SEV))]
cam <- fread(file.path(S, "SC3_celltype_GRG_camera.tsv"))
dc[, celltype := factor(celltype, levels = cam[order(P), celltype])]
cl <- cam[!is.na(P), .(celltype = factor(celltype, levels = levels(dc$celltype)), lab = paste0(fifelse(direction == "Down", "\u2193", "\u2191"), stars(BH)))]
s12c <- ggplot(dc, aes(group, GRG_mean, fill = group)) + geom_boxplot(width = 0.6, outlier.shape = NA, linewidth = 0.25) +
  geom_point(position = position_jitter(width = 0.12, seed = SEED), shape = 21, size = 0.6, stroke = 0.15) +
  geom_text(data = cl, aes(x = 2.5, y = Inf, label = lab), vjust = 1.3, size = GT, family = "Arial", inherit.aes = FALSE) +
  scale_fill_manual(values = COL_SEV, guide = "none") + scale_y_continuous(expand = expansion(mult = c(0.05, 0.22))) +
  facet_wrap(~ celltype, nrow = 2, scales = "free_y") + labs(x = NULL, y = "Donor mean GRG score",
  title = "Per-cell-type pseudobulk camera (BH across 9 cell types)") +
  theme(axis.text.x = element_text(angle = 40, hjust = 1), strip.text = element_text(face = "bold"))
save_rows(list(list(height = 72, cells = list(list(width = 88, obj = s12a, tag = "a"), list(width = 86, obj = s12b, tag = "b"))),
               list(height = 96, cells = list(list(width = 174, obj = s12c, tag = "c")))), "FigS12")

## ---------- S13 关键基因：细胞类型特异性变化 + 逐核表达 + 供者分布 ----------
g2 <- fread(file.path(S, "SC2_celltype_keygene_logFC.tsv"))[tested == TRUE]
ord_g <- g2[celltype == "Hepatocyte"][order(logFC), gene]
g2[, `:=`(gene = factor(gene, levels = rev(ord_g)), celltype = factor(celltype, levels = names(COL_CELLTYPE)), lab = fifelse(BH < 0.05, stars(BH), ""))]
s13a <- ggplot(g2, aes(celltype, gene, fill = pmax(pmin(logFC, 0.5), -0.5))) + geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = lab), size = GT, family = "Arial") +
  scale_fill_gradient2(low = COL_DN, mid = "white", high = COL_UP, limits = c(-0.5, 0.5), name = "log2FC per\nseverity step") +
  scale_x_discrete(drop = FALSE) + labs(x = NULL, y = NULL, title = "Donor pseudobulk per cell type") +
  theme(axis.text.x = element_text(angle = 40, hjust = 1), axis.text.y = element_text(face = "italic"), axis.line = element_blank(), axis.ticks = element_blank())
suppressPackageStartupMessages(library(hexbin))
hc <- fread(file.path(U, "U1_hep_cells.tsv.gz"))
xr <- quantile(hc$umap_1, c(0.005, 0.995)); yr <- quantile(hc$umap_2, c(0.005, 0.995))
hc <- hc[group %in% c("ctrl", "aMASH") & umap_1 %between% xr & umap_2 %between% yr]
hl <- melt(hc, id.vars = c("umap_1", "umap_2", "group"), measure.vars = c("GCDH", "AASS", "ACADSB"), variable.name = "gene")
hl[, group := factor(group, levels = c("ctrl", "aMASH"))]
hl[, v := value / quantile(value[value > 0], 0.98), by = gene]  # 每个基因按非零值 98 分位缩放到同一色标，仅作展示
s13b <- ggplot(hl, aes(umap_1, umap_2, z = v)) + stat_summary_hex(fun = mean, bins = 55, colour = NA) +
  scale_fill_gradientn(colours = PAL_SC, limits = c(0, 0.6), oob = scales::squish, name = "Mean scaled
expression") +
  facet_grid(group ~ gene) + coord_equal() + labs(x = NULL, y = NULL) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), strip.text.x = element_text(face = "bold.italic"),
        strip.text.y = element_text(face = "bold"), legend.key.height = unit(4, "mm"), legend.key.width = unit(2, "mm"))
dn <- fread(file.path(R, "20_figsrc/F2_hep_donor_scores.tsv")); dn[, group := factor(group, levels = names(COL_SEV))]
bh <- fread(file.path(R, "13_hep_mash/13_Shep_GRG_genes.tsv"))
gl <- c("GCDH", "AASS", "ACADSB", "SARDH", "ACOT13", "ACAA2")
dl <- melt(dn, id.vars = c("sample", "group"), measure.vars = gl, variable.name = "gene")
st <- bh[gene %in% gl, .(gene = factor(gene, levels = gl), lab = stars(BH))]
s13c <- ggplot(dl, aes(group, value, fill = group)) + stat_slab(scale = 0.5, alpha = 0.85, colour = NA, position = position_nudge(x = 0.12), normalize = "xy") +
  geom_boxplot(width = 0.14, outlier.shape = NA, linewidth = 0.25, position = position_nudge(x = -0.04)) +
  geom_point(position = position_nudge(x = -0.2), shape = 21, size = 0.6, stroke = 0.15) +
  geom_text(data = st, aes(x = 2.5, y = Inf, label = lab), vjust = 1.3, size = GT, family = "Arial", inherit.aes = FALSE) +
  scale_fill_manual(values = COL_SEV, guide = "none") + scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  facet_wrap(~ gene, nrow = 1, scales = "free_y") + labs(x = NULL, y = "Hepatocyte pseudobulk log2 CPM") +
  theme(axis.text.x = element_text(angle = 40, hjust = 1), strip.text = element_text(face = "bold.italic"))
save_rows(list(list(height = 70, cells = list(list(width = 72, obj = s13a, tag = "a"), list(width = 102, obj = s13b, tag = "b"))),
               list(height = 56, cells = list(list(width = 174, obj = s13c, tag = "c")))), "FigS13")

## ---------- S14 肝细胞内关键 GRG 基因的供者级共调控 ----------
cr <- fread(file.path(S, "SC5_hep_keygene_spearman.tsv")); M <- as.matrix(cr[, -1]); rownames(M) <- cr$gene
ho <- hclust(as.dist(1 - M))$order; cl5 <- as.data.table(as.table(M)); setnames(cl5, c("g1", "g2", "rho"))
cl5[, `:=`(g1 = factor(g1, levels = rownames(M)[ho]), g2 = factor(g2, levels = rev(rownames(M)[ho])))]
cl5[, up := as.integer(g1) > (nrow(M) + 1 - as.integer(g2))]
s14 <- ggplot(cl5, aes(g1, g2)) + geom_tile(data = cl5[up == FALSE], aes(fill = rho), colour = "white", linewidth = 0.4) +
  geom_point(data = cl5[up == TRUE], aes(size = abs(rho), colour = rho)) +
  geom_text(data = cl5[up == FALSE & as.character(g1) != as.character(g2)], aes(label = sprintf("%.2f", rho)), size = GT, family = "Arial") +
  scale_fill_gradient2(low = COL_DN, mid = "white", high = COL_UP, limits = c(-1, 1), name = "Spearman \u03c1") +
  scale_colour_gradient2(low = COL_DN, mid = "white", high = COL_UP, limits = c(-1, 1), guide = "none") + scale_size_continuous(range = c(0.5, 4), guide = "none") +
  labs(x = NULL, y = NULL, title = sprintf("Hepatocyte pseudobulk, %d donors", 48)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"), axis.text.y = element_text(face = "italic"), axis.line = element_blank(), axis.ticks = element_blank())
save_rows(list(list(height = 120, cells = list(list(width = 150, obj = s14)))), "FigS14")
