# Fig. 3 Visium 空间分带（只读盘画图；n = 4 供者，探索性）
source(file.path(Sys.getenv("GLUT_PROJ", "."), "scripts/fig_config.R"))
sp <- fread(file.path(R, "20_figsrc/F3_spatial_spots.tsv")); zb <- fread(file.path(R, "20_figsrc/F3_zonation_deciles.tsv"))
zs <- fread(file.path(R, "16_batch2/B_zonation_per_section.tsv"))
DON <- c(JBO14 = "H35", JBO15 = "H35", JBO18 = "H36", JBO19 = "H37", JBO22 = "H38")

## a/b 代表切片：H&E 低分辨率底图 + spot 评分（坐标 = fullres 像素 × tissue_lowres_scalef）
VD <- "D:/cc/生信代码/Paper_Lactylation_MASH/01_数据/GSE192741_visium"
GSM <- c(JBO19 = "GSM5764427", JBO18 = "GSM5764426")
he <- function(s) { f <- list.files(VD, pattern = paste0(GSM[[s]], "_tissue_lowres_image"), full.names = TRUE)
  tf <- tempfile(fileext = ".png"); con <- gzfile(f, "rb"); writeBin(readBin(con, "raw", 5e7), tf); close(con)
  sf <- jsonlite::fromJSON(readLines(gzfile(list.files(VD, pattern = paste0(GSM[[s]], "_scalefactors"), full.names = TRUE)), warn = FALSE))
  list(img = png::readPNG(tf), sf = sf$tissue_lowres_scalef) }
rng <- quantile(sp[sample %in% c("JBO19", "JBO18"), GRG_full], c(0.02, 0.98))
PAL_G <- c("#FFF5EB", "#FDBE85", "#FD8D3C", "#D94701", "#7F2704"); PAL_Z <- c("#2166AC", "#F7F7F7", "#B2182B")
smap <- function(s, var, ttl, lims, pal, nm, show_he = TRUE) {
  h <- he(s); d <- sp[sample == s]; d[, `:=`(x = pxl_col * h$sf, y = pxl_row * h$sf)]
  pad <- 25; xr <- range(d$x) + c(-pad, pad); yr <- range(d$y) + c(-pad, pad)
  W <- dim(h$img)[2]; Hh <- dim(h$img)[1]
  p <- ggplot(d, aes(x, y))
  if (show_he) p <- p + annotation_raster(h$img, xmin = 0, xmax = W, ymin = -Hh, ymax = 0)
  p + geom_point(aes(x = x, y = -y, fill = get(var)), shape = 21, size = 0.9, stroke = 0, alpha = if (show_he) 0.85 else 1) +
    scale_fill_gradientn(colours = pal, limits = lims, oob = scales::squish, name = nm) +
    coord_equal(xlim = xr, ylim = -rev(yr), expand = FALSE) + labs(title = ttl) + theme_void(base_family = "Arial") +
    theme(plot.title = element_text(size = 7.5, hjust = 0.5, lineheight = 0.9), legend.text = element_text(size = 7), legend.title = element_text(size = 7),
          legend.position = "right", legend.key.width = unit(2, "mm"), legend.key.height = unit(3.5, "mm"), legend.box.spacing = unit(0.5, "mm"), plot.tag = element_text(size = 10, face = "bold", family = "Arial"))
}
he_only <- function(s, ttl) { h <- he(s); d <- sp[sample == s]; d[, `:=`(x = pxl_col * h$sf, y = pxl_row * h$sf)]
  pad <- 25; xr <- range(d$x) + c(-pad, pad); yr <- range(d$y) + c(-pad, pad); W <- dim(h$img)[2]; Hh <- dim(h$img)[1]
  ggplot() + annotation_raster(h$img, xmin = 0, xmax = W, ymin = -Hh, ymax = 0) + geom_blank(data = d, aes(x, -y)) +
    coord_equal(xlim = xr, ylim = -rev(yr), expand = FALSE) + labs(title = ttl) + theme_void(base_family = "Arial") +
    theme(plot.title = element_text(size = 7.5, hjust = 0.5, lineheight = 0.9), plot.tag = element_text(size = 10, face = "bold", family = "Arial")) }
pa <- (he_only("JBO19", "JBO19, steatotic
H&E") + labs(tag = "a")) + smap("JBO19", "GRG_full", "GRG score", rng, PAL_G, "GRG", FALSE) +
      smap("JBO19", "zonation", "Zonation
(portal → central)", c(0, 1), PAL_Z, "Zon.", FALSE) 
pb <- (he_only("JBO18", "JBO18, healthy
H&E") + labs(tag = "b")) + smap("JBO18", "GRG_full", "GRG score", rng, PAL_G, "GRG", FALSE) +
      smap("JBO18", "zonation", "Zonation
(portal → central)", c(0, 1), PAL_Z, "Zon.", FALSE)

## c 按 zonation 十分位的评分（切片内 z）
zb[, variable := factor(variable, levels = c("GRG_full", "MITO_ref200"), labels = c("GRG", "Mitochondrial reference (200 genes)"))]
zb[, section := paste0(sample, " (", diet, ")")]
pc <- ggplot(zb, aes(zdec, z_within, colour = variable, group = interaction(sample, variable))) +
  geom_line(linewidth = 0.4, alpha = 0.8) + geom_point(size = 0.8) +
  scale_colour_manual(values = c(GRG = COL_GRG, `Mitochondrial reference (200 genes)` = COL_REF), name = NULL) +
  scale_x_continuous(breaks = c(1, 5, 10), labels = c("1\nportal", "5", "10\ncentral")) +
  labs(x = "Zonation decile", y = "Median spot score (z within section)") +
  theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1), legend.background = element_blank())

## d 每切片分带相关（控制肝细胞占比的偏相关）
dz <- zs[feature %in% c("GRG_full", "GRG_enzyme_axis", "MITO_ref200")]
dz[, feature := factor(feature, levels = c("MITO_ref200", "GRG_enzyme_axis", "GRG_full"), labels = c("Mitochondrial\nreference", "Enzyme axis", "GRG"))]
dz[, section := factor(paste0(sample, " (", DON[sample], ", ", diet, ")"))]
pd <- ggplot(dz, aes(rho_partial, section)) + geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey50") +
  geom_line(aes(group = section), colour = "grey75", linewidth = 0.4) +
  geom_point(aes(fill = feature), shape = 21, size = 2, stroke = 0.3) +
  scale_fill_manual(values = c(GRG = COL_GRG, `Enzyme axis` = "#F4A582", `Mitochondrial\nreference` = COL_REF), name = NULL) +
  labs(x = "Partial Spearman \u03c1 with zonation\n(adjusted for hepatocyte fraction)", y = NULL) + theme(legend.position = "right")

rows <- list(list(height = 50, cells = list(list(width = 174, obj = pa))), list(height = 50, cells = list(list(width = 174, obj = pb))),
             list(height = 58, cells = list(list(width = 80, obj = pc, tag = "c"), list(width = 94, obj = pd, tag = "d"))))
save_rows(rows, "Fig3")
