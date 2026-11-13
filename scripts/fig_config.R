# 出图配置（改自 Paper_Lactylation_MASH/00_config.R + scripts/P00_prod_theme.R）
# Springer 成图：174 mm 宽、≤234 mm 高、Arial、正文 8 pt、最小 7 pt；面板标签小写 a/b/c（Clin Exp Med 须知）
# 显著性：图上只打星号（全局规则 10）：* P<0.05  ** P<0.01  *** P<0.001  **** P<0.0001  ns；精确 P 在 results/*.tsv
if (nzchar(Sys.getenv("R_LIBS_GLUT"))) .libPaths(Sys.getenv("R_LIBS_GLUT"))
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
if (.Platform$OS.type == "windows") windowsFonts(Arial = windowsFont("Arial"))
PROJ <- Sys.getenv("GLUT_PROJ", "."); R <- file.path(Sys.getenv("GLUT_OUT", PROJ), "results")
FIGP <- file.path(Sys.getenv("GLUT_OUT", PROJ), "figures"); dir.create(FIGP, showWarnings = FALSE)
PAGE_W <- 174; PAGE_HMAX <- 234; MIN_PT <- 7; GT <- 2.5
SEED <- 20261113

COL_SEV <- c(ctrl = "#4393C3", MASL = "#92C5DE", eMASH = "#F4A582", aMASH = "#B2182B")
COL_GRG <- "#B2182B"; COL_REF <- "grey60"; COL_UP <- "#D6604D"; COL_DN <- "#2166AC"
COL_CELLTYPE <- c(Hepatocyte = "#D6604D", Endothelial = "#4393C3", Macrophage = "#F4A582", Tcell_NK = "#5AAE61",
                  Stellate = "#9970AB", Cholangiocyte = "#E7B800", Bcell = "#80B1D3", Plasma = "#2166AC", Cycling = "#999999")

stars <- function(p) as.character(cut(p, c(-Inf, 1e-4, 1e-3, 1e-2, 0.05, Inf), labels = c("****", "***", "**", "*", "ns"), right = FALSE))

theme_prod <- function(base_size = 8) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(axis.line = element_line(linewidth = 0.35, colour = "black"), axis.ticks = element_line(linewidth = 0.35, colour = "black"),
          axis.text = element_text(size = 7, colour = "black"), axis.title = element_text(size = 8),
          legend.text = element_text(size = 7), legend.title = element_text(size = 7.5), legend.key.size = unit(3, "mm"),
          strip.text = element_text(size = 7.5), strip.background = element_blank(),
          plot.title = element_text(size = 8, face = "bold", hjust = 0.5), panel.grid = element_blank(),
          plot.tag = element_text(size = 10, face = "bold", family = "Arial"))
}
theme_set(theme_prod())
TAG_THEME <- theme(plot.tag = element_text(size = 10, face = "bold", family = "Arial"))

check_svg_fonts <- function(stem) {
  s <- paste(readLines(file.path(FIGP, paste0(stem, ".svg")), warn = FALSE), collapse = " ")
  fs <- as.numeric(regmatches(s, gregexpr("(?<=font-size: )[0-9.]+(?=px)", s, perl = TRUE))[[1]])
  out <- data.table(figure = stem, n_text = length(fs), min_pt = min(fs), max_pt = max(fs), n_below_min = sum(fs < MIN_PT - 0.05))
  fwrite(out, file.path(FIGP, paste0(stem, "_fontcheck.tsv")), sep = "\t"); print(out); invisible(out)
}
draw_rows <- function(rows, gap_mm = 2) {
  grid::grid.newpage(); y <- 0
  for (r in rows) { x <- 0
    for (cl in r$cells) {
      vp <- grid::viewport(x = grid::unit(x, "mm"), y = grid::unit(1, "npc") - grid::unit(y, "mm"), width = grid::unit(cl$width, "mm"),
                           height = grid::unit(r$height, "mm"), just = c("left", "top"))
      grid::pushViewport(vp); o <- cl$obj
      if (is.function(o)) { gridGraphics::grid.echo(o, newpage = FALSE)
        if (!is.null(cl$tag)) grid::grid.text(cl$tag, x = grid::unit(0.5, "mm"), y = grid::unit(1, "npc") - grid::unit(0.5, "mm"), just = c("left", "top"),
                                              gp = grid::gpar(fontsize = 10, fontface = "bold", fontfamily = "Arial"))
      } else { if (!is.null(cl$tag)) o <- o + labs(tag = cl$tag); print(o + TAG_THEME, newpage = FALSE) }
      grid::popViewport(); x <- x + cl$width }
    y <- y + r$height + gap_mm }
}
save_rows <- function(rows, stem, gap_mm = 2) {
  h <- sum(vapply(rows, `[[`, numeric(1), "height")) + gap_mm * (length(rows) - 1)
  for (r in rows) stopifnot(sum(vapply(r$cells, `[[`, numeric(1), "width")) <= PAGE_W + 0.01)
  stopifnot(h <= PAGE_HMAX + 0.01)
  f <- file.path(FIGP, stem); w <- PAGE_W / 25.4; hh <- h / 25.4; dr <- function() draw_rows(rows, gap_mm)
  svglite::svglite(paste0(f, ".svg"), width = w, height = hh); dr(); dev.off()
  grDevices::cairo_pdf(paste0(f, ".pdf"), width = w, height = hh, family = "Arial"); dr(); dev.off()
  grDevices::tiff(paste0(f, ".tiff"), width = w, height = hh, units = "in", res = 600, type = "cairo", compression = "lzw"); dr(); dev.off()
  grDevices::png(paste0(f, ".png"), width = w, height = hh, units = "in", res = 300, type = "cairo"); dr(); dev.off()
  cat("saved", stem, sprintf("%.0f x %.0f mm\n", PAGE_W, h)); check_svg_fonts(stem)
}
# 单细胞/空间连续评分统一色标（用户 11-13 从 6 种候选中定：Navy-to-gold）
PAL_SC <- c("#1B2A49", "#2D4B7A", "#4E7FA8", "#9DC3C9", "#F2E3B3", "#F5B14C", "#E2701E")
