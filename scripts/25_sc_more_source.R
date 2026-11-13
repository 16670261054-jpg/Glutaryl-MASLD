# 单细胞补充图第三批源数据（预注册变更记录 SC1–SC5；只报告）
if (nzchar(Sys.getenv("R_LIBS_GLUT"))) .libPaths(Sys.getenv("R_LIBS_GLUT"))
suppressPackageStartupMessages({ library(data.table); library(Matrix); library(limma); library(edgeR) })
PROJ <- Sys.getenv("GLUT_PROJ", "."); RES <- file.path(Sys.getenv("GLUT_OUT", PROJ), "results/25_sc_more"); dir.create(RES, showWarnings = FALSE)
LAC <- Sys.getenv("GLUT_SN", "data/snRNA"); SM <- Sys.getenv("GLUT_SM", "data/GSE289173_series_matrix.txt.gz")
write_tab <- function(df, name) { p <- file.path(RES, name); tmp <- paste0(p, ".tmp"); fwrite(as.data.frame(df), tmp, sep = "\t"); file.rename(tmp, p); message("written: ", p) }
grg <- readLines(file.path(PROJ, "geneset/GRG_full.txt")); hub <- readLines(file.path(Sys.getenv("GLUT_OUT", PROJ), "results/10_gates/02_hub.txt"))
hsig <- fread(file.path(Sys.getenv("GLUT_OUT", PROJ), "results/13_hep_mash/13_Shep_GRG_genes.tsv"))[BH < 0.05, gene]
KEY <- unique(c("GCDH", "AASS", "SIRT5", hub, hsig))
mito_ref <- readLines(file.path(Sys.getenv("GLUT_OUT", PROJ), "results/16_batch2/mito_reference_200.txt"))
PP <- c("ASS1", "PCK1", "CPS1", "HAL", "SDS", "ARG1", "ASGR2"); PC <- c("GLUL", "CYP2E1", "CYP1A2", "CYP3A4", "OAT", "RHBG", "ADH4")
ov <- intersect(c(PP, PC), grg); PPx <- setdiff(PP, grg); PCx <- setdiff(PC, grg); grg_x <- setdiff(grg, c(PP, PC))
write_tab(data.table(item = c("zonation markers overlapping GRG (removed)", "periportal used", "pericentral used", "GRG genes after removing markers"),
                     value = c(paste(ov, collapse = ","), paste(PPx, collapse = ","), paste(PCx, collapse = ","), length(grg_x))), "00_marker_overlap.tsv")

## ---- 供者伪 bulk（SC2 / SC3 / SC5）----
ref <- readRDS(file.path(LAC, "sc_data/deconv_ref_pseudobulk.rds"))
lab <- ifelse(grepl("^Hep[0-9]+$", ref$label), "Hepatocyte", ref$label); key <- paste(ref$sample, lab, sep = "|")
agg <- t(rowsum(t(ref$mat), key)); nc <- tapply(ref$ncell, key, sum)[colnames(agg)]
sm <- readLines(gzfile(SM)); ch <- sm[grepl("^!Sample_characteristics_ch1", sm)]
row_of <- function(pat) gsub('"', '', strsplit(grep(pat, sm, value = TRUE)[1], "\t")[[1]][-1])
getch <- function(k) { ln <- ch[grepl(k, ch)][1]; trimws(sub(".*:", "", gsub('"', '', strsplit(ln, "\t")[[1]][-1]))) }
meta <- data.table(sample = row_of("^!Sample_title"), group = getch("disease group"), sex = getch("Sex: "))
meta[, severity := c(ctrl = 0, MASL = 1, eMASH = 2, aMASH = 3)[group]]; meta[, batch := ifelse(grepl("^FL", sample), "FL", "ctrl_series")]
cts <- sort(unique(sub("^[^|]+\\|", "", colnames(agg)))); out_g <- list(); out_c <- list(); hepE <- NULL
for (ct in cts) {
  k <- sub("^[^|]+\\|", "", colnames(agg)) == ct; cn <- round(agg[, k, drop = FALSE]); colnames(cn) <- sub("\\|.*", "", colnames(cn))
  md <- meta[match(colnames(cn), sample)]; md[, n_nuc := as.numeric(nc[k])]; kk <- md$n_nuc >= 20; cn <- cn[, kk, drop = FALSE]; md <- md[kk]
  if (ncol(cn) < 12 || uniqueN(md$batch) < 2) { out_c[[ct]] <- data.table(celltype = ct, n_donors = ncol(cn)); next }
  des <- model.matrix(~ factor(sex) + factor(batch) + severity, md)
  dge <- DGEList(cn); dge <- dge[filterByExpr(dge, des), , keep.lib.sizes = FALSE]; dge <- calcNormFactors(dge); v <- voom(dge, des)
  tt <- topTable(eBayes(lmFit(v, des)), coef = "severity", n = Inf, sort.by = "none")
  out_g[[ct]] <- data.table(celltype = ct, gene = KEY, logFC = tt[KEY, "logFC"], BH = tt[KEY, "adj.P.Val"], tested = KEY %in% rownames(tt))
  cm <- camera(v, list(GRG = which(rownames(v) %in% grg)), des, contrast = ncol(des))
  out_c[[ct]] <- data.table(celltype = ct, n_donors = ncol(cn), n_genes = cm$NGenes, direction = cm$Direction, P = cm$PValue)
  if (ct == "Hepatocyte") hepE <- v$E
}
G <- rbindlist(out_g); C <- rbindlist(out_c, fill = TRUE); C[!is.na(P), BH := p.adjust(P, "BH")]
write_tab(G, "SC2_celltype_keygene_logFC.tsv"); write_tab(C, "SC3_celltype_GRG_camera.tsv")
k5 <- intersect(KEY, rownames(hepE)); cr <- cor(t(hepE[k5, ]), method = "spearman")
write_tab(data.table(gene = rownames(cr), cr), "SC5_hep_keygene_spearman.tsv")

## ---- SC1 全图谱（大对象，最后载入）----
suppressPackageStartupMessages(library(Seurat))
so <- readRDS(file.path(LAC, "sc_data/GSE289173_annotated.rds"))
um <- Embeddings(so, "umap"); X <- GetAssayData(so, layer = "data")
a <- data.table(cell = colnames(so), umap_1 = um[, 1], umap_2 = um[, 2], celltype = so$celltype, group = so$group, sample = so$sample,
                GRG = Matrix::colMeans(X[intersect(grg, rownames(X)), ]))
for (g in c("GCDH", "AASS", "ACADSB", "SARDH")) a[[g]] <- as.numeric(X[g, ])
rm(so, X); invisible(gc())
write_tab(a, "SC1_atlas_cells.tsv.gz")
write_tab(a[, .(GRG_mean = mean(GRG), n = .N), by = .(sample, group, celltype)][n >= 20], "SC1_donor_celltype_GRG.tsv")
cat("DONE\n"); print(C)
