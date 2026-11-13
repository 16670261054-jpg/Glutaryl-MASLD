# D3 FXR 激动剂干预检验（预注册变更记录 D3，看数据前写死）
if (nzchar(Sys.getenv("R_LIBS_GLUT"))) .libPaths(Sys.getenv("R_LIBS_GLUT"))
suppressPackageStartupMessages({ library(data.table); library(limma) })
PROJ <- Sys.getenv("GLUT_PROJ", "."); RES <- file.path(Sys.getenv("GLUT_OUT", PROJ), "results/18_fxr"); dir.create(RES, showWarnings = FALSE)
DL <- file.path(PROJ, "data/fxr"); GS <- file.path(PROJ, "geneset")
write_tab <- function(df, name) { p <- file.path(RES, name); tmp <- paste0(p, ".tmp")
  fwrite(as.data.frame(df), tmp, sep = "\t", quote = FALSE); file.rename(tmp, p); message("written: ", p) }
grg <- readLines(file.path(GS, "GRG_full.txt")); POS <- c("NR0B2", "SLC51B", "ABCB11")
read_sm <- function(f) { sm <- readLines(gzfile(f)); i0 <- grep("^!series_matrix_table_begin", sm); i1 <- grep("^!series_matrix_table_end", sm)
  ex <- fread(text = sm[(i0 + 1):(i1 - 1)]); id <- as.character(ex[[1]]); ex <- as.matrix(ex[, -1]); rownames(ex) <- id
  ch <- sm[grepl("^!Sample_characteristics_ch1", sm)]
  getch <- function(k) { ln <- ch[grepl(k, ch)][1]; trimws(sub("^[^:]*:", "", gsub('"', '', strsplit(ln, "\t")[[1]][-1]))) }
  ttl <- gsub('"', '', strsplit(grep("^!Sample_title", sm, value = TRUE), "\t")[[1]][-1])
  list(ex = ex, getch = getch, title = ttl) }
annot <- function(f, ex) { ln <- readLines(gzfile(f)); ln <- ln[!grepl("^[!#^]", ln)]; a <- fread(text = ln, sep = "\t", quote = "", select = c("ID", "Gene symbol"))
  setnames(a, c("ID", "sym")); a <- a[sym != "" & !grepl("///", sym)][as.character(ID) %in% rownames(ex)]
  e <- ex[as.character(a$ID), ]; o <- order(-rowMeans(e)); e <- e[o, ]; s <- a$sym[o]; k <- !duplicated(s); e <- e[k, ]; rownames(e) <- s[k]; e }
one <- function(E, des, coef, label) {
  fit <- eBayes(lmFit(E, des)); tt <- topTable(fit, coef = coef, n = Inf, sort.by = "none")
  cm <- camera(E, list(GRG_full = which(rownames(E) %in% grg)), des, contrast = coef)
  pc <- data.table(gene = POS, logFC = tt[POS, "logFC"], P = tt[POS, "P.Value"])
  pc_ok <- sum(pc$logFC > 0 & pc$P < 0.05, na.rm = TRUE) >= 2
  list(sum = data.table(dataset = label, n = ncol(E), n_GRG = cm$NGenes, direction = cm$Direction, camera_P = cm$PValue,
    posctrl_hits = sum(pc$logFC > 0 & pc$P < 0.05, na.rm = TRUE), posctrl_ok = pc_ok,
    verdict = fifelse(!pc_ok, "VOID_posctrl", fifelse(cm$PValue < 0.05, paste0("GRG_", cm$Direction), "GRG_ns"))),
    pc = cbind(dataset = label, pc), genes = data.table(dataset = label, gene = rownames(tt), logFC = tt$logFC, P = tt$P.Value)[gene %in% grg]) }
out <- list()
## 主检 GSE76161
s1 <- read_sm(file.path(DL, "GSE76161_series_matrix.txt.gz")); E1 <- annot(file.path(DL, "GPL11532.annot.gz"), s1$ex)
trt <- factor(ifelse(grepl("OCA", s1$getch("treatment")), "OCA", "DMSO"), levels = c("DMSO", "OCA")); don <- factor(s1$getch("patient_id"))
out$main <- one(E1, model.matrix(~ don + trt), "trtOCA", "GSE76161_PCLS_OCA_vs_DMSO")
## 复现 GSE124053（HepaRG）
s2 <- read_sm(file.path(DL, "GSE124053_series_matrix.txt.gz")); E2 <- annot(file.path(DL, "GPL6244.annot.gz"), s2$ex)
lig <- s2$getch("ligand treatment"); tp <- s2$getch("time point")
for (L in c("FXR-L", "CDCA")) for (T in c("4h", "24h")) { k <- (lig %in% c(L, "Control")) & tp == T
  g <- factor(ifelse(lig[k] == L, "trt", "ctrl"), levels = c("ctrl", "trt"))
  out[[paste(L, T)]] <- one(E2[, k], model.matrix(~ g), "gtrt", paste0("GSE124053_HepaRG_", L, "_", T)) }
S <- rbindlist(lapply(out, `[[`, "sum"))
write_tab(S, "18_fxr_summary.tsv"); write_tab(rbindlist(lapply(out, `[[`, "pc")), "18_fxr_posctrl.tsv")
write_tab(rbindlist(lapply(out, `[[`, "genes")), "18_fxr_GRG_genes.tsv"); print(S); print(fread(file.path(RES, "18_fxr_posctrl.tsv")))
