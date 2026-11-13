# FXR-agonist datasets: GRG change versus mitochondria-matched random gene sets
if (nzchar(Sys.getenv("R_LIBS_GLUT"))) .libPaths(Sys.getenv("R_LIBS_GLUT"))
suppressPackageStartupMessages({ library(data.table); library(limma); library(org.Hs.eg.db); library(AnnotationDbi) })
PROJ <- Sys.getenv("GLUT_PROJ", "."); RES <- file.path(Sys.getenv("GLUT_OUT", PROJ), "results/26_editor_checks"); dir.create(RES, showWarnings = FALSE)
LAC <- Sys.getenv("GLUT_SN", "data/snRNA"); N <- 1000
write_tab <- function(df, name) { p <- file.path(RES, name); tmp <- paste0(p, ".tmp"); fwrite(as.data.frame(df), tmp, sep = "\t"); file.rename(tmp, p); message("written: ", p) }
grg <- readLines(file.path(PROJ, "geneset/GRG_full.txt")); KEY <- c("GCDH", "AASS", "ACADSB", "SARDH", "ACOT13", "ACAA2")
mito <- unique(na.omit(suppressMessages(AnnotationDbi::select(org.Hs.eg.db, "GO:0005739", "SYMBOL", "GOALL"))$SYMBOL))

mk_draw <- function(mu, pool, nb, g) { br <- unique(quantile(mu[pool], 0:nb / nb)); b <- findInterval(mu[pool], br, rightmost.closed = TRUE, all.inside = TRUE); names(b) <- pool
  hb <- findInterval(mu[g], br, rightmost.closed = TRUE, all.inside = TRUE); tb <- table(hb)
  function() unlist(lapply(names(tb), function(x) { cand <- names(b)[b == as.integer(x)]; sample(cand, tb[[x]], replace = length(cand) < tb[[x]]) })) }
## ---------- FXR mito null ----------
FD <- file.path(PROJ, "data/fxr")
read_sm <- function(f) { sm <- readLines(gzfile(f)); i0 <- grep("^!series_matrix_table_begin", sm); i1 <- grep("^!series_matrix_table_end", sm)
  ex <- fread(text = sm[(i0 + 1):(i1 - 1)]); id <- as.character(ex[[1]]); ex <- as.matrix(ex[, -1]); rownames(ex) <- id
  ch <- sm[grepl("^!Sample_characteristics_ch1", sm)]
  getch <- function(k) { ln <- ch[grepl(k, ch)][1]; trimws(sub("^[^:]*:", "", gsub('"', '', strsplit(ln, "\t")[[1]][-1]))) }
  list(ex = ex, getch = getch) }
annot <- function(f, ex) { ln <- readLines(gzfile(f)); ln <- ln[!grepl("^[!#^]", ln)]; a <- fread(text = ln, sep = "\t", quote = "", select = c("ID", "Gene symbol"))
  setnames(a, c("ID", "sym")); a <- a[sym != "" & !grepl("///", sym)][as.character(ID) %in% rownames(ex)]
  e <- ex[as.character(a$ID), ]; o <- order(-rowMeans(e)); e <- e[o, ]; s <- a$sym[o]; k <- !duplicated(s); e <- e[k, ]; rownames(e) <- s[k]; e }
nullN3 <- function(E, des, coef, label) {
  fit <- eBayes(lmFit(E, des)); tv <- setNames(topTable(fit, coef = coef, n = Inf, sort.by = "none")$t, rownames(E))
  g <- rownames(E)[rownames(E) %in% grg]; T0 <- mean(tv[g]); m <- rowMeans(E); pool <- setdiff(intersect(rownames(E), mito), grg)
  set.seed(20261122); dr <- mk_draw(m, pool, 5, g); Tn <- replicate(N, mean(tv[dr()]))
  data.table(dataset = label, n_GRG = length(g), mean_t_GRG = T0, null_median_N3 = median(Tn), null_P_N3 = (sum(abs(Tn) >= abs(T0)) + 1) / (N + 1)) }
out <- list()
s1 <- read_sm(file.path(FD, "GSE76161_series_matrix.txt.gz")); E1 <- annot(file.path(FD, "GPL11532.annot.gz"), s1$ex)
trt <- factor(ifelse(grepl("OCA", s1$getch("treatment")), "OCA", "DMSO"), levels = c("DMSO", "OCA")); don <- factor(s1$getch("patient_id"))
out$pcls <- nullN3(E1, model.matrix(~ don + trt), "trtOCA", "GSE76161_PCLS_OCA_vs_DMSO")
s2 <- read_sm(file.path(FD, "GSE124053_series_matrix.txt.gz")); E2 <- annot(file.path(FD, "GPL6244.annot.gz"), s2$ex)
lig <- s2$getch("ligand treatment"); tp <- s2$getch("time point")
for (Lg in c("FXR-L", "CDCA")) for (Tm in c("4h", "24h")) { k <- (lig %in% c(Lg, "Control")) & tp == Tm
  gg <- factor(ifelse(lig[k] == Lg, "trt", "ctrl"), levels = c("ctrl", "trt"))
  out[[paste(Lg, Tm)]] <- nullN3(E2[, k], model.matrix(~ gg), "ggtrt", paste0("GSE124053_HepaRG_", Lg, "_", Tm)) }
S2 <- rbindlist(out); write_tab(S2, "E2_fxr_mito_null.tsv"); print(S2)
