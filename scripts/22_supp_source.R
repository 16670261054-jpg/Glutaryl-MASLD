# 补充图源数据：S2 三套 bulk 中 4 个 hub 按纤维化分期；S5 病理分项散点所需表型（只取数，不做新检验）
if (nzchar(Sys.getenv("R_LIBS_GLUT"))) .libPaths(Sys.getenv("R_LIBS_GLUT")); suppressPackageStartupMessages(library(data.table))
PROJ <- Sys.getenv("GLUT_PROJ", "."); RES <- file.path(Sys.getenv("GLUT_OUT", PROJ), "results/22_supp"); dir.create(RES, showWarnings = FALSE)
STIFF <- Sys.getenv("GLUT_BULK", "data/bulk"); LAC <- Sys.getenv("GLUT_SN", "data/snRNA")
write_tab <- function(df, name) { p <- file.path(RES, name); tmp <- paste0(p, ".tmp"); fwrite(as.data.frame(df), tmp, sep = "\t"); file.rename(tmp, p); message("written: ", p) }
hub <- readLines(file.path(Sys.getenv("GLUT_OUT", PROJ), "results/10_gates/02_hub.txt"))
per <- fread(file.path(STIFF, "R1rep_scores_per_sample.tsv")); per[cohort == "GSE135251_NASH", cohort := "GSE135251"]
S2 <- rbindlist(lapply(c("GSE135251", "GSE162694", "GSE130970"), function(cid) { L <- readRDS(file.path(STIFF, sprintf("%s_log2cpm.rds", cid)))
  P <- per[cohort == cid & !is.na(fibrosis) & gsm %in% colnames(L)]; Z <- t(scale(t(L[hub, P$gsm])))
  melt(data.table(cohort = cid, gsm = P$gsm, fibrosis = P$fibrosis, t(Z)), id.vars = c("cohort", "gsm", "fibrosis"), variable.name = "gene", value.name = "z") }))
write_tab(S2, "S2_hub_by_stage.tsv")
## S5 表型
soft <- function(f) { ln <- readLines(f, warn = FALSE); gsm <- NA; rows <- list()
  for (x in ln) { if (grepl("^\\^SAMPLE = ", x)) { gsm <- sub("^\\^SAMPLE = ", "", x); rows[[gsm]] <- list(gsm = gsm) }
    else if (!is.na(gsm) && grepl("^!Sample_characteristics_ch1 = ", x)) { kv <- sub("^!Sample_characteristics_ch1 = ", "", x)
      rows[[gsm]][[trimws(sub(":.*", "", kv))]] <- trimws(sub("^[^:]*:", "", kv)) } }
  rbindlist(rows, fill = TRUE) }
p130 <- soft("D:/cc/MASH博士课题基因筛选/_tools/raw/cohort_scan/soft_GSE130970.txt")
p130 <- p130[, .(gsm, steatosis = as.numeric(`steatosis grade`), inflammation = as.numeric(`lobular inflammation grade`),
                 ballooning = as.numeric(`cytological ballooning grade`), fibrosis = as.numeric(`fibrosis stage`), NAS = as.numeric(`nafld activity score`))]
sm <- readLines(gzfile(file.path(LAC, "01_数据/GSE135251/GSE135251_series_matrix.txt.gz")))
g <- gsub('"', '', strsplit(grep("^!Sample_geo_accession", sm, value = TRUE), "\t")[[1]][-1])
gv <- function(k) as.numeric(trimws(sub(".*:", "", gsub('"', '', strsplit(grep(paste0("^!Sample_characteristics_ch1.*", k), sm, value = TRUE)[1], "\t")[[1]][-1]))))
p135 <- data.table(gsm = g, NAS = gv("nas score"), fibrosis = gv("fibrosis stage"))
sc <- fread(file.path(Sys.getenv("GLUT_OUT", PROJ), "results/14_batch1/C1_scores.tsv"))
S5 <- rbind(melt(merge(sc[cohort == "GSE130970"], p130, by = "gsm"), id.vars = c("gsm", "cohort", "score")),
            melt(merge(sc[cohort == "GSE135251"], p135, by = "gsm"), id.vars = c("gsm", "cohort", "score")))
write_tab(S5, "S5_pathology_long.tsv"); print(S5[, .N, by = .(cohort, variable)])
