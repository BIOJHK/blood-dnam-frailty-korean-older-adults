#############################################################################
# 03. Replication of previously reported CpGs: Mak et al. (2024) Aging Cell, doi:10.1111/acel.14135
# Manuscript section: Overlap with previously reported frailty-associated CpGs
#   - ST3: continuous FI (129 CpGs)          <-> primary EWAS in 01 (continuous FI)
#   - ST4: categorical, frail vs non-frail   <-> secondary EWAS in 01 (Frailtyfrail vs Non-frail)
#   Replication criterion: our P < rep_p  &  same direction of effect
#
# Inputs: outputs of 01_preprocessing_EWAS.R (out_dir)
#############################################################################

suppressPackageStartupMessages({
  library(readxl)
})


## ===========================================================================
## 0. Settings — edit paths here only
## ===========================================================================
v2_dir <- "Results_directory"   # out_dir of 01_preprocessing_EWAS.R

cfg <- list(
  primary_file   = file.path(v2_dir, "EWAS_primary_FI_continuous.txt"),
  secondary_file = file.path(v2_dir, "EWAS_secondary_Frailtyfrail_vs_Non-frail.txt"),
  rep_xlsx       = "replication_supplement_xlsx",   # Mak et al. 2024 Supplement 2 (ST3, ST4)
  out_dir        = file.path(v2_dir, "replication"),

  rep_p      = 0.05,                       # replication criterion (nominal P)
  fdr_cut    = 0.1,                        # our EWAS significance threshold (for flagging core CpGs)
  n_expected = c(ST3 = 129, ST4 = 538)     # numbers reported in the main text of the paper
)
dir.create(cfg$out_dir, showWarnings = FALSE, recursive = TRUE)
out <- function(f) file.path(cfg$out_dir, f)


## ===========================================================================
## 1. Read the published supplementary tables (ST3, ST4)
## ===========================================================================
col_names_mak <- c(
  "cpg", "gene", "chr", "position", "region", "relation_to_island",
  "beta_SATSA450K", "se_SATSA450K", "p_SATSA450K", "fdr_SATSA450K", "blank1",
  "beta_SATSAEPIC", "se_SATSAEPIC", "p_SATSAEPIC", "fdr_SATSAEPIC", "blank2",
  "beta_LSADT1997", "se_LSADT1997", "p_LSADT1997", "fdr_LSADT1997", "blank3",
  "beta_LSADT2007", "se_LSADT2007", "p_LSADT2007", "fdr_LSADT2007", "blank4",
  "beta_meta", "se_meta", "p_meta", "fdr_meta",
  "dir_SATSA450K", "dir_SATSAEPIC", "dir_LSADT1997", "dir_LSADT2007", "I2"
)

read_mak <- function(sheet) {
  x <- as.data.frame(read_excel(cfg$rep_xlsx, sheet = sheet, skip = 3, col_names = FALSE))[, 1:35]
  colnames(x) <- col_names_mak
  x <- x[!is.na(x$cpg), c("cpg", "gene", "beta_meta", "p_meta", "fdr_meta")]
  x$beta_meta     <- as.numeric(x$beta_meta)
  x$direction_rep <- sign(x$beta_meta)
  if (nrow(x) != cfg$n_expected[[sheet]])
    stop(sheet, ": ", nrow(x), " CpGs (does not match ", cfg$n_expected[[sheet]], " in the paper)")
  x
}

mak_s3 <- read_mak("ST3")   # continuous FI
mak_s4 <- read_mak("ST4")   # frail vs non-frail
message("ST3 & ST4 overlapping CpGs: ", length(intersect(mak_s3$cpg, mak_s4$cpg)))


## ===========================================================================
## 2. Read our EWAS results (01)
## ===========================================================================
read_ours <- function(file) {
  x <- read.delim(file, check.names = FALSE)
  rownames(x) <- x$CpG
  x
}
primary   <- read_ours(cfg$primary_file)
secondary <- read_ours(cfg$secondary_file)


## ===========================================================================
## 3. Comparison functions
## ===========================================================================
compare_rep <- function(rep, ours) {
  m <- merge(rep,
             data.frame(cpg = ours$CpG, logFC = ours$logFC,
                        P.Value = ours$P.Value, adj.P.Val = ours$adj.P.Val),
             by = "cpg")
  m$direction_ours  <- sign(m$logFC)
  m$direction_match <- m$direction_rep == m$direction_ours
  m$replicated      <- (m$P.Value < cfg$rep_p & m$direction_match) %in% TRUE
  m[order(m$P.Value), ]
}

## Expected proportion under the null = rep_p / 2 (P < rep_p and same direction by chance)
## Note: ignores correlation between CpGs, so the binomial test is approximate
summ_row <- function(label, n_reported, present, replicated, m = NULL) {
  k <- length(replicated); n <- length(present)
  data.frame(
    set        = label,
    reported   = n_reported,
    present    = n,
    nomP       = if (is.null(m)) NA else sum(m$P.Value < cfg$rep_p, na.rm = TRUE),
    same_dir   = if (is.null(m)) NA else sum(m$direction_match, na.rm = TRUE),
    replicated = k,
    pct        = round(100 * k / n, 1),
    binom_p    = binom.test(k, n, p = cfg$rep_p / 2, alternative = "greater")$p.value
  )
}


## ===========================================================================
## 4. Run comparisons
## ===========================================================================
merged_s3 <- compare_rep(mak_s3, primary)     # continuous FI
merged_s4 <- compare_rep(mak_s4, secondary)   # frail vs non-frail

rep_s3 <- merged_s3$cpg[merged_s3$replicated]
rep_s4 <- merged_s4$cpg[merged_s4$replicated]

## ST3 union ST4 (CpGs present in both are counted once)
present_union    <- union(merged_s3$cpg, merged_s4$cpg)
replicated_union <- union(rep_s3, rep_s4)

summary_tab <- rbind(
  summ_row("ST3_continuous",  nrow(mak_s3), merged_s3$cpg, rep_s3, merged_s3),
  summ_row("ST4_categorical", nrow(mak_s4), merged_s4$cpg, rep_s4, merged_s4),
  summ_row("union",           length(union(mak_s3$cpg, mak_s4$cpg)),
           present_union, replicated_union)
)
print(summary_tab)


## ===========================================================================
## 5. Core CpGs replicated in both ST3 and ST4
## ===========================================================================
core <- intersect(rep_s3, rep_s4)
core_tab <- data.frame(
  cpg              = core,
  gene             = mak_s3$gene[match(core, mak_s3$cpg)],
  in_primary_FDR   = core %in% primary$CpG[primary$adj.P.Val < cfg$fdr_cut],
  in_secondary_FDR = core %in% secondary$CpG[secondary$adj.P.Val < cfg$fdr_cut]
)
message("CpGs replicated in both ST3 and ST4: ", length(core))
print(core_tab)


## ===========================================================================
## 6. Save
## ===========================================================================
write.csv(merged_s3, out("concordance_ST3_continuous.csv"),  row.names = FALSE)
write.csv(merged_s4, out("concordance_ST4_categorical.csv"), row.names = FALSE)
write.csv(core_tab,  out("core_replicated_CpGs.csv"),        row.names = FALSE)
write.table(summary_tab, out("replication_summary.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
