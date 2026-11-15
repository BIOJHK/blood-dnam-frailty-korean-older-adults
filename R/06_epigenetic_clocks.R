#############################################################################
# 06. Epigenetic clocks
# Manuscript section: Estimation of epigenetic age
#
# Methods: Horvath, Hannum, Skin & Blood, PhenoAge, GrimAge, GrimAge v2
#   - Main clocks : methylclock + meffonym
#       methylclock : Horvath, Hannum, Skin & Blood (skinHorvath), PhenoAge (Levine)
#       meffonym    : GrimAge, GrimAge v2 (+ Horvath/Hannum/Skin/PhenoAge as a cross-check)
#   - PC clocks   : methylCIPHER::calcPCClocks (Higgins-Chen et al. 2022)
#   - DunedinPACE (v0.99.0): pace of aging, original values used (no age residualization)
#   - EAA         : residual from regressing each clock's predicted age on chronological age
#   - Association of FI with each aging measure: Pearson correlation
#
# Inputs: outputs of 01_preprocessing_EWAS.R (out_dir)
#   - beta_ssNoob_BMIQ_filtered.rds : ssNoob -> BMIQ beta (CpG x sample)
#   - pheno_with_cellprop.rds       : FI, Age, Sex (F/M)
#
# Package installation (once; not part of this script):
#   BiocManager::install("methylclock")
#   devtools::install_github("perishky/meffonym")
#   devtools::install_github("MorganLevineLab/methylCIPHER")
#   devtools::install_github("danbelsky/DunedinPACE")
#   The PC clock coefficient file (CalcAllPCClocks.RData) is downloaded separately into pcclock_dir
#############################################################################

suppressPackageStartupMessages({
  library(methylclock)
  library(meffonym)
  library(methylCIPHER)
  library(DunedinPACE)
})


## ===========================================================================
## 0. Settings — edit paths here only
## ===========================================================================
v2_dir <- "Results_directory"   # out_dir of 01_preprocessing_EWAS.R

cfg <- list(
  beta_file   = file.path(v2_dir, "beta_ssNoob_BMIQ_filtered.rds"),
  pheno_file  = file.path(v2_dir, "pheno_with_cellprop.rds"),
  pcclock_dir = "pcclock_directory",            # folder containing CalcAllPCClocks.RData
  out_dir     = file.path(v2_dir, "epigenetic_clocks"),

  ## methylclock clock name -> output name
  mc_clocks = c(Horvath = "Horvath", Hannum = "Hannum",
                skinHorvath = "SkinBlood", Levine = "PhenoAge"),

  ## meffonym model name -> output name (calibrate: kept from the original script)
  mf_models = c(horvath = "Horvath", hannum = "Hannum", skin = "SkinBlood",
                phenoage = "PhenoAge", grimage = "GrimAge", grimagev2 = "GrimAge2"),
  mf_calibrate = c("horvath"),
  mf_need_covars = c("grimage", "grimagev2"),   # models that need age / female rows

  ## DunedinPACE: proportion of required probes (kept from the original script)
  dp_prop = 0.7
)
dir.create(cfg$out_dir, showWarnings = FALSE, recursive = TRUE)
out <- function(f) file.path(cfg$out_dir, f)


## ===========================================================================
## 1. Read & align data
## ===========================================================================
beta  <- readRDS(cfg$beta_file)        # CpG x sample
pheno <- readRDS(cfg$pheno_file)
stopifnot(identical(rownames(pheno), colnames(beta)))

FI     <- pheno$FI
age    <- pheno$Age
female <- as.integer(pheno$Sex == "F")  # 01: Sex factor (F/M) -> female = 1
stopifnot(!anyNA(age), !anyNA(female))
message("Samples: ", ncol(beta), " (female ", sum(female), ", male ", sum(female == 0), ")")

clocks <- data.frame(Sentrix_id  = colnames(beta),
                     Sample_Name = pheno$Sample_Name,
                     FI = FI, Age = age, Female = female)


## ===========================================================================
## 2. methylclock: Horvath, Hannum, Skin & Blood, PhenoAge
## ===========================================================================
invisible(checkClocks(beta))   # report missing CpGs per clock (CpGs removed by QC filtering)

mc <- DNAmAge(beta, clocks = names(cfg$mc_clocks), normalize = FALSE, cell.count = FALSE)
mc <- as.data.frame(mc)
stopifnot(identical(as.character(mc$id), colnames(beta)))
for (k in names(cfg$mc_clocks)) clocks[[paste0("mc_", cfg$mc_clocks[[k]])]] <- mc[[k]]


## ===========================================================================
## 3. meffonym: GrimAge, GrimAge v2 (+ clocks for cross-checking)
## ===========================================================================
missing_models <- setdiff(names(cfg$mf_models), meffonym.models())
if (length(missing_models)) stop("Models not available in meffonym: ", paste(missing_models, collapse = ", "))

beta_cov <- rbind(beta, female = female, age = age)   # input for GrimAge-type models

mf_coverage <- data.frame()
for (m in names(cfg$mf_models)) {
  x   <- if (m %in% cfg$mf_need_covars) beta_cov else beta
  res <- meffonym.score(x, m, calibrate = m %in% cfg$mf_calibrate)
  clocks[[paste0("mf_", cfg$mf_models[[m]])]] <- res$score

  sites <- setdiff(names(meffonym.get.model(m)$coefficients), c("(Intercept)", "intercept", "age", "female"))
  mf_coverage <- rbind(mf_coverage, data.frame(model = m, n_CpG = length(sites),
                                               n_present = sum(sites %in% rownames(beta))))
}
print(mf_coverage)


## ===========================================================================
## 4. methylCIPHER: PC clocks
## ===========================================================================
## methylCIPHER input: sample x CpG; datPheno needs Age / Female
## Missing CpGs are imputed internally with training means
datPheno <- data.frame(Age = age, Female = female, row.names = colnames(beta))
pc <- calcPCClocks(path_to_PCClocks_directory = cfg$pcclock_dir,
                   datMeth  = t(beta),
                   datPheno = datPheno)
pc_cols <- grep("^PC", colnames(pc), value = TRUE)
clocks[pc_cols] <- pc[colnames(beta), pc_cols]


## ===========================================================================
## 5. DunedinPACE (pace of aging; 1.0 = one year of biological aging per chronological year)
## ===========================================================================
dp_req <- getRequiredProbes(backgroundList = FALSE)$DunedinPACE   # 173 model CpGs
dp_bg  <- getRequiredProbes(backgroundList = TRUE)$DunedinPACE    # 20,000 incl. normalization background
dp_coverage <- data.frame(set = c("model_173", "with_background_20000"),
                          n_CpG = c(length(dp_req), length(dp_bg)),
                          n_present = c(sum(dp_req %in% rownames(beta)), sum(dp_bg %in% rownames(beta))))
print(dp_coverage)

dp <- PACEProjector(beta, proportionOfProbesRequired = cfg$dp_prop)$DunedinPACE
if (!is.null(names(dp))) dp <- dp[colnames(beta)]   # align sample order if named
stopifnot(length(dp) == ncol(beta))
clocks$DunedinPACE <- as.numeric(dp)                # original (non-residualized) values


## ===========================================================================
## 6. EAA (residual of clock ~ chronological age)
##    DunedinPACE is not residualized (original values used)
## ===========================================================================
## Only age-scale PC clocks (PCDNAmTL and GrimAge components such as PCPACKYRS are excluded)
pc_age <- intersect(c("PCHorvath1", "PCHorvath2", "PCHannum", "PCPhenoAge", "PCGrimAge"), pc_cols)
age_clocks <- c(grep("^mc_|^mf_", colnames(clocks), value = TRUE), pc_age)
for (k in age_clocks) {
  ok <- !is.na(clocks[[k]])
  clocks[[paste0(k, "_Accel")]] <- NA_real_
  clocks[[paste0(k, "_Accel")]][ok] <- resid(lm(clocks[[k]][ok] ~ age[ok]))
}


## ===========================================================================
## 7. QC: correlation with age & agreement between packages (methylclock vs meffonym)
## ===========================================================================
age_cor <- data.frame(clock = age_clocks,
                      r_age = sapply(age_clocks, function(k) cor(clocks[[k]], age, use = "complete.obs")),
                      mean  = sapply(age_clocks, function(k) mean(clocks[[k]], na.rm = TRUE)),
                      row.names = NULL)
print(age_cor)

shared <- intersect(cfg$mc_clocks, cfg$mf_models)
pkg_agree <- do.call(rbind, lapply(shared, function(k) {
  a <- clocks[[paste0("mc_", k)]]; b <- clocks[[paste0("mf_", k)]]
  data.frame(clock = k, r = cor(a, b, use = "complete.obs"),
             mean_diff_mc_minus_mf = mean(a - b, na.rm = TRUE))
}))
print(pkg_agree)


## ===========================================================================
## 8. Pearson correlation of FI with epigenetic aging measures
##    EAA (each clock) + DunedinPACE (original values); predicted ages also recorded for reference
## ===========================================================================
fi_measures <- data.frame(
  measure = c(paste0(age_clocks, "_Accel"), "DunedinPACE", age_clocks),
  type    = c(rep("EAA", length(age_clocks)), "DunedinPACE", rep("predicted_age", length(age_clocks))),
  stringsAsFactors = FALSE
)
fi_cor <- do.call(rbind, lapply(seq_len(nrow(fi_measures)), function(i) {
  k  <- fi_measures$measure[i]
  ct <- cor.test(clocks$FI, clocks[[k]], method = "pearson")
  data.frame(measure = k, type = fi_measures$type[i],
             n = sum(complete.cases(clocks$FI, clocks[[k]])),
             r = unname(ct$estimate),
             CI_low = ct$conf.int[1], CI_high = ct$conf.int[2],
             p.value = ct$p.value)
}))
print(fi_cor)


## ===========================================================================
## 9. Save
## ===========================================================================
write.table(clocks,      out("epigenetic_clocks.txt"),          sep = "\t", quote = FALSE, row.names = FALSE)
write.table(age_cor,     out("QC_clock_age_correlation.txt"),   sep = "\t", quote = FALSE, row.names = FALSE)
write.table(pkg_agree,   out("QC_methylclock_vs_meffonym.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(mf_coverage, out("QC_meffonym_CpG_coverage.txt"),   sep = "\t", quote = FALSE, row.names = FALSE)
write.table(dp_coverage, out("QC_DunedinPACE_CpG_coverage.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(fi_cor,      out("FI_vs_epigenetic_aging_pearson.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(clocks, out("epigenetic_clocks.rds"))

## Check versions against the Methods
## (R 4.6.1, methylclock 1.18.0, meffonym 1.2.3, methylCIPHER 0.2.0, DunedinPACE 0.99.0)
versions <- c(R = paste(R.version$major, R.version$minor, sep = "."),
              sapply(c("methylclock", "meffonym", "methylCIPHER", "DunedinPACE"),
                     function(p) as.character(packageVersion(p))))
print(versions)
writeLines(paste(names(versions), versions, sep = "\t"), out("package_versions.txt"))
