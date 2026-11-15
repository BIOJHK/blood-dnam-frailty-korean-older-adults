#############################################################################
# 04. Genomic distribution of FI-associated CpGs (Reference / Hyper / Hypo)
# Manuscript section: Genomic distribution and enrichment analyses
#
# Input: EWAS_primary_FI_continuous.txt from 01_preprocessing_EWAS.R
#   - Significant CpGs : primary EWAS FDR < fdr_cut
#   - Hyper / Hypo     : FI coefficient (logFC) > 0 / < 0
#   - Reference        : all EPIC manifest probes (ref_set = "manifest")
#                        or all CpGs tested in the EWAS (ref_set = "tested")
# Features: relation to CpG island, DNase hypersensitivity, location relative to gene
# Test: Hyper vs Hypo, 2x2 per category (expected count < 5 -> Fisher, otherwise chi-square), BH correction
#############################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
})


## ===========================================================================
## 0. Settings — edit paths here only
## ===========================================================================
v2_dir <- "Results_directory"   # out_dir of 01_preprocessing_EWAS.R

cfg <- list(
  primary_file = file.path(v2_dir, "EWAS_primary_FI_continuous.txt"),
  out_dir      = file.path(v2_dir, "genomic_features"),

  fdr_cut = 0.1,
  ref_set = "manifest",   # "manifest": all EPIC probes / "tested": CpGs tested in the EWAS (= all.cpg in 05)

  island_order = c("Island", "N_Shelf", "N_Shore", "OpenSea", "S_Shelf", "S_Shore"),
  dnase_order  = c("DNase Sensitivity", "DNase Insensitivity"),
  group_order  = c("TSS1500", "TSS200", "5'UTR", "Exon", "Body", "3'UTR", "Intergenic"),

  colors = c(Reference = "#548235", Hyper = "#4472C4", Hypo = "#ED7D31"),
  plot_w = 7, plot_h = 4.5, dpi = 600
)
group_levels <- names(cfg$colors)   # Reference -> Hyper -> Hypo
dir.create(cfg$out_dir, showWarnings = FALSE, recursive = TRUE)
out <- function(f) file.path(cfg$out_dir, f)


## ===========================================================================
## 1. EWAS results -> significant CpGs + Hyper / Hypo
## ===========================================================================
ewas  <- read.delim(cfg$primary_file, check.names = FALSE)
valid <- !is.na(ewas$P.Value) & !is.na(ewas$adj.P.Val)
tested_cpg <- ewas$CpG[valid]

sig <- ewas[valid & ewas$adj.P.Val < cfg$fdr_cut, c("CpG", "logFC")]
sig$Direction <- ifelse(sig$logFC > 0, "Hyper", ifelse(sig$logFC < 0, "Hypo", NA))
sig <- sig[!is.na(sig$Direction), ]

message("FDR < ", cfg$fdr_cut, " CpGs: ", nrow(sig), "  <- check against the manuscript (930)",
        " (Hyper ", sum(sig$Direction == "Hyper"), " / Hypo ", sum(sig$Direction == "Hypo"), ")")


## ===========================================================================
## 2. Annotation (ilm10b4.hg19, same as 01) -> feature variables
## ===========================================================================
anno <- as.data.frame(getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19))
stopifnot(all(sig$CpG %in% rownames(anno)))

## Location relative to gene: one representative location (exclusive categories -> bar plot / 2x2 tests)
primary_group <- function(x) {
  x <- trimws(x); x[is.na(x)] <- ""
  first <- sub(";.*$", "", x)
  first[first == ""] <- "Intergenic"
  first[first %in% c("1stExon", "ExonBnd")] <- "Exon"
  first
}
dnase_status <- function(x) {
  ifelse(is.na(x) | trimws(x) == "", "DNase Insensitivity", "DNase Sensitivity")
}

make_feat <- function(cpgs, direction) {
  a <- anno[cpgs, ]
  data.frame(
    CpG                   = cpgs,
    Direction             = direction,
    Relation_to_Island    = a$Relation_to_Island,
    UCSC_RefGene_Group    = a$UCSC_RefGene_Group,
    DNase                 = dnase_status(a$DNase_Hypersensitivity_NAME),
    RefGene_Group_primary = primary_group(a$UCSC_RefGene_Group),
    stringsAsFactors      = FALSE
  )
}

ref_cpg   <- if (cfg$ref_set == "manifest") rownames(anno) else intersect(tested_cpg, rownames(anno))
ref_label <- if (cfg$ref_set == "manifest") "EPIC manifest (all probes)" else "CpGs tested in EWAS"
message("Reference: ", ref_label, ", ", length(ref_cpg), " probes (includes ",
        sum(ref_cpg %in% sig$CpG), " significant CpGs)")

ann_all <- rbind(make_feat(ref_cpg, "Reference"),
                 make_feat(sig$CpG, sig$Direction))
ann_all$Direction <- factor(ann_all$Direction, levels = group_levels)
print(table(ann_all$Direction))


## ===========================================================================
## 3. Category proportions (%) by group — for plotting
## ===========================================================================
## multi_label = TRUE counts every semicolon-separated location (use with raw UCSC_RefGene_Group)
get_membership <- function(x, multi_label) {
  if (!multi_label) return(as.character(x))
  x <- trimws(x); x[is.na(x)] <- ""
  lapply(strsplit(x, ";"), function(p) {
    if (length(p) == 0) return("Intergenic")
    p[p %in% c("1stExon", "ExonBnd")] <- "Exon"
    unique(p)
  })
}
in_cat <- function(membership, k) {
  if (is.list(membership)) vapply(membership, function(m) k %in% m, logical(1))
  else membership %in% k
}

summarize_category_pct <- function(df, var, level_order, multi_label = FALSE) {
  g   <- df$Direction
  mem <- get_membership(df[[var]], multi_label)
  n_g <- as.integer(table(g))
  res <- do.call(rbind, lapply(level_order, function(k) {
    data.frame(Category = k, Group = levels(g), n = n_g,
               n_in = as.integer(tapply(in_cat(mem, k), g, sum)))
  }))
  res$Percent  <- 100 * res$n_in / res$n
  res$Category <- factor(res$Category, levels = level_order)
  res$Group    <- factor(res$Group, levels = group_levels)
  res
}


## ===========================================================================
## 4. Hyper vs Hypo tests (2x2 per category, BH correction)
##    OR: Haldane-Anscombe correction (+0.5) if any cell is zero
## ===========================================================================
sig_stars <- function(p) {
  as.character(cut(p, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                   labels = c("***", "**", "*", "ns"), right = TRUE))
}

test_category_enrichment <- function(df, var, level_order,
                                     group1 = "Hyper", group2 = "Hypo", multi_label = FALSE) {
  sub <- df[df$Direction %in% c(group1, group2), ]
  gv  <- factor(as.character(sub$Direction), levels = c(group1, group2))
  mem <- get_membership(sub[[var]], multi_label)
  n1  <- sum(gv == group1); n2 <- sum(gv == group2)

  res <- do.call(rbind, lapply(level_order, function(k) {
    ic  <- factor(in_cat(mem, k), levels = c(TRUE, FALSE))
    tab <- table(gv, ic)
    x1  <- tab[group1, "TRUE"]; x2 <- tab[group2, "TRUE"]

    p <- NA_real_; test <- NA_character_; OR <- NA_real_
    if (x1 + x2 > 0 && x1 + x2 < n1 + n2) {            # not testable if one column is all zero
      exp_low <- any(suppressWarnings(chisq.test(tab)$expected) < 5)
      p    <- if (exp_low) fisher.test(tab)$p.value else chisq.test(tab, correct = TRUE)$p.value
      test <- if (exp_low) "Fisher" else "Chisq"
      cc   <- c(a = x1, b = tab[group1, "FALSE"], c = x2, d = tab[group2, "FALSE"])
      if (any(cc == 0)) cc <- cc + 0.5
      OR   <- unname((cc["a"] / cc["b"]) / (cc["c"] / cc["d"]))
    }
    data.frame(Category = k, n1 = n1, n2 = n2, x1 = x1, x2 = x2,
               pct1 = 100 * x1 / n1, pct2 = 100 * x2 / n2,
               OR = OR, p.value = p, test = test)
  }))

  res$p.adj    <- p.adjust(res$p.value, method = "BH")
  res$stars    <- sig_stars(res$p.adj)
  res$Category <- factor(res$Category, levels = level_order)
  names(res)[match(c("x1", "x2", "pct1", "pct2"), names(res))] <-
    c(paste0("n_", group1), paste0("n_", group2), paste0("pct_", group1), paste0("pct_", group2))
  rownames(res) <- NULL
  res
}


## ===========================================================================
## 5. Three-group dodged bar plot + Hyper/Hypo bracket & stars
## ===========================================================================
plot_grouped_pct <- function(pct_df, star_df, title, star_groups = c("Hyper", "Hypo"),
                             dodge_width = 0.7, bar_width = 0.6) {
  k      <- length(group_levels)
  offset <- function(gname) (match(gname, group_levels) - (k + 1) / 2) * (dodge_width / k)
  x1 <- offset(star_groups[1]); x2 <- offset(star_groups[2])
  y_max <- max(pct_df$Percent, na.rm = TRUE) * 1.25

  lab <- aggregate(Percent ~ Category, data = pct_df[pct_df$Group %in% star_groups, ], FUN = max)
  lab <- merge(lab, star_df[, c("Category", "stars")], by = "Category")
  lab$stars[is.na(lab$stars)] <- ""
  xi            <- as.numeric(lab$Category)
  lab$x_start   <- xi + x1
  lab$x_end     <- xi + x2
  lab$x_mid     <- xi + (x1 + x2) / 2
  lab$y_bracket <- lab$Percent + y_max * 0.035
  lab$y_tick    <- lab$y_bracket - y_max * 0.015
  lab$y_star    <- lab$y_bracket + y_max * 0.015

  ggplot(pct_df, aes(x = Category, y = Percent, fill = Group)) +
    geom_col(position = position_dodge(width = dodge_width), width = bar_width,
             color = "white", linewidth = 0.2) +
    geom_segment(data = lab, aes(x = x_start, xend = x_end, y = y_bracket, yend = y_bracket),
                 inherit.aes = FALSE, linewidth = 0.4) +
    geom_segment(data = lab, aes(x = x_start, xend = x_start, y = y_tick, yend = y_bracket),
                 inherit.aes = FALSE, linewidth = 0.4) +
    geom_segment(data = lab, aes(x = x_end, xend = x_end, y = y_tick, yend = y_bracket),
                 inherit.aes = FALSE, linewidth = 0.4) +
    geom_text(data = lab, aes(x = x_mid, y = y_star, label = stars),
              inherit.aes = FALSE, size = 4, fontface = "bold", vjust = 0) +
    scale_fill_manual(values = cfg$colors, breaks = group_levels) +
    scale_y_continuous(labels = function(x) paste0(x, "%"),
                       expand = expansion(mult = c(0, 0.08)), limits = c(0, y_max)) +
    labs(title = title,
         subtitle = paste0("Significance: ", star_groups[1], " vs ", star_groups[2],
                           " (BH-adjusted); Reference = ", ref_label),
         x = NULL, y = "Proportion of CpGs (%)", fill = "Group") +
    theme_classic(base_size = 13) +
    theme(panel.grid      = element_blank(),
          axis.line.x     = element_line(color = "black"),
          plot.subtitle   = element_text(size = 9, color = "grey35"),
          legend.position = "top")
}


## ===========================================================================
## 6. Run & save
## ===========================================================================
run_comparison <- function(var, level_order, title, label, multi_label = FALSE) {
  pct  <- summarize_category_pct(ann_all, var, level_order, multi_label)
  stat <- test_category_enrichment(ann_all, var, level_order, multi_label = multi_label)
  p    <- plot_grouped_pct(pct, stat, title)

  message("===== ", title, " : Hyper vs Hypo =====")
  print(stat[, c("Category", "n_Hyper", "n_Hypo", "pct_Hyper", "pct_Hypo",
                 "OR", "p.value", "p.adj", "stars", "test")])

  write.csv(pct,  out(paste0(label, "_percent.csv")),          row.names = FALSE)
  write.csv(stat, out(paste0(label, "_Hyper_vs_Hypo.csv")),    row.names = FALSE)
  ggsave(out(paste0(label, "_RefHyperHypo.tiff")), p,
         width = cfg$plot_w, height = cfg$plot_h, dpi = cfg$dpi)
  list(pct = pct, stat = stat, plot = p)
}

island_result  <- run_comparison("Relation_to_Island",    cfg$island_order,
                                 "Relation to CpG Island",        "island")
dnase_result   <- run_comparison("DNase",                 cfg$dnase_order,
                                 "DNase Hypersensitivity",        "dnase")
refgene_result <- run_comparison("RefGene_Group_primary", cfg$group_order,
                                 "CpG Location Relative to Gene", "refgene")

write.table(as.data.frame(table(Direction = ann_all$Direction)), out("group_counts.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
