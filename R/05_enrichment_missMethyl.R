#############################################################################
# 05. Gene set enrichment of FI-associated CpGs (missMethyl::gsameth)
# Manuscript section: Genomic distribution and enrichment analyses
#
# Input: EWAS_primary_FI_continuous.txt from 01_preprocessing_EWAS.R
#   - sig.cpg : primary EWAS FDR < fdr_cut
#   - all.cpg : all CpGs tested in the EWAS (background)
# Settings: CpGs in promoter regions (TSS1500 / TSS200 / 5'UTR / 1stExon),
#           prior.prob = TRUE (probe-number bias correction),
#           results recorded for both fract.counts = TRUE and FALSE
# Libraries: Reactome 2024, GO Biological Process 2026 (GMT)
#############################################################################

suppressPackageStartupMessages({
  library(missMethyl)
  library(org.Hs.eg.db)
  library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  library(ggplot2)
})


## ===========================================================================
## 0. Settings — edit paths here only
## ===========================================================================
v2_dir <- "Results_directory"   # out_dir of 01_preprocessing_EWAS.R

cfg <- list(
  primary_file = file.path(v2_dir, "EWAS_primary_FI_continuous.txt"),
  reactome_gmt = "reactome_gmt_file",   # Reactome Pathways 2024 (GMT)
  gobp_gmt     = "gobp_gmt_file",       # GO Biological Process 2026 (GMT)
  old_file     = NULL,                  # (optional) set a path only to compare with earlier EWAS results, e.g. "old_660k_file"
  out_dir      = file.path(v2_dir, "enrichment"),

  fdr_cut      = 0.1,
  promoter_grp = c("TSS1500", "TSS200", "5'UTR", "1stExon"),
  top_k        = 20,   # number of top pathways whose gene tables are saved
  plot_n       = 10    # number of top pathways drawn in the dot plot
)
dir.create(cfg$out_dir, showWarnings = FALSE, recursive = TRUE)
out <- function(f) file.path(cfg$out_dir, f)


## ===========================================================================
## 1. Read EWAS results -> sig.cpg / all.cpg
## ===========================================================================
ewas <- read.delim(cfg$primary_file, check.names = FALSE)
rownames(ewas) <- ewas$CpG

valid   <- !is.na(ewas$P.Value) & !is.na(ewas$adj.P.Val)
all.cpg <- rownames(ewas)[valid]
sig.cpg <- rownames(ewas)[valid & ewas$adj.P.Val < cfg$fdr_cut]
message("Background CpGs: ", length(all.cpg),
        " | FDR < ", cfg$fdr_cut, " CpGs: ", length(sig.cpg), "  <- check against the manuscript (930)")

## (Optional) overlap with earlier EWAS results
if (!is.null(cfg$old_file)) {
  old     <- read.table(cfg$old_file)
  old_sig <- rownames(old)[old$adj.P.Val < cfg$fdr_cut]
  print(c(new              = length(sig.cpg),
          old              = length(old_sig),
          overlap          = length(intersect(sig.cpg, old_sig)),
          universe_overlap = length(intersect(all.cpg, rownames(old)))))
}

## Number of significant / background genes (promoter-based; for Methods and Results)
mapped <- getMappedEntrezIDs(sig.cpg = sig.cpg, all.cpg = all.cpg, array.type = "EPIC",
                             genomic.features = cfg$promoter_grp)
map_n <- c(sig_CpG        = length(sig.cpg),
           all_CpG        = length(all.cpg),
           sig_genes      = length(mapped$sig.eg),
           universe_genes = length(mapped$universe))
print(map_n)
write.table(data.frame(item = names(map_n), n = map_n), out("mapped_gene_counts_promoter.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)


## ===========================================================================
## 2. Helper functions
## ===========================================================================
ann <- as.data.frame(getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19))

read_gmt <- function(file) {                    # read GMT (strip whitespace / hidden characters)
  ln <- strsplit(readLines(file, warn = FALSE), "\t")
  sets <- lapply(ln, function(x) {
    g <- trimws(x[-(1:2)]); g <- gsub(",.*", "", g)
    unique(g[g != ""])
  })
  names(sets) <- trimws(sapply(ln, `[`, 1))
  sets
}

sym_to_eg <- function(s) {                      # symbol -> Entrez (SYMBOL, then ALIAS if that fails)
  s <- unique(s)
  m <- suppressMessages(mapIds(org.Hs.eg.db, keys = s, column = "ENTREZID",
                               keytype = "SYMBOL", multiVals = "first"))
  miss <- names(m)[is.na(m)]
  if (length(miss) > 0) {
    a <- tryCatch(suppressMessages(mapIds(org.Hs.eg.db, keys = miss, column = "ENTREZID",
                                          keytype = "ALIAS", multiVals = "first")),
                  error = function(e) NULL)
    if (!is.null(a)) m[miss] <- a[miss]
  }
  m
}

make_pairs <- function(cpgs) {                  # CpG - gene - region pairs
  a <- ann[intersect(cpgs, rownames(ann)), c("UCSC_RefGene_Name", "UCSC_RefGene_Group")]
  g <- strsplit(a$UCSC_RefGene_Name, ";")
  f <- strsplit(a$UCSC_RefGene_Group, ";")
  ok <- lengths(g) > 0 & lengths(g) == lengths(f)
  data.frame(cpg   = rep(rownames(a)[ok], lengths(g[ok])),
             gene  = unlist(g[ok]),
             group = unlist(f[ok]))
}


## ===========================================================================
## 3. Promoter probe counts per gene
## ===========================================================================
all_prom <- make_pairs(all.cpg)
all_prom <- unique(all_prom[all_prom$group %in% cfg$promoter_grp, c("cpg", "gene")])
sig_prom <- all_prom[all_prom$cpg %in% sig.cpg, ]

n_probe <- table(all_prom$gene)
n_sig   <- table(sig_prom$gene)
sig_prom$n_genes_on_cpg <- as.integer(table(sig_prom$cpg)[sig_prom$cpg])

gene_tab <- data.frame(gene = names(n_probe), n_probe = as.integer(n_probe))
gene_tab$n_sig_probe <- as.integer(n_sig[gene_tab$gene])
gene_tab$n_sig_probe[is.na(gene_tab$n_sig_probe)] <- 0
gene_tab$max_genes_on_sig_cpg <- as.integer(tapply(sig_prom$n_genes_on_cpg,
                                                   sig_prom$gene, max)[gene_tab$gene])
gene_tab$sig_cpgs <- tapply(sig_prom$cpg, sig_prom$gene, paste, collapse = ";")[gene_tab$gene]
gene_tab$entrez   <- sym_to_eg(gene_tab$gene)[gene_tab$gene]

gene_tab <- gene_tab[order(-gene_tab$n_sig_probe, -gene_tab$n_probe), ]
write.csv(gene_tab, out("gene_probe_counts_promoter.csv"), row.names = FALSE)


## ===========================================================================
## 4. gsameth per library (fract.counts TRUE / FALSE)
## ===========================================================================
run_library <- function(gmt_file, label, min_size = 1, max_size = Inf) {
  sets_sym <- read_gmt(gmt_file)
  m <- sym_to_eg(unlist(sets_sym))
  message(label, " symbol conversion rate: ", round(mean(!is.na(m)) * 100, 1), "%")
  sets_eg <- lapply(sets_sym, function(s) unique(na.omit(m[s])))
  sets_eg <- sets_eg[lengths(sets_eg) >= min_size & lengths(sets_eg) <= max_size]

  args <- list(sig.cpg = sig.cpg, all.cpg = all.cpg, collection = sets_eg,
               array.type = "EPIC", genomic.features = cfg$promoter_grp,
               prior.prob = TRUE, sig.genes = TRUE, plot.bias = FALSE)
  r_frac   <- topGSA(do.call(gsameth, c(args, fract.counts = TRUE)),  number = Inf)
  r_nofrac <- topGSA(do.call(gsameth, c(args, fract.counts = FALSE)), number = Inf)

  terms <- rownames(r_frac)                     # ordered by P_frac
  summ <- data.frame(
    Term       = terms,
    N          = r_frac$N,
    DE_frac    = r_frac$DE,
    P_frac     = r_frac$P.DE,
    FDR_frac   = r_frac$FDR,
    DE_nofrac  = r_nofrac[terms, "DE"],
    P_nofrac   = r_nofrac[terms, "P.DE"],
    FDR_nofrac = r_nofrac[terms, "FDR"]
  )
  summ$n_sig_cpg <- sapply(terms, function(t) {
    cp <- gene_tab$sig_cpgs[gene_tab$entrez %in% sets_eg[[t]] & gene_tab$n_sig_probe > 0]
    length(unique(unlist(strsplit(na.omit(cp), ";"))))
  })
  summ$genes_per_sig_cpg <- round(summ$DE_nofrac / pmax(summ$n_sig_cpg, 1), 2)
  summ$median_probe <- sapply(terms, function(t)
    median(gene_tab$n_probe[gene_tab$entrez %in% sets_eg[[t]]]))
  summ$SigGenes <- r_frac$SigGenesInSet

  write.csv(summ, out(paste0(label, "_pathway_summary.csv")), row.names = FALSE)
  list(summary = summ, sets_eg = sets_eg)
}

## Per-gene probe table for one pathway
pathway_genes <- function(res, term) {
  x <- gene_tab[gene_tab$entrez %in% res$sets_eg[[term]],
                c("gene", "n_probe", "n_sig_probe", "max_genes_on_sig_cpg", "sig_cpgs")]
  x[order(-x$n_sig_probe, -x$n_probe), ]
}

## Save the significant-gene tables of the top k pathways in one file
save_top_genes <- function(res, label, k = cfg$top_k) {
  tt  <- head(res$summary$Term, k)
  tab <- do.call(rbind, lapply(tt, function(t) {
    x <- pathway_genes(res, t); x <- x[x$n_sig_probe > 0, ]
    if (nrow(x) == 0) return(NULL)
    cbind(Term = t, x)
  }))
  write.csv(tab, out(paste0(label, "_top", k, "_pathway_genes.csv")), row.names = FALSE)
  tab
}


## ===========================================================================
## 5. Run: Reactome 2024 / GO BP 2026
## ===========================================================================
reac <- run_library(cfg$reactome_gmt, "Reactome2024")
print(head(reac$summary[, 1:11], 10))
save_top_genes(reac, "Reactome2024")
saveRDS(reac, out("Reactome2024_gsameth.rds"))

gobp <- run_library(cfg$gobp_gmt, "GOBP2026")
print(head(gobp$summary[, 1:11], 20))
save_top_genes(gobp, "GOBP2026")
saveRDS(gobp, out("GOBP2026_gsameth.rds"))

## Example — genes of a specific pathway: pathway_genes(reac, reac$summary$Term[1])


## ===========================================================================
## 6. Dot plot (top plot_n pathways, x = -log10 P_frac)
## ===========================================================================
plot_dot <- function(dat, size_var, size_title, title, xlim = NULL) {
  ggplot(dat, aes(x = logP, y = reorder(Term, logP))) +
    geom_point(aes(colour = logP, size = .data[[size_var]])) +
    scale_colour_gradient(low = "#3E6DB5", high = "#E0281E", guide = "none") +
    scale_size(range = c(3, 8), name = size_title) +
    scale_x_continuous(limits = xlim, labels = function(x) sprintf("%.1f", x)) +
    labs(title = title,
         x = expression(-log[10](italic("p-value"))),
         y = NULL) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.y         = element_text(colour = "black"),
      axis.text.x         = element_text(colour = "black", size = 12),
      axis.title.x        = element_text(size = 18),
      plot.title          = element_text(size = 16),
      plot.title.position = "plot"
    )
}

## xlim = NULL gives an automatic range. To fix it, e.g. xlim = c(1.8, 3.6)
## (points outside the range are dropped from the plot, so check against new results)
save_dotplots <- function(res, label, title, xlim = NULL) {
  df <- head(res$summary, cfg$plot_n)
  df$Term <- sub(" \\(GO:.*\\)$", "", df$Term)   # strip GO IDs
  df$logP <- -log10(df$P_frac)

  p_de <- plot_dot(df, "DE_frac", "Gene count",    title, xlim)   # size = significant genes in term
  p_n  <- plot_dot(df, "N",       "Gene set size", title, xlim)   # size = total genes in term
  ggsave(out(paste0(label, "_dotplot_DEfrac.tiff")), p_de, width = 8, height = 3.8, dpi = 600)
  ggsave(out(paste0(label, "_dotplot_N.tiff")),      p_n,  width = 8, height = 3.8, dpi = 600)
  invisible(list(DE = p_de, N = p_n))
}

save_dotplots(reac, "Reactome2024", "Reactome Pathway")
save_dotplots(gobp, "GOBP2026",     "Gene Ontology: Biological Process")
