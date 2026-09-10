setwd("/path/to/your/working/directory")

carriers_file <- "/path/to/your/data_file.tsv"
meta_file     <- "/path/to/your/metadata.xlsx"
out_dir       <- "RDA_ARG_carriers_genus"

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

N_PERM <- 999

suppressPackageStartupMessages({
  library(readr)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(vegan)
  library(ggplot2)
})

fmt_sci <- function(x) ifelse(is.na(x), NA_character_, sprintf("%.2E", as.numeric(x)))


cat("[1] Reading ARG carrier table...\n")
if (!file.exists(carriers_file)) stop("File not found: ", carriers_file)

carriers <- readr::read_tsv(carriers_file, show_col_types = FALSE)

need_cols <- c("ARG_GENE","RESISTANCE","sample","sequence",
               "TPM_prop","Phylum","Class","Order","Family","Genus","Species")
miss <- setdiff(need_cols, names(carriers))
if (length(miss) > 0) {
  warning("Carrier table missing columns: ", paste(miss, collapse = ", "),
          ". Some later steps might not use these columns.")
}

carriers <- carriers %>%
  mutate(
    sample   = as.character(sample),
    Genus    = as.character(Genus),
    TPM_prop = suppressWarnings(as.numeric(TPM_prop))
  ) %>%
  filter(!is.na(TPM_prop))

cat("  Record count: ", nrow(carriers), "\n")


cat("[2] Reading metadata and aligning samples...\n")
if (!file.exists(meta_file)) stop("Metadata not found: ", meta_file)

meta_raw <- read_excel(meta_file) %>%
  mutate(across(everything(), ~ if (is.character(.)) trimws(.) else .))

need_meta <- c("Sample_id","run_accession","category","continent","origin")
miss_meta <- setdiff(need_meta, names(meta_raw))
if (length(miss_meta) > 0) {
  stop("Metadata missing columns: ", paste(miss_meta, collapse = ", "))
}

meta_raw <- meta_raw %>%
  mutate(
    Sample_id     = as.character(Sample_id),
    run_accession = as.character(run_accession),
    category      = as.character(category),
    continent     = as.character(continent),
    origin        = as.character(origin)
  )

carriers2 <- carriers %>%
  left_join(meta_raw %>% select(Sample_id, run_accession, category, continent, origin),
            by = c("sample" = "run_accession")) %>%
  filter(!is.na(Sample_id))

cat("  Records matching Sample_id: ", nrow(carriers2), "\n")


cat("[3] Building genus-level matrix (summing TPM_prop)...\n")

genus_tab <- carriers2 %>%
  filter(!is.na(Genus) & Genus != "") %>%
  group_by(Sample_id, Genus) %>%
  summarise(TPM = sum(TPM_prop, na.rm = TRUE), .groups = "drop")

cat("  Unique samples: ", n_distinct(genus_tab$Sample_id), "\n")
cat("  Unique genera: ",   n_distinct(genus_tab$Genus), "\n")

genus_wide <- genus_tab %>%
  tidyr::pivot_wider(
    id_cols     = Sample_id,
    names_from  = Genus,
    values_from = TPM,
    values_fill = 0
  )

genus_mat <- genus_wide %>%
  column_to_rownames("Sample_id") %>%
  as.matrix()

row_sum <- rowSums(genus_mat, na.rm = TRUE)
col_sum <- colSums(genus_mat, na.rm = TRUE)

keep_row <- row_sum > 0
keep_col <- col_sum > 0

genus_mat <- genus_mat[keep_row, keep_col, drop = FALSE]

cat("  Filtered matrix dimensions: ", nrow(genus_mat), " samples x ", ncol(genus_mat), " genera\n")
if (nrow(genus_mat) < 3) stop("Filtered sample count < 3, cannot perform RDA / PERMANOVA.")


meta_use <- meta_raw %>%
  filter(Sample_id %in% rownames(genus_mat)) %>%
  distinct(Sample_id, .keep_all = TRUE)

meta_use <- meta_use[match(rownames(genus_mat), meta_use$Sample_id), , drop = FALSE]

meta_use <- meta_use %>%
  mutate(
    category  = factor(category),
    continent = factor(continent),
    origin    = factor(origin)
  )

stopifnot(identical(rownames(genus_mat), meta_use$Sample_id))

cat("[INFO] Number of samples for RDA/PERMANOVA: ", nrow(meta_use), "\n")


cat("[5] Hellinger transformation of genus matrix...\n")
Y_hell <- decostand(genus_mat, method = "hellinger")


cat("[6] Performing RDA (Genus ~ origin + category + continent)...\n")

rda_genus <- rda(Y_hell ~ origin + category + continent, data = meta_use)
print(rda_genus)

eig <- rda_genus$CCA$eig
prop_explained <- eig / sum(eig)
prop_tbl <- tibble(
  axis = names(eig),
  eigenvalue = as.numeric(eig),
  prop = as.numeric(prop_explained),
  prop_pct = prop * 100
)

write.table(prop_tbl,
            file = file.path(out_dir, "RDA_genus_axis_eig.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("Explanation variance of top RDA constrained axes (%):\n")
print(round(prop_tbl$prop_pct[1:min(5, nrow(prop_tbl))], 2))


cat("[7] Performing overall ANOVA and by-terms ANOVA on RDA...\n")

anova_overall <- anova(rda_genus, permutations = N_PERM)
anova_terms   <- anova(rda_genus, permutations = N_PERM, by = "terms")

write.table(as.data.frame(anova_overall),
            file = file.path(out_dir, "RDA_genus_anova_overall.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)
write.table(as.data.frame(anova_terms),
            file = file.path(out_dir, "RDA_genus_anova_terms.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)

cat("overall ANOVA:\n")
print(anova_overall)
cat("by-terms ANOVA:\n")
print(anova_terms)


cat("[8] Extracting RDA1/RDA2 site & species scores and plotting...\n")

site_scores <- scores(rda_genus, display = "sites", choices = 1:2) %>%
  as.data.frame() %>%
  rownames_to_column("Sample_id") %>%
  left_join(meta_use, by = "Sample_id")

species_scores <- scores(rda_genus, display = "species", choices = 1:2) %>%
  as.data.frame() %>%
  rownames_to_column("Genus")

write.table(site_scores,
            file = file.path(out_dir, "RDA_genus_site_scores.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(species_scores,
            file = file.path(out_dir, "RDA_genus_species_scores.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

rda1_pct <- ifelse(length(prop_explained) >= 1, prop_explained[1] * 100, NA)
rda2_pct <- ifelse(length(prop_explained) >= 2, prop_explained[2] * 100, NA)

p_rda <- ggplot(site_scores, aes(x = RDA1, y = RDA2)) +
  geom_point(aes(color = origin, shape = category), alpha = 0.7, size = 1.8) +
  theme_bw(base_size = 12) +
  labs(
    x = paste0("RDA1 (", ifelse(is.na(rda1_pct), "NA", sprintf("%.2f", rda1_pct)), "%)"),
    y = paste0("RDA2 (", ifelse(is.na(rda2_pct), "NA", sprintf("%.2f", rda2_pct)), "%)"),
    title = "RDA on ARG-carrying bacteria (Genus, TPM_prop, Hellinger)"
  ) +
  theme(
    legend.position = "right",
    panel.grid = element_blank()
  )

pdf(file.path(out_dir, "RDA_genus_TPMprop_hellinger.pdf"),
    width = 7, height = 6)
print(p_rda)
dev.off()


cat("[9] PERMANOVA (adonis2): Comparing variance explained by origin / category / continent at genus level...\n")

D_hell <- dist(Y_hell, method = "euclidean")


per_multi <- adonis2(
  D_hell ~ origin + category + continent,
  data = meta_use,
  permutations = N_PERM,
  by = "margin"
)

per_multi_df <- as.data.frame(per_multi)
per_multi_df$term <- rownames(per_multi_df)


per_multi_terms <- per_multi_df %>%
  filter(!(term %in% c("Residual","Total"))) %>%
  transmute(
    model        = "multi",
    term         = term,
    Df           = Df,
    SumOfSqs     = SumOfSqs,
    R2           = R2,
    F            = F,
    p_value      = `Pr(>F)`,
    p_value_sci  = fmt_sci(`Pr(>F)`),
    n_samples    = attr(D_hell, "Size")
  )


simple_list <- list()
for (g in c("origin","category","continent")) {
  form <- as.formula(paste("D_hell ~", g))
  per_g <- adonis2(form, data = meta_use, permutations = N_PERM)
  df_g  <- as.data.frame(per_g)
  df_g$term <- rownames(df_g)
  
  df_g <- df_g[df_g$term == g, , drop = FALSE]
  if (nrow(df_g) == 0) next
  
  simple_list[[g]] <- tibble(
    model        = "single",
    term         = g,
    Df           = df_g$Df,
    SumOfSqs     = df_g$SumOfSqs,
    R2           = df_g$R2,
    F            = df_g$F,
    p_value      = df_g$`Pr(>F)`,
    p_value_sci  = fmt_sci(df_g$`Pr(>F)`),
    n_samples    = attr(D_hell, "Size")
  )
}

per_simple <- bind_rows(simple_list)

per_all <- bind_rows(per_multi_terms, per_simple)

write.table(per_all,
            file = file.path(out_dir, "PERMANOVA_genus_hellinger_all.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\n[PERMANOVA Results (Key Info)]\n")
print(
  per_all %>%
    arrange(model, term) %>%
    select(model, term, R2, F, p_value, p_value_sci, n_samples)
)

cat("\n[OK] Genus-level RDA + PERMANOVA complete! Output directory: ", normalizePath(out_dir), "\n")