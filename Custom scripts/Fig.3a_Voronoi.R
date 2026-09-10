library(WeightedTreemaps)
library(dplyr)
library(stringr)
library(grid)

setwd("/path/to/your/working/directory")

data_file <- "/path/to/your/data_file.tsv"

dat <- read.table(
  data_file,
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  quote = "",
  comment.char = ""
)

stopifnot(ncol(dat) == 13)

dat <- dat %>%
  mutate(
    Genus      = str_squish(as.character(Genus)),
    Species    = str_squish(as.character(Species)),
    Species_XY = str_squish(as.character(Species_XY)),
    Pathogen   = str_squish(as.character(Pathogen))
  ) %>%
  mutate(
    Species_final = case_when(
      !is.na(Species_XY) & Species_XY != "" ~ Species_XY,
      !is.na(Species)    & Species    != "" ~ Species,
      TRUE                                  ~ NA_character_
    )
  )

target_genera <- c(
  "Thermoactinomyces",
  "Lentibacillus",
  "Streptomyces",
  "Mammaliicoccus",
  "Kroppenstedtia",
  "Bacillus",
  "Enterobacter",
  "Klebsiella",
  "Acinetobacter",
  "Enterococcus",
  "Streptococcus",
  "Escherichia",
  "Staphylococcus"
)

dat_row <- dat %>%
  mutate(
    Genus_group = if_else(Genus %in% target_genera, Genus, "Others"),
    Species_block = if_else(
      Pathogen == "Unclassified",
      "Unclassified",
      Species_final
    )
  ) %>%
  filter(!is.na(Genus_group),
         !is.na(Species_block),
         Species_block != "")

df_tm <- dat_row %>%
  count(Genus_group, Species_block, name = "size")

df_genus_count <- df_tm %>%
  group_by(Genus_group) %>%
  summarise(n_genus = sum(size), .groups = "drop")

dat_plot <- df_tm %>%
  left_join(df_genus_count, by = "Genus_group") %>%
  mutate(
    Genus_lab   = paste0(Genus_group, "\n", n_genus),
    Species_lab = paste0(Species_block, "\n", size)
  )

levels_vec      <- c("Genus_lab", "Species_lab")
positioning_vec <- c("clustered_by_area", "clustered_by_area")

tm <- tryCatch(
  {
    voronoiTreemap(
      data         = dat_plot,
      levels       = levels_vec,
      cell_size    = "size",
      sort         = FALSE,
      filter       = 0,
      shape        = "circle",
      maxIteration = 200,
      error_tol    = 0.02,
      convergence  = "slow",
      seed         = 1,
      positioning  = positioning_vec,
      verbose      = TRUE
    )
  },
  error = function(e) {
    message("Circle shape failed, switching to rounded_rect: ", conditionMessage(e))
    voronoiTreemap(
      data         = dat_plot,
      levels       = levels_vec,
      cell_size    = "size",
      sort         = FALSE,
      filter       = 0,
      shape        = "rounded_rect",
      maxIteration = 200,
      error_tol    = 0.02,
      convergence  = "slow",
      seed         = 1,
      positioning  = positioning_vec,
      verbose      = TRUE
    )
  }
)

genus_levels <- unique(dat_plot$Genus_lab)
genus_names  <- sub("\n.*$", "", genus_levels)

genus_colors <- sapply(
  genus_names,
  function(g) {
    if (g == "Thermoactinomyces") "#37B2C8"
    else if (g == "Lentibacillus") "#5CC0CF"
    else if (g == "Streptomyces") "#80CED6"
    else if (g == "Mammaliicoccus") "#A5DCE0"
    else if (g == "Kroppenstedtia") "#C1E7E9"
    else if (g == "Bacillus") "#D0EEF1"
    else if (g == "Enterobacter") "#F8D6E0"
    else if (g == "Klebsiella") "#EBA2AF"
    else if (g == "Acinetobacter") "#DC8390"
    else if (g == "Enterococcus") "#C26571"
    else if (g == "Streptococcus") "#A84757"
    else if (g == "Escherichia") "#8D2A3D"
    else if (g == "Staphylococcus") "#710C23"
    else "#E9E9E9"
  }
)

cells_df  <- as.data.frame(tm)
poly_list <- get_polygons(tm)

centers <- t(vapply(poly_list, function(p) {
  coords <- p[[1]]
  poly_centroid(coords[, 1], coords[, 2])
}, numeric(2)))

centers <- as.data.frame(centers)
colnames(centers) <- c("x", "y")

cells_df$x <- centers$x
cells_df$y <- centers$y

species2genus <- dat_plot %>%
  distinct(Species_lab, Genus_group)

cells_df <- cells_df %>%
  left_join(species2genus, by = c("name" = "Species_lab"))

label_big <- subset(cells_df, level == 1)
label_big$label <- as.character(label_big$name)

label_small <- subset(cells_df, level == 2 & Genus_group != "Others")
if (nrow(label_small) > 0) {
  label_small$label <- as.character(label_small$name)
}

grid.newpage()
pushViewport(
  viewport(
    x = 0.5, y = 0.5,
    width  = 0.9,
    height = 0.9,
    xscale = c(0, 2000),
    yscale = c(0, 2000),
    name   = "treemap_vp"
  )
)

drawTreemap(
  tm,
  add          = TRUE,
  levels       = 1:2,
  color_type   = "categorical",
  color_level  = 1,
  color_palette = unname(genus_colors),
  border_level = c(1, 2),
  border_size  = 6,
  border_color = "white",
  label_level  = NULL,
  legend       = FALSE
)

if (nrow(label_big) > 0) {
  grid.text(
    label_big$label,
    x = label_big$x,
    y = label_big$y,
    default.units = "native",
    gp = gpar(
      col        = "black",
      fontsize   = 7,
      fontface   = "bold",
      lineheight = 0.9
    )
  )
}

if (nrow(label_small) > 0) {
  grid.text(
    label_small$label,
    x = label_small$x,
    y = label_small$y,
    default.units = "native",
    gp = gpar(
      col        = "black",
      fontsize   = 5,
      fontface   = "plain",
      lineheight = 0.9
    )
  )
}

upViewport()