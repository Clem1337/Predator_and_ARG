library(tidyverse)
library(readxl)
library(janitor)
library(pheatmap)
library(stringr)

setwd("/path/to/your/working/directory")

cat("Reading abundance table...\n")
abund_raw <- read_xlsx("/path/to/abundance_table.xlsx") %>%
  clean_names()  

sample_cols <- setdiff(names(abund_raw), c("mode", "genus"))
cat("Successfully identified", length(sample_cols), "sample columns.\n")


cat("Reading metadata...\n")
meta_path <- "/path/to/metadata.xlsx"

meta <- read_xlsx(meta_path) %>%
  clean_names()

meta_small <- meta %>%
  mutate(
    sample_id = str_trim(sample_id),
    origin = str_trim(origin)
  ) %>%
  select(sample_id, category, continent, origin)

cat("Data reading and environment initialization complete!\n")

animal_samples <- meta_small$sample_id[meta_small$origin == "animal" & meta_small$sample_id %in% colnames(mat_raw)]
plant_samples  <- meta_small$sample_id[meta_small$origin == "plant" & meta_small$sample_id %in% colnames(mat_raw)]

total_animal <- length(animal_samples)
total_plant  <- length(plant_samples)

summary_df <- data.frame(genus = rownames(mat_raw)) %>%
  mutate(
    n_animal_present = rowSums(mat_raw[, animal_samples, drop = FALSE] > 0, na.rm = TRUE),
    n_plant_present  = rowSums(mat_raw[, plant_samples, drop = FALSE] > 0, na.rm = TRUE),
    
    n_animal_absent  = pmax(total_animal - n_animal_present, 0),
    n_plant_absent   = pmax(total_plant - n_plant_present, 0)
  )

df1 <- summary_df %>%
  select(genus,
         animal       = n_animal_present,
         total_animal = n_animal_absent,
         plant        = n_plant_present,
         total_plant  = n_plant_absent) %>%
  pivot_longer(cols = -genus, names_to = "slice", values_to = "value") %>%
  mutate(panel = "All Samples")

df2 <- summary_df %>%
  select(genus,
         animal       = n_animal_present,
         total_animal = n_animal_absent) %>%
  pivot_longer(cols = -genus, names_to = "slice", values_to = "value") %>%
  mutate(panel = "All Animal Samples")

df3 <- summary_df %>%
  select(genus,
         plant       = n_plant_present,
         total_plant = n_plant_absent) %>%
  pivot_longer(cols = -genus, names_to = "slice", values_to = "value") %>%
  mutate(panel = "All Plant Samples")

df4 <- summary_df %>%
  select(genus,
         animal = n_animal_present,
         plant  = n_plant_present) %>%
  pivot_longer(cols = -genus, names_to = "slice", values_to = "value") %>%
  mutate(panel = "Source Proportion in Detected Samples")


pie_df <- bind_rows(df1, df2, df3, df4)

genus_levels <- rev(target_taxa)

pie_df <- pie_df %>%
  mutate(
    genus = factor(genus, levels = genus_levels),
    panel = factor(panel, levels = c("All Samples", "All Animal Samples", "All Plant Samples", "Source Proportion in Detected Samples")),
    slice = factor(slice, levels = c("animal", "plant", "total_animal", "total_plant"))
  )

fill_colors <- c(
  animal       = "#C2842E",  
  plant        = "#89D0C4",  
  total_animal = "#F2E2CD",  
  total_plant  = "#D2ECE8"   
)

pie_plot <- ggplot(pie_df, aes(x = "", y = value, fill = slice)) +
  geom_col(width = 1, color = NA, position = "fill") + 
  coord_polar(theta = "y") +
  facet_grid(genus ~ panel) +
  scale_fill_manual(values = fill_colors, drop = FALSE) +
  theme_void() +
  theme(
    strip.text.x    = element_text(size = 10, face = "bold", margin = margin(b = 5)),
    strip.text.y    = element_text(size = 9, angle = 0, hjust = 0, face = "italic"), 
    legend.position = "bottom",
    legend.title    = element_blank(),
    panel.spacing   = unit(0.2, "lines")
  )

print(pie_plot)