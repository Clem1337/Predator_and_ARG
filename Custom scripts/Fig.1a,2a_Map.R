setwd("/path/to/your/working/directory")

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(stringr)
  library(sf)
  library(rnaturalearth)
  library(ggplot2)
})

xlsx_path <- "your_data_file.xlsx"
dat_iso <- readxl::read_excel(xlsx_path, sheet = "country") %>%
  rename(iso3 = country,
         arg_avg = avg_ARG_copies_per_cell) %>%
  mutate(
    iso3    = str_to_upper(str_trim(as.character(iso3))),
    arg_avg = suppressWarnings(as.numeric(arg_avg))
  ) %>%
  filter(!is.na(iso3), !is.na(arg_avg)) %>%
  group_by(iso3) %>%                      
  summarise(arg_avg = mean(arg_avg, na.rm = TRUE), .groups = "drop")

# dat_iso <- dat_iso %>% mutate(iso3 = dplyr::recode(iso3, "ALA" = "FIN", .default = iso3))

world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") %>%
  dplyr::filter(continent != "Antarctica")
world_robin <- sf::st_transform(world, "+proj=robin")

world_dat <- world_robin %>%
  left_join(dat_iso, by = c("iso_a3" = "iso3"))

p <- ggplot(world_dat) +
  geom_sf(aes(fill = arg_avg), color = "white", linewidth = 0.24) +
  scale_fill_gradient(
    name = NULL,
    low = "#F7D48E", high = "#903A27",
    na.value = "grey85"
  ) +
  theme_void(base_size = 12) +
  theme(
    legend.position = c(0.12, 0.32),        
    legend.key.height = unit(4, "mm"),
    legend.key.width  = unit(4, "mm")
  )

print(p)