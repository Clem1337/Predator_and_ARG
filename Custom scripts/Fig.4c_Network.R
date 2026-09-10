setwd("/path/to/your/working/directory")

library(readxl)
library(dplyr)
library(tidygraph)
library(ggraph)
library(ggplot2)
library(ggrepel)
library(grid)

IN_XLSX <- "your_data_file.xlsx"
MODES   <- c("predator", "symbiotroph", "phototroph", "unknown")

read_clean_sheet <- function(path, sheet_name) {
  if (!file.exists(path)) stop(paste("Error: File not found", path))
  df <- read_xlsx(path, sheet = sheet_name)
  
  required_cols <- c("source", "target", "sign")
  miss <- setdiff(required_cols, names(df))
  if (length(miss) > 0) stop(paste0("Error: sheet '", sheet_name, "' is missing columns: ", paste(miss, collapse = ", ")))
  
  df %>%
    mutate(
      source = trimws(as.character(source)),
      target = trimws(as.character(target)),
      sign   = tolower(trimws(as.character(sign)))
    ) %>%
    filter(!is.na(source) & !is.na(target)) %>%
    filter(source != "" & target != "") %>%
    filter(source != "NA" & target != "NA") %>%
    mutate(
      sign = if_else(sign == "positive", "positive", "negative")
    )
}

mode_dfs <- lapply(MODES, function(m) read_clean_sheet(IN_XLSX, m))
names(mode_dfs) <- MODES

edges_all <- bind_rows(lapply(MODES, function(m) {
  df <- mode_dfs[[m]]
  if (nrow(df) == 0) return(df[0, , drop = FALSE])
  
  df %>%
    mutate(
      mode = m,
      from = paste0(m, "::", source),
      to   = target
    ) %>%
    select(from, to, sign, mode, everything())
}))

edges_all <- edges_all %>%
  filter(!is.na(from) & !is.na(to) & from != "" & to != "")

if (nrow(edges_all) == 0) stop("Error: Edge table is empty. Please check the Excel content.")

carrier_names <- sort(unique(edges_all$to))
carrier_df <- data.frame(name = carrier_names, stringsAsFactors = FALSE) %>%
  mutate(
    group = "Carrier",
    node_type = "Carrier",
    label = name
  )

source_names <- sort(unique(edges_all$from))
source_df <- data.frame(name = source_names, stringsAsFactors = FALSE) %>%
  mutate(
    group = sub("::.*$", "", name),
    node_type = "Protist",
    label = sub("^.*::", "", name)
  )

nodes_all <- bind_rows(source_df, carrier_df)

nodes_all$group <- factor(
  nodes_all$group,
  levels = c("predator", "symbiotroph", "phototroph", "unknown", "Carrier")
)
nodes_all$node_type <- factor(nodes_all$node_type, levels = c("Protist", "Carrier"))

protist_sign_stat <- edges_all %>%
  group_by(from) %>%
  summarise(n_sign = n_distinct(sign), .groups = "drop") %>%
  mutate(show_protist_label = (n_sign == 1))

nodes_all <- nodes_all %>%
  left_join(protist_sign_stat, by = c("name" = "from")) %>%
  mutate(
    label_carrier = if_else(node_type == "Carrier", label, NA_character_),
    label_protist = if_else(node_type == "Protist" & coalesce(show_protist_label, FALSE), label, NA_character_),
    protist_filled = if_else(node_type == "Protist" & !is.na(label_protist), TRUE, FALSE)
  )

set.seed(123)

fill_circle <- function(n, center, r) {
  if (n == 0) return(data.frame(x = numeric(0), y = numeric(0)))
  theta <- seq(0, 2*pi, length.out = n + 1)[-(n + 1)]
  data.frame(
    x = center[1] + r * cos(theta),
    y = center[2] + r * sin(theta)
  )
}

nodes_all$x <- 0
nodes_all$y <- 0

centers <- list(
  predator     = c(-4.5,  3),
  symbiotroph  = c(-4.5, -3),
  phototroph   = c( 4.5,  3),
  unknown      = c( 4.5, -3),
  Carrier      = c( 0,  0)
)

radii <- list(
  predator     = 1.6,
  symbiotroph  = 1.6,
  phototroph   = 1.6,
  unknown      = 1.6,
  Carrier      = 1.5
)

for (g in levels(nodes_all$group)) {
  idx <- which(nodes_all$group == g)
  if (length(idx) > 0) {
    xy <- fill_circle(length(idx), centers[[as.character(g)]], radii[[as.character(g)]])
    nodes_all[idx, c("x", "y")] <- xy
  }
}

nodes_all <- nodes_all %>%
  mutate(
    shape_plot = case_when(
      node_type == "Carrier" ~ "Carrier",
      node_type == "Protist" & protist_filled ~ "Protist_filled",
      node_type == "Protist" & !protist_filled ~ "Protist_hollow",
      TRUE ~ "Other"
    )
  )

shape_vals <- c(
  Carrier        = 15,  
  Protist_filled = 16,  
  Protist_hollow = 1,   
  Other          = 3
)

graph_net <- tbl_graph(
  nodes = nodes_all,
  edges = edges_all %>% select(from, to, sign, mode),
  directed = FALSE,
  node_key = "name"
)

p_net <- ggraph(graph_net, layout = "manual", x = x, y = y) +
  geom_edge_link(aes(color = sign), alpha = 0.4, width = 0.7) +
  scale_edge_colour_manual(values = c(positive = "#6F81BF", negative = "#FFDA76")) +
  geom_node_point(aes(shape = shape_plot), size = 3, stroke = 0.8, colour = "black") +
  scale_shape_manual(values = shape_vals) +
  geom_node_text(
    aes(label = label_carrier),
    na.rm = TRUE,
    size = 2.8,
    repel = TRUE,
    point.padding = unit(0.12, "lines")
  ) +
  geom_node_text(
    aes(label = label_protist),
    na.rm = TRUE,
    size = 2.8,
    repel = TRUE,
    point.padding = unit(0.12, "lines")
  ) +
  theme_void() +
  coord_equal()   

print(p_net)

library(dplyr)
library(tidyr)

target_carriers <- NULL

carriers_use <- if (is.null(target_carriers)) {
  sort(unique(edges_all$to))
} else {
  target_carriers
}

carrier_neg_pos <- edges_all %>%
  filter(to %in% carriers_use) %>%
  mutate(sign = if_else(sign == "positive", "positive", "negative")) %>%
  count(to, sign, name = "n_edges") %>%
  tidyr::complete(to = carriers_use, sign = c("negative","positive"), fill = list(n_edges = 0)) %>%
  pivot_wider(names_from = sign, values_from = n_edges) %>%
  mutate(
    total = negative + positive,
    negative_prop = if_else(total > 0, negative / total, NA_real_),
    positive_prop = if_else(total > 0, positive / total, NA_real_),
    neg_to_pos_ratio = if_else(positive > 0, negative / positive, NA_real_)  
  ) %>%
  arrange(desc(negative_prop), desc(total))

print(carrier_neg_pos)

write.csv(carrier_neg_pos, "13carriers_negative_positive_ratio.csv",
          row.names = FALSE, fileEncoding = "UTF-8")